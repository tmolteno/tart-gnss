// Copyright (c) 2026 Tim Molteno <tim@elec.ac.nz>
// SPDX-License-Identifier: GPL-3.0

//! Galileo E1 Open Service pilot (E1-C) signal acquisition via FFT-based
//! circular cross-correlation.
//!
//! The E1-C pilot channel uses 4092-chip primary memory codes (4 ms period
//! at 1.023 Mcps), with a 25-chip secondary code overlay.

#[path = "galileo_codes.rs"]
mod galileo_codes;

use crate::correlate::correlate_code;
use crate::observation::Observation;
use galileo_codes::{GALILEO_E1_C_CODES, GALILEO_E1_CHIPS, GALILEO_E1_CODE_PERIOD};
pub use galileo_codes::GALILEO_E1_NUM_SATS;
use num_complex::Complex;
use rayon::prelude::*;
use rustfft::FftPlanner;
use serde::Serialize;
use std::f64::consts::TAU;
use std::sync::atomic::{AtomicUsize, Ordering};

/// Result of a Galileo E1 signal acquisition attempt for a single satellite.
#[derive(Debug, Clone, Serialize)]
pub struct GalileoAcquisitionResult {
    pub prn: usize,
    /// Peak correlation magnitude (normalised).
    pub signal_strength: f64,
    /// Second-highest correlation magnitude at the same frequency,
    /// at least 1 chip away from the peak (for ACR C/N0 estimation).
    pub second_peak: f64,
    /// Code-phase offset in fractions of a code period [0, 1).
    pub codephase_frac: f64,
    /// Doppler frequency offset in Hz (relative to centre frequency).
    pub frequency: f64,
}

/// Per-SV Galileo acquisition result with per-antenna measurements.
#[derive(Debug, Clone, Serialize)]
pub struct GalileoPrnResult {
    /// Constellation label, e.g. "GSAT03".
    pub sv: String,
    /// Per-antenna signal strengths.
    pub strengths: Vec<f64>,
    /// Per-antenna code-phase offsets (fraction of a code period).
    pub phases: Vec<f64>,
    /// Per-antenna Doppler frequency offsets (Hz).
    pub freqs: Vec<f64>,
    /// ACR C/N0 estimate per antenna (dB-Hz), if available.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cn0_acr: Option<Vec<f64>>,
    /// Median phase across antennas (only when >1 antenna).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub phase_median: Option<f64>,
    /// MAD of phase across antennas (only when >1 antenna).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub phase_mad: Option<f64>,
    /// Median frequency across antennas (only when >1 antenna).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub freq_median: Option<f64>,
    /// MAD of frequency across antennas (only when >1 antenna).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub freq_mad: Option<f64>,
}

/// Collection of acquisition results for all Galileo SVs, grouped by PRN.
#[derive(Debug, Clone, Serialize)]
pub struct GalileoAllAcquisitionOutput {
    pub antenna_numbers: Vec<usize>,
    pub results: Vec<GalileoPrnResult>,
}

// ---------------------------------------------------------------------------
// Code generation
// ---------------------------------------------------------------------------

/// Decode a hex string (MSB-first per hex digit) into ±1 bipolar chips.
///
/// Each hex character encodes 4 chips.  `hex_char` 0 → chips [-1,-1,-1,-1],
/// F → [+1,+1,+1,+1].
fn hex_to_chips(hex: &str) -> [f64; GALILEO_E1_CHIPS] {
    let mut chips = [0.0f64; GALILEO_E1_CHIPS];
    for (i, ch) in hex.chars().enumerate() {
        let val = ch.to_digit(16).expect("invalid hex char in Galileo code");
        for bit in 0..4 {
            let chip_idx = i * 4 + bit;
            if chip_idx < GALILEO_E1_CHIPS {
                chips[chip_idx] = if (val >> (3 - bit)) & 1 == 1 { 1.0 } else { -1.0 };
            }
        }
    }
    chips
}

/// Generate the 4092-chip Galileo E1-C primary code for a given PRN (1-based).
///
/// Returns a `[f64; GALILEO_E1_CHIPS]` array of ±1 values.
pub fn generate_e1c_code(prn: usize) -> [f64; GALILEO_E1_CHIPS] {
    assert!(
        prn >= 1 && prn <= GALILEO_E1_NUM_SATS,
        "Galileo PRN must be in 1..={GALILEO_E1_NUM_SATS}, got {prn}"
    );
    hex_to_chips(GALILEO_E1_C_CODES[prn - 1])
}

/// Resample the E1-C code to `samples_per_code` samples per 4-ms period,
/// repeating for `epochs` full periods.
pub fn e1c_code_resampled(samples_per_code: f64, prn: usize, epochs: f64) -> Vec<f64> {
    let ca = generate_e1c_code(prn);
    let samples_per_chip = samples_per_code / GALILEO_E1_CHIPS as f64;
    let num_samples = (samples_per_code * epochs).floor() as usize;

    (0..num_samples)
        .map(|n| {
            let idx = ((n as f64 / samples_per_chip).floor() as usize) % GALILEO_E1_CHIPS;
            ca[idx]
        })
        .collect()
}

// ---------------------------------------------------------------------------
// FFT-based acquisition (single PRN)
// ---------------------------------------------------------------------------

/// Perform parallel code-phase search (FFT circular cross-correlation) for a
/// single Galileo E1-C PRN over a frequency search band.
///
/// `signal_f32` and `phasepoints` are pre-computed by the caller to avoid
/// redundant allocations across PRN/antenna iterations.
pub fn acquire_galileo_single(
    signal_f32: &[f32],
    phasepoints: &[f64],
    sampling_freq: f64,
    center_freq: f64,
    search_band: f64,
    prn: usize,
    samples_per_code_period: usize,
) -> GalileoAcquisitionResult {
    let num_samples = signal_f32.len();

    // Resample the code to exactly `num_samples` (including a partial final
    // period), so its length matches the FFT size. Flooring to whole periods
    // here made the code shorter than the FFT when the signal length is not a
    // multiple of the code period, triggering a rustfft buffer-too-small panic.
    let epochs_available = num_samples as f64 / samples_per_code_period as f64;

    // --- Frequency bins ----------------------------------------------------
    let freq_bin_size: f64 = 300.0;
    let n_freq_bins = (2.0 * search_band / freq_bin_size).round() as usize + 1;
    let fc: Vec<f64> = (0..n_freq_bins)
        .map(|i| {
            center_freq - search_band
                + (i as f64) * (2.0 * search_band) / (n_freq_bins as f64 - 1.0)
        })
        .collect();

    // --- FFT planner (thread-local) ----------------------------------------
    thread_local! {
        static PLANNER: std::cell::RefCell<FftPlanner<f32>> =
            std::cell::RefCell::new(FftPlanner::new());
    }

    let (fft, ifft) = PLANNER.with(|p| {
        let mut planner = p.borrow_mut();
        (
            planner.plan_fft_forward(num_samples),
            planner.plan_fft_inverse(num_samples),
        )
    });

    // --- Local code replica (FFT + conjugate once) -------------------------
    let code = e1c_code_resampled(samples_per_code_period as f64, prn, epochs_available);
    let mut code_complex: Vec<Complex<f32>> = code
        .iter()
        .map(|&v| Complex::new(v as f32, 0.0))
        .collect();
    fft.process(&mut code_complex);
    for c in &mut code_complex {
        *c = c.conj();
    }

    // --- Parallel frequency-bin correlation --------------------------------
    let peak = correlate_code(
        signal_f32,
        phasepoints,
        &code_complex,
        &fc,
        num_samples,
        samples_per_code_period,
        sampling_freq,
        fft,
        ifft,
    );

    let codephase_in_samples = peak.best_codephase % samples_per_code_period;
    let codephase_frac = codephase_in_samples as f64 / samples_per_code_period as f64;
    let frequency = fc[peak.best_freq_idx] - center_freq;

    GalileoAcquisitionResult {
        prn,
        signal_strength: peak.best_peak as f64,
        second_peak: peak.second_peak as f64,
        codephase_frac,
        frequency,
    }
}

// ---------------------------------------------------------------------------
// All-SV search
// ---------------------------------------------------------------------------

/// Search for all 50 Galileo E1-C SVs across selected antennas.
///
/// If `ant_filter` is `Some(antennas)`, only those antennas are used; otherwise all
/// antennas are searched.  PRN processing is parallelised via rayon.
///
/// For each PRN, per-antenna signal strengths, code-phase offsets, and
/// Doppler frequency offsets are collected.
pub fn acquire_all_galileo(
    obs: &Observation,
    center_freq: f64,
    search_band: f64,
    ant_filter: Option<Vec<usize>>,
    prn_filter: Option<&[usize]>,
    debug: bool,
    cn0: bool,
) -> GalileoAllAcquisitionOutput {
    let sampling_freq = obs.get_sampling_rate();
    let n_ant = obs.config.num_antenna();
    let num_samples_per_code_period = (sampling_freq * GALILEO_E1_CODE_PERIOD) as usize;
    // Use 8 ms of data (2 code periods) per antenna
    let num_samples = 2 * num_samples_per_code_period;

    // Which antennas to search
    let ant_indices: Vec<usize> = ant_filter.unwrap_or_else(|| (0..n_ant).collect());

    // Pre-extract and de-mean all antenna data, convert to f32 once
    let ant_data: Vec<(Vec<f32>, Vec<f64>)> = ant_indices
        .iter()
        .map(|&ant_idx| {
            let bipolar = obs.get_antenna(ant_idx);
            let mean = bipolar.iter().sum::<f64>() / bipolar.len() as f64;
            let raw: Vec<f64> = bipolar.iter().map(|&v| v - mean).collect();
            let len = num_samples.min(raw.len());
            let signal_f32: Vec<f32> = raw[..len].iter().map(|&v| v as f32).collect();
            let phase_const = TAU / sampling_freq;
            let phasepoints: Vec<f64> =
                (0..len).map(|i| phase_const * i as f64).collect();
            (signal_f32, phasepoints)
        })
        .collect();

    // Parallel PRN search ----------------------------------------------------
    let prn_list: Vec<usize> = if let Some(filter) = prn_filter {
        filter.to_vec()
    } else {
        (1..=GALILEO_E1_NUM_SATS).collect()
    };
    let total = prn_list.len();
    let counter = AtomicUsize::new(0);

    let mut results: Vec<GalileoPrnResult> = prn_list
        .into_par_iter()
        .map(|prn| {
            let mut strengths = Vec::with_capacity(ant_indices.len());
            let mut second_peaks: Vec<f64> = if cn0 {
                Vec::with_capacity(ant_indices.len())
            } else {
                Vec::new()
            };
            let mut phases = Vec::with_capacity(ant_indices.len());
            let mut freqs = Vec::with_capacity(ant_indices.len());

            for (i, (signal_f32, phasepoints)) in ant_data.iter().enumerate() {
                let result = acquire_galileo_single(
                    signal_f32,
                    phasepoints,
                    sampling_freq,
                    center_freq,
                    search_band,
                    prn,
                    num_samples_per_code_period,
                );

                if debug {
                    eprintln!(
                        "  galileo PRN {:2} ant {:2}: strength={:.3}  phase={:.6}  freq={:.1} Hz",
                        prn,
                        ant_indices[i],
                        result.signal_strength,
                        result.codephase_frac,
                        result.frequency
                    );
                }

                strengths.push(result.signal_strength);
                if cn0 {
                    second_peaks.push(result.second_peak);
                }
                phases.push(result.codephase_frac);
                freqs.push(result.frequency);
            }

            // Compute ACR C/N0 per antenna
            let cn0_acr: Option<Vec<f64>> = if cn0 {
                let cn0s: Vec<f64> = strengths
                    .iter()
                    .zip(second_peaks.iter())
                    .filter_map(|(&v_m, &v_s)| {
                        if v_s > 0.0 && v_m > v_s {
                            let r_a = (v_m / v_s).powi(2);
                            crate::acr::estimate_cn0(r_a)
                        } else {
                            None
                        }
                    })
                    .collect();
                if cn0s.is_empty() { None } else { Some(cn0s) }
            } else {
                None
            };

            let n = counter.fetch_add(1, Ordering::Relaxed) + 1;
            eprintln!("  galileo [{n}/{total}]");
            if debug {
                eprintln!("  galileo [{n}/{total}] PRN {prn:02} complete");
            }

            let (phase_median, phase_mad, freq_median, freq_mad) =
                if phases.len() > 1 {
                    let pm = crate::stats::median(&phases);
                    let fm = crate::stats::median(&freqs);
                    (
                        Some(pm),
                        Some(crate::stats::mad(&phases, pm)),
                        Some(fm),
                        Some(crate::stats::mad(&freqs, fm)),
                    )
                } else {
                    (None, None, None, None)
                };

            GalileoPrnResult {
                sv: format!("GSAT{prn:02}"),
                strengths,
                phases,
                freqs,
                cn0_acr,
                phase_median,
                phase_mad,
                freq_median,
                freq_mad,
            }
        })
        .collect();

    results.sort_by(|a, b| a.sv.cmp(&b.sv));

    GalileoAllAcquisitionOutput { antenna_numbers: ant_indices, results }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::{expected_recovered_sample, phasepoints, synth_signal};

    #[test]
    fn test_generate_e1c_code_length() {
        let code = generate_e1c_code(1);
        assert_eq!(code.len(), GALILEO_E1_CHIPS);
    }

    #[test]
    fn test_generate_e1c_code_bipolar() {
        let code = generate_e1c_code(5);
        for &v in &code {
            assert!(v == 1.0 || v == -1.0, "unexpected value {v}");
        }
    }

    #[test]
    fn test_e1c_code_resampled_length() {
        let samples_per_code = 4092.0;
        let code = e1c_code_resampled(samples_per_code, 1, 2.0);
        assert_eq!(code.len(), 8184);
    }

    #[test]
    fn test_all_prns_generate() {
        for prn in 1..=GALILEO_E1_NUM_SATS {
            let code = generate_e1c_code(prn);
            assert_eq!(code.len(), GALILEO_E1_CHIPS, "PRN {prn} wrong length");
        }
    }

    #[test]
    fn test_acquire_galileo_single_recovers_delay() {
        let fs = 1.023e6; // one sample per chip; 4092 chips per 4 ms period
        let period = GALILEO_E1_CHIPS;
        let code = e1c_code_resampled(period as f64, 1, 2.0); // 2 epochs

        let delay = 900usize;
        let sig = synth_signal(&code, period, delay, 0.0, fs, 0.05, 11);
        let r = acquire_galileo_single(
            &sig, &phasepoints(fs, sig.len()), fs, 0.0, 3000.0, 1, period,
        );
        let recovered = (r.codephase_frac * period as f64).round() as usize;
        assert_eq!(recovered, expected_recovered_sample(period, delay));
        assert!(r.codephase_frac >= 0.0 && r.codephase_frac < 1.0);
        assert!(r.signal_strength > 0.0);
    }

    #[test]
    fn test_acquire_galileo_single_noise_contrast() {
        let fs = 1.023e6;
        let period = GALILEO_E1_CHIPS;
        let code = e1c_code_resampled(period as f64, 1, 2.0);
        let n = code.len();

        let injected = synth_signal(&code, period, 600, 0.0, fs, 0.05, 12);
        let r_inj = acquire_galileo_single(
            &injected, &phasepoints(fs, n), fs, 0.0, 3000.0, 1, period,
        );
        let pure_noise: Vec<f32> = (0..n).map(|i| ((i as f64 * 2.7).sin() * 0.05) as f32).collect();
        let r_noise = acquire_galileo_single(
            &pure_noise, &phasepoints(fs, n), fs, 0.0, 3000.0, 1, period,
        );
        assert!(r_inj.signal_strength > 10.0 * r_noise.signal_strength);
    }

    #[test]
    fn test_acquire_non_period_multiple_no_panic() {
        // Signal length is NOT a whole multiple of the code period (2.5
        // periods). Regression: rustfft used to panic with "Provided FFT
        // buffer was too small" because the local code replica was resampled
        // to fewer samples than the FFT size.
        let period = GALILEO_E1_CHIPS;
        let fs = 1.023e6;
        let n = (2.5 * period as f64) as usize;
        let code = e1c_code_resampled(period as f64, 1, 2.5); // length n, partial period
        let delay = 1000usize;
        let sig: Vec<f32> = (0..n).map(|i| code[(i + delay) % period] as f32).collect();
        let res = acquire_galileo_single(
            &sig, &phasepoints(fs, n), fs, 0.0, 3000.0, 1, period,
        );
        assert!(res.codephase_frac >= 0.0 && res.codephase_frac < 1.0);
        assert!(res.signal_strength > 0.0);
    }
}
