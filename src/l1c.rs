// Copyright (c) 2026 Tim Molteno <tim@elec.ac.nz>
// SPDX-License-Identifier: GPL-3.0

//! GPS L1C pilot (L1Cp) signal acquisition via FFT-based
//! circular cross-correlation.
//!
//! The L1C pilot channel uses 10230-chip Weil primary codes (10 ms period
//! at 1.023 Mcps), with BOC(1,1) modulation applied for acquisition.
//! Reference: IS-GPS-800H.

#[path = "l1c_codes.rs"]
mod l1c_codes;

// Re-export for use by main.rs
pub use l1c_codes::{L1C_NUM_SATS, l1c_code_resampled};

use crate::correlate::correlate_code;
use crate::observation::Observation;
use l1c_codes::L1C_CODE_PERIOD;
use rayon::prelude::*;
use rustfft::FftPlanner;
use serde::Serialize;
use std::f64::consts::TAU;
use std::sync::atomic::{AtomicUsize, Ordering};

/// IF frequency for GPS L1 band (Hz).
pub const L1C_IF: f64 = 4.092e6;

/// Doppler search bandwidth (Hz).
pub const L1C_SEARCH_BAND: f64 = 6000.0;

/// Result of a GPS L1C signal acquisition attempt for a single satellite.
#[derive(Debug, Clone, Serialize)]
pub struct L1CAcquisitionResult {
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

/// Per-SV GPS L1C acquisition result with per-antenna measurements.
#[derive(Debug, Clone, Serialize)]
pub struct L1CPrnResult {
    /// Constellation label, e.g. "GPSL1C03".
    pub sv: String,
    /// Per-antenna signal strengths.
    pub strengths: Vec<f64>,
    /// Per-antenna code-phase offsets (fraction of a code period).
    pub phases: Vec<f64>,
    /// Per-antenna Doppler frequency offsets (Hz).
    pub freqs: Vec<f64>,
    /// ACR C/N0 estimate per antenna (dB-Hz), aligned with the per-antenna
    /// vectors above; `null` where the estimate failed (e.g. a dead channel).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cn0_acr: Option<Vec<Option<f64>>>,
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

/// Collection of acquisition results for all GPS L1C SVs, grouped by PRN.
#[derive(Debug, Clone, Serialize)]
pub struct L1CAllAcquisitionOutput {
    pub antenna_numbers: Vec<usize>,
    pub results: Vec<L1CPrnResult>,
}

// ---------------------------------------------------------------------------
// FFT-based acquisition (single PRN)
// ---------------------------------------------------------------------------

/// Perform parallel code-phase search (FFT circular cross-correlation) for a
/// single GPS L1C PRN over a frequency search band.
///
/// Uses the BOC(1,1)-modulated pilot code as the local replica.
///
/// `signal_f32` and `phasepoints` are pre-computed by the caller to avoid
/// redundant allocations across PRN/antenna iterations.
pub fn acquire_l1c_single(
    signal_f32: &[f32],
    phasepoints: &[f64],
    sampling_freq: f64,
    center_freq: f64,
    search_band: f64,
    prn: usize,
    samples_per_code_period: usize,
) -> L1CAcquisitionResult {
    let num_samples = signal_f32.len();

    // --- Frequency bins ----------------------------------------------------
    let freq_bin_size: f64 = 300.0;
    let n_freq_bins = (2.0 * search_band / freq_bin_size).round() as usize + 1;
    let fc: Vec<f64> = if n_freq_bins == 1 {
        vec![center_freq]
    } else {
        (0..n_freq_bins)
            .map(|i| {
                center_freq - search_band
                    + (i as f64) * (2.0 * search_band) / (n_freq_bins as f64 - 1.0)
            })
            .collect()
    };

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

    // --- Local code replicas (FFT + conjugate once each) -------------------
    // Resample the code to exactly `num_samples` so its length matches the
    // FFT size (a partial final period is allowed).  Two hypotheses cover a
    // possible sign flip between the two integrated periods (1800-chip
    // overlay code): `[c, c]` and `[c, -c]`.
    let code = l1c_code_resampled(samples_per_code_period as f64, prn, num_samples);
    let (code_complex, code_complex_alt) =
        crate::correlate::code_spectra(&code, samples_per_code_period, &fft);

    // --- Parallel frequency-bin correlation --------------------------------
    let peak = correlate_code(
        signal_f32,
        phasepoints,
        &code_complex,
        &code_complex_alt,
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

    L1CAcquisitionResult {
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

/// Search for all 63 GPS L1C SVs across selected antennas.
///
/// If `ant_filter` is `Some(antennas)`, only those antennas are used; otherwise all
/// antennas are searched.  PRN processing is parallelised via rayon.
///
/// For each PRN, per-antenna signal strengths, code-phase offsets, and
/// Doppler frequency offsets are collected.
pub fn acquire_all_l1c(
    obs: &Observation,
    center_freq: f64,
    search_band: f64,
    ant_filter: Option<Vec<usize>>,
    prn_filter: Option<&[usize]>,
    debug: bool,
    cn0: bool,
) -> L1CAllAcquisitionOutput {
    let sampling_freq = obs.get_sampling_rate();
    let n_ant = obs.config.num_antenna();
    let num_samples_per_code_period = (sampling_freq * L1C_CODE_PERIOD) as usize;
    // Use 20 ms of data (2 code periods) per antenna
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
        (1..=L1C_NUM_SATS).collect()
    };
    let total = prn_list.len();
    let counter = AtomicUsize::new(0);

    let mut results: Vec<L1CPrnResult> = prn_list
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
                let result =
                    acquire_l1c_single(signal_f32, phasepoints, sampling_freq, center_freq, search_band, prn, num_samples_per_code_period);

                if debug {
                    eprintln!(
                        "  l1c PRN {:2} ant {:2}: strength={:.3}  phase={:.6}  freq={:.1} Hz",
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

            // Compute ACR C/N0 per antenna, keeping one entry per antenna
            // (None where the peak ratio is unusable) so the vector stays
            // aligned with the other per-antenna results.
            let cn0_acr: Option<Vec<Option<f64>>> = if cn0 {
                let cn0s: Vec<Option<f64>> = strengths
                    .iter()
                    .zip(second_peaks.iter())
                    .map(|(&v_m, &v_s)| {
                        if v_s > 0.0 && v_m > v_s {
                            let r_a = (v_m / v_s).powi(2);
                            crate::acr::estimate_cn0(r_a, crate::acr::GPS_L1C_ACR_TABLE)
                        } else {
                            None
                        }
                    })
                    .collect();
                if cn0s.iter().all(Option::is_none) { None } else { Some(cn0s) }
            } else {
                None
            };

            let n = counter.fetch_add(1, Ordering::Relaxed) + 1;
            eprintln!("  l1c [{n}/{total}]");
            if debug {
                eprintln!("  l1c [{n}/{total}] PRN {prn:02} complete");
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

            L1CPrnResult {
                sv: format!("GPSL1C{prn:02}"),
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

    // Sort by SV label to maintain ordering
    results.sort_by(|a, b| a.sv.cmp(&b.sv));

    L1CAllAcquisitionOutput { antenna_numbers: ant_indices, results }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::{
        expected_recovered_sample, phasepoints, synth_signal, synth_signal_cn0,
    };
    use l1c_codes::{generate_l1c_pilot_code, L1C_CHIPS};

    #[test]
    fn test_generate_l1c_pilot_code_length() {
        let code = generate_l1c_pilot_code(1);
        assert_eq!(code.len(), L1C_CHIPS * 2);
    }

    #[test]
    fn test_generate_l1c_pilot_code_bipolar() {
        let code = generate_l1c_pilot_code(5);
        for &v in &code {
            assert!(v == 1.0 || v == -1.0, "unexpected value {v}");
        }
    }

    #[test]
    fn test_l1c_code_resampled_length() {
        let samples_per_code = (L1C_CHIPS * 2) as f64;
        let code = l1c_code_resampled(samples_per_code, 1, L1C_CHIPS * 4);
        assert_eq!(code.len(), L1C_CHIPS * 4); // 2 periods of BOC code
    }

    #[test]
    fn test_all_prns_generate() {
        for prn in 1..=L1C_NUM_SATS {
            let code = generate_l1c_pilot_code(prn);
            assert_eq!(
                code.len(),
                L1C_CHIPS * 2,
                "PRN {prn} wrong pilot code length"
            );
        }
    }

    #[test]
    fn test_acquire_l1c_single_recovers_delay() {
        let fs = 1.023e6; // one sample per chip; 10230 chips per 10 ms period
        let period = L1C_CHIPS;
        let code = l1c_code_resampled(period as f64, 1, 2 * period);

        let delay = 2500usize;
        let sig = synth_signal(&code, period, delay, 0.0, fs, 0.05, 31);
        let r = acquire_l1c_single(
            &sig, &phasepoints(fs, sig.len()), fs, 0.0, 3000.0, 1, period,
        );
        let recovered = (r.codephase_frac * period as f64).round() as usize;
        assert_eq!(recovered, expected_recovered_sample(period, delay));
        assert!(r.codephase_frac >= 0.0 && r.codephase_frac < 1.0);
    }

    #[test]
    fn test_acquire_l1c_single_noise_contrast() {
        let fs = 1.023e6;
        let period = L1C_CHIPS;
        let code = l1c_code_resampled(period as f64, 1, 2 * period);
        let n = code.len();
        let injected = synth_signal(&code, period, 700, 0.0, fs, 0.05, 32);
        let r_inj = acquire_l1c_single(
            &injected, &phasepoints(fs, n), fs, 0.0, 3000.0, 1, period,
        );
        let pure_noise: Vec<f32> = (0..n).map(|i| ((i as f64 * 1.3).sin() * 0.05) as f32).collect();
        let r_noise = acquire_l1c_single(
            &pure_noise, &phasepoints(fs, n), fs, 0.0, 3000.0, 1, period,
        );
        assert!(r_inj.signal_strength > 10.0 * r_noise.signal_strength);
    }

    #[test]
    fn test_acquire_non_period_multiple_no_panic() {
        // Signal lengths are NOT whole multiples of the code period.
        // Regression: rustfft used to panic with "Provided FFT buffer was
        // too small" because the local code replica could be resampled to
        // fewer samples than the FFT size.
        let period = L1C_CHIPS;
        let fs = 1.023e6;
        let delay = 1000usize;

        // Sweep signal lengths around whole code-period boundaries, including
        // lengths that are NOT whole multiples of the code period. This sweep
        // includes lengths (e.g. period + 6) that previously produced a code
        // replica one sample shorter than the FFT size.
        let mut lengths: Vec<usize> = Vec::new();
        lengths.extend((period - 3)..=(period + 12));
        lengths.extend((2 * period - 3)..=(2 * period + 4));

        for &n in &lengths {
            let code = l1c_code_resampled(period as f64, 1, n); // length n, partial period
            assert_eq!(code.len(), n, "code length {n}");
            // % n: code has length n; index within the code itself so the
            // signal stays in bounds even when n < period.
            let sig: Vec<f32> = (0..n).map(|i| code[(i + delay) % n] as f32).collect();
            let res = acquire_l1c_single(
                &sig, &phasepoints(fs, n), fs, 0.0, 3000.0, 1, period,
            );
            assert!(res.codephase_frac >= 0.0 && res.codephase_frac < 1.0, "n={n}");
            assert!(res.signal_strength > 0.0, "n={n}");
        }
    }

    #[test]
    fn test_cn0_acr_end_to_end() {
        // Median ACR C/N0 estimate over trials validates the 20 ms L1C
        // table against the real acquisition chain (synthesis → acquisition
        // → estimate_cn0).  Both replica hypotheses are exercised: no flip
        // (primary `[c, c]` wins, ±1 dB) and an overlay-code flip at the
        // code epoch (`[c, -c]` wins).  A flipped window is inherently
        // biased low by 0 to -6 dB depending on where the flip lands within
        // the period (median -2.5 dB), so the flip case asserts the median
        // in the [-4.0, -0.5] dB band below the injected value.
        let fs = 1.023e6; // 1 sample per chip
        let period = L1C_CHIPS;
        let f0 = 100e3;
        let cn0_true = 45.0;
        let code = l1c_code_resampled(period as f64, 1, 2 * period);
        let pp = phasepoints(fs, 2 * period);

        for (flip, n_trials, lo, hi) in [(false, 50usize, -1.0, 1.0), (true, 120usize, -4.0, -0.5)] {
            let mut ests: Vec<f64> = Vec::with_capacity(n_trials);
            for t in 0..n_trials {
                let delay = (t * 431) % period;
                let sig = synth_signal_cn0(
                    &code, period, delay, f0, fs, cn0_true, 0x1C1C + t as u64, flip,
                );
                let r = acquire_l1c_single(&sig, &pp, fs, f0, 0.0, 1, period);
                let r_a = (r.signal_strength / r.second_peak).powi(2);
                if let Some(c) = crate::acr::estimate_cn0(r_a, crate::acr::GPS_L1C_ACR_TABLE) {
                    ests.push(c);
                }
            }
            assert!(!ests.is_empty(), "no usable estimates (flip={flip})");
            let med = crate::stats::median(&ests);
            let d = med - cn0_true;
            assert!(
                d >= lo && d <= hi,
                "median C/N0 {med:.2} vs injected {cn0_true} (flip={flip}, n={})",
                ests.len()
            );
        }
    }
}
