// Copyright (c) 2026 Tim Molteno <tim@elec.ac.nz>
// SPDX-License-Identifier: GPL-3.0

//! GPS L1 C/A code generation and signal acquisition via FFT-based
//! circular cross-correlation.
//!
//! Ported from `tart/tart/operation/acquisition.py`.

use num_complex::Complex;
use rayon::prelude::*;
use rustfft::FftPlanner;
use std::f64::consts::TAU;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::correlate::correlate_code;
use crate::observation::Observation;

/// Result of a GPS signal acquisition attempt for a single antenna.
#[derive(Debug, Clone, serde::Serialize)]
pub struct AcquisitionResult {
    pub prn: usize,
    /// Peak correlation magnitude (normalised).
    pub signal_strength: f64,
    /// Second-highest correlation magnitude at the same frequency,
    /// at least 1 chip away from the peak (for ACR C/N0 estimation).
    pub second_peak: f64,
    /// Code-phase offset in fractions of a millisecond [0, 1).
    pub codephase_frac: f64,
    /// Doppler frequency offset in Hz (relative to centre frequency).
    pub frequency: f64,
}

/// Per-SV GPS acquisition result with per-antenna measurements.
#[derive(Debug, Clone, serde::Serialize)]
pub struct GpsPrnResult {
    /// Constellation label, e.g. "GPS03".
    pub sv: String,
    /// Per-antenna signal strengths.
    pub strengths: Vec<f64>,
    /// Per-antenna code-phase offsets (fraction of a millisecond).
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

/// Collection of GPS acquisition results for all PRNs, grouped by PRN.
#[derive(Debug, Clone, serde::Serialize)]
pub struct GpsAllAcquisitionOutput {
    pub antenna_numbers: Vec<usize>,
    pub results: Vec<GpsPrnResult>,
}

/// Number of GPS PRNs with defined C/A codes.
pub const GPS_NUM_SATS: usize = 38;

/// GPS L1 intermediate frequency in Hz.
pub const GPS_IF: f64 = 4.092e6;

/// Default search bandwidth in Hz.
pub const GPS_SEARCH_BAND: f64 = 6000.0;

// ---------------------------------------------------------------------------
// GPS C/A code generation
// ---------------------------------------------------------------------------

/// Code-phase delay table (G2 shift) for GPS PRNs 1–38.
///
/// From the Danish GPS Centre lecture notes.
const CODE_DELAY_TABLE: [usize; 38] = [
    5, 6, 7, 8, 17, 18, 139, 140, 141, 251, 252, 254, 255, 256, 257, 258, 469, 470, 471, 472,
    473, 474, 509, 512, 513, 514, 515, 516, 859, 860, 861, 862, 863, 950, 947, 948, 950, 0,
];

/// Number of chips in one GPS L1 C/A code period.
pub const CA_CHIPS: usize = 1023;

/// Generate the 1023-chip C/A gold code for an arbitrary G2 shift.
///
/// `g2shift` is the number of chips to delay G2 relative to G1 (0..1022).
/// Used by both GPS and SBAS code generators.
pub fn generate_ca_code_from_delay(g2shift: usize) -> [f64; CA_CHIPS] {
    // --- G1 ----------------------------------------------------------------
    let mut g1 = [0.0f64; CA_CHIPS];
    let mut lfsr = [-1.0f64; 10]; // initialised to all -1 (bipolar)
    for i in 0..CA_CHIPS {
        g1[i] = lfsr[9];
        let save_bit = lfsr[2] * lfsr[9];
        lfsr.copy_within(0..9, 1);
        lfsr[0] = save_bit;
    }

    // --- G2 ----------------------------------------------------------------
    let mut g2 = [0.0f64; CA_CHIPS];
    lfsr = [-1.0f64; 10];
    for i in 0..CA_CHIPS {
        g2[i] = lfsr[9];
        let save_bit =
            lfsr[1] * lfsr[2] * lfsr[5] * lfsr[7] * lfsr[8] * lfsr[9];
        lfsr.copy_within(0..9, 1);
        lfsr[0] = save_bit;
    }

    // --- Shift G2 and combine ----------------------------------------------
    g2.rotate_right(g2shift);
    let mut ca = [0.0f64; CA_CHIPS];
    for i in 0..CA_CHIPS {
        ca[i] = -(g1[i] * g2[i]);
    }
    ca
}

/// Generate the 1023-chip GPS C/A gold code for a given PRN (1-based).
///
/// Returns a `[f64; CA_CHIPS]` array of ±1 values.
pub fn generate_ca_code(prn: usize) -> [f64; CA_CHIPS] {
    assert!(prn >= 1 && prn <= 38, "PRN must be in 1..=38, got {prn}");
    let g2shift = CODE_DELAY_TABLE[prn - 1];
    generate_ca_code_from_delay(g2shift)
}

/// Resample a pre-generated C/A code array to `samples_per_code` samples
/// per period, repeating for `epochs` full periods.
pub fn gold_code_from_ca(
    samples_per_code: f64,
    ca: &[f64; CA_CHIPS],
    epochs: f64,
) -> Vec<f64> {
    let samples_per_chip = samples_per_code / CA_CHIPS as f64;
    let num_samples = (samples_per_code * epochs).floor() as usize;

    (0..num_samples)
        .map(|n| {
            let idx = ((n as f64 / samples_per_chip).floor() as usize) % CA_CHIPS;
            ca[idx]
        })
        .collect()
}

/// Resample the C/A code to `samples_per_code` samples per period,
/// repeating for `epochs` full periods.
pub fn gold_code(samples_per_code: f64, prn: usize, epochs: f64) -> Vec<f64> {
    let ca = generate_ca_code(prn);
    gold_code_from_ca(samples_per_code, &ca, epochs)
}

// ---------------------------------------------------------------------------
// FFT-based acquisition
// ---------------------------------------------------------------------------

/// Perform parallel code-phase search (FFT circular cross-correlation) for a
/// single GPS PRN over a frequency search band.
///
/// `signal_f32` and `phasepoints` are pre-computed by the caller.
/// `samples_per_chunk` is the number of samples per 1-ms code period.
pub fn acquire_full(
    signal_f32: &[f32],
    phasepoints: &[f64],
    sampling_freq: f64,
    center_freq: f64,
    search_band: f64,
    prn: usize,
    samples_per_chunk: usize,
) -> AcquisitionResult {
    let num_samples = signal_f32.len();
    let samples_per_ms = samples_per_chunk as f64;
    let epochs_available = num_samples as f64 / samples_per_ms;

    // --- Frequency bins ----------------------------------------------------
    let freq_bin_size: f64 = 300.0;
    let n_freq_bins =
        (2.0 * search_band / freq_bin_size).round() as usize + 1;
    let fc: Vec<f64> = (0..n_freq_bins)
        .map(|i| center_freq - search_band + (i as f64) * (2.0 * search_band) / (n_freq_bins as f64 - 1.0))
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
    let code = gold_code(samples_per_ms, prn, epochs_available);
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
        samples_per_chunk,
        sampling_freq,
        fft,
        ifft,
    );

    let codephase_in_samples = peak.best_codephase % samples_per_chunk;
    let codephase_frac = codephase_in_samples as f64 / samples_per_ms;
    let frequency = fc[peak.best_freq_idx] - center_freq;

    AcquisitionResult {
        prn,
        signal_strength: peak.best_peak as f64,
        second_peak: peak.second_peak as f64,
        codephase_frac,
        frequency,
    }
}

// ---------------------------------------------------------------------------
// All-PRN search
// ---------------------------------------------------------------------------

/// Search for all GPS L1 C/A PRNs across selected antennas.
///
/// If `ant_filter` is `Some(antennas)`, only those antennas are used; otherwise all
/// antennas are searched.  PRN processing is parallelised via rayon.
///
/// For each PRN, per-antenna signal strengths, code-phase offsets, and
/// Doppler frequency offsets are collected.
pub fn acquire_all_gps(
    obs: &Observation,
    center_freq: f64,
    search_band: f64,
    ant_filter: Option<Vec<usize>>,
    prn_filter: Option<&[usize]>,
    debug: bool,
    cn0: bool,
) -> GpsAllAcquisitionOutput {
    let sampling_freq = obs.get_sampling_rate();
    let n_ant = obs.config.num_antenna();
    let samples_per_ms = sampling_freq / 1000.0;
    let num_samples_per_ms = samples_per_ms as usize;
    // Use 2 ms of data per antenna (matches the single-PRN mode).
    let num_samples = 2 * num_samples_per_ms;

    let ant_indices: Vec<usize> = ant_filter
        .unwrap_or_else(|| (0..n_ant).collect());

    // Pre-extract and de-mean all antenna data, convert to f32 once.
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

    let prn_list: Vec<usize> = if let Some(filter) = prn_filter {
        filter.to_vec()
    } else {
        (1..=GPS_NUM_SATS).collect()
    };
    let total = prn_list.len();
    let counter = AtomicUsize::new(0);

    let mut results: Vec<GpsPrnResult> = prn_list
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
                let result = acquire_full(signal_f32, phasepoints, sampling_freq, center_freq, search_band, prn, num_samples_per_ms);

                if debug {
                    eprintln!(
                        "  gps PRN {:2} ant {:2}: strength={:.3}  phase={:.6}  freq={:.1} Hz",
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
            eprintln!("  gps [{n}/{total}]");
            if debug {
                eprintln!("  gps [{n}/{total}] PRN {prn:02} complete");
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

            GpsPrnResult {
                sv: format!("GPS{prn:02}"),
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

    GpsAllAcquisitionOutput { antenna_numbers: ant_indices, results }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_ca_code_length() {
        let ca = generate_ca_code(1);
        assert_eq!(ca.len(), CA_CHIPS);
    }

    #[test]
    fn test_generate_ca_code_bipolar() {
        let ca = generate_ca_code(5);
        for &v in &ca {
            assert!(v == 1.0 || v == -1.0, "unexpected value {v}");
        }
    }

    #[test]
    fn test_gold_code_length() {
        let code = gold_code(1023.0, 1, 2.0);
        assert_eq!(code.len(), 2046);
    }
}
