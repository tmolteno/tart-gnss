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
pub use l1c_codes::L1C_NUM_SATS;

use crate::observation::Observation;
use l1c_codes::{l1c_code_resampled, L1C_CODE_PERIOD};
use num_complex::Complex;
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
    pub results: Vec<L1CPrnResult>,
}

// ---------------------------------------------------------------------------
// FFT-based acquisition (single PRN)
// ---------------------------------------------------------------------------

/// Perform parallel code-phase search (FFT circular cross-correlation) for a
/// single GPS L1C PRN over a frequency search band.
///
/// Uses the BOC(1,1)-modulated pilot code as the local replica.
pub fn acquire_l1c_single(
    x: &[f64],
    sampling_freq: f64,
    center_freq: f64,
    search_band: f64,
    prn: usize,
) -> L1CAcquisitionResult {
    let sampling_period = 1.0 / sampling_freq;
    // L1C code period is 10 ms
    let samples_per_code_period = (sampling_freq * L1C_CODE_PERIOD) as usize;

    let epochs_available = (x.len() as f64 / samples_per_code_period as f64).floor();
    // Use at least 1 full code period of samples for the local replica,
    // but never exceed the actual available signal samples.
    let epochs = epochs_available.max(1.0);
    let ideal_samples = (epochs * samples_per_code_period as f64) as usize;
    let num_samples = ideal_samples.min(x.len());
    let effective_epochs = num_samples as f64 / samples_per_code_period as f64;

    // --- Frequency bins ----------------------------------------------------
    let freq_bin_size: f64 = 300.0;
    let n_freq_bins = (2.0 * search_band / freq_bin_size).round() as usize + 1;
    let fc: Vec<f64> = (0..n_freq_bins)
        .map(|i| {
            center_freq - search_band
                + (i as f64) * (2.0 * search_band) / (n_freq_bins as f64 - 1.0)
        })
        .collect();

    // --- FFT planner -------------------------------------------------------
    let mut planner: FftPlanner<f32> = FftPlanner::new();
    let fft = planner.plan_fft_forward(num_samples);
    let ifft = planner.plan_fft_inverse(num_samples);

    // --- Local code replica ------------------------------------------------
    // The BOC(1,1) code has 2 * L1C_CHIPS = 20460 samples per period
    let code = l1c_code_resampled(samples_per_code_period as f64, prn, effective_epochs);
    let mut code_complex: Vec<Complex<f32>> = code
        .iter()
        .map(|&v| Complex::new(v as f32, 0.0))
        .collect();
    fft.process(&mut code_complex);
    for c in &mut code_complex {
        *c = c.conj();
    }

    // --- Pre-compute phase ramp --------------------------------------------
    let phase_const = TAU * sampling_period;
    let phasepoints: Vec<f64> = (0..num_samples)
        .map(|i| phase_const * i as f64)
        .collect();

    // --- Per-frequency-bin correlation -------------------------------------
    let mut best_peak: f32 = f32::NEG_INFINITY;
    let mut best_freq_idx: usize = 0;
    let mut best_codephase: usize = 0;

    let signal_f32: Vec<f32> = x[..num_samples].iter().map(|&v| v as f32).collect();

    for (fi, &freq) in fc.iter().enumerate() {
        let mut iq: Vec<Complex<f32>> = phasepoints
            .iter()
            .zip(signal_f32.iter())
            .map(|(&p, &s)| {
                let phase = (p * freq) as f32;
                Complex::new(phase.cos(), phase.sin()) * Complex::new(s, 0.0)
            })
            .collect();

        fft.process(&mut iq);

        for (iq_elem, &code_elem) in iq.iter_mut().zip(code_complex.iter()) {
            *iq_elem *= code_elem;
        }

        ifft.process(&mut iq);

        let scale = 1.0 / (num_samples as f32).sqrt();
        let corr: Vec<f32> = iq.iter().map(|c| c.norm() * scale).collect();

        let (peak_idx, &peak_val) = corr
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap();

        if peak_val > best_peak {
            best_peak = peak_val;
            best_freq_idx = fi;
            best_codephase = peak_idx;
        }
    }

    let codephase_in_samples = best_codephase % samples_per_code_period;
    let codephase_frac = codephase_in_samples as f64 / samples_per_code_period as f64;
    let frequency = fc[best_freq_idx] - center_freq;

    L1CAcquisitionResult {
        prn,
        signal_strength: best_peak as f64,
        codephase_frac,
        frequency,
    }
}

// ---------------------------------------------------------------------------
// All-SV search
// ---------------------------------------------------------------------------

/// Search for all 63 GPS L1C SVs across selected antennas.
///
/// If `ant_filter` is `Some(idx)`, only that antenna is used; otherwise all
/// antennas are searched.  PRN processing is parallelised via rayon.
///
/// For each PRN, per-antenna signal strengths, code-phase offsets, and
/// Doppler frequency offsets are collected.
pub fn acquire_all_l1c(
    obs: &Observation,
    center_freq: f64,
    search_band: f64,
    ant_filter: Option<usize>,
    debug: bool,
) -> L1CAllAcquisitionOutput {
    let sampling_freq = obs.get_sampling_rate();
    let n_ant = obs.config.num_antenna();
    let num_samples_per_code_period = (sampling_freq * L1C_CODE_PERIOD) as usize;
    // Use 20 ms of data (2 code periods) per antenna
    let num_samples = 2 * num_samples_per_code_period;

    // Which antennas to search
    let ant_indices: Vec<usize> = if let Some(ant) = ant_filter {
        vec![ant]
    } else {
        (0..n_ant).collect()
    };

    // Pre-extract and de-mean all antenna data (avoids repeated HDF5 reads)
    let ant_data: Vec<Vec<f64>> = ant_indices
        .iter()
        .map(|&ant_idx| {
            let bipolar = obs.get_antenna(ant_idx);
            let mean = bipolar.iter().sum::<f64>() / bipolar.len() as f64;
            let raw: Vec<f64> = bipolar.iter().map(|&v| v - mean).collect();
            raw[..num_samples.min(raw.len())].to_vec()
        })
        .collect();

    // Parallel PRN search ----------------------------------------------------
    let total = L1C_NUM_SATS;
    let counter = AtomicUsize::new(0);

    let mut results: Vec<L1CPrnResult> = (1..=L1C_NUM_SATS)
        .into_par_iter()
        .map(|prn| {
            let mut strengths = Vec::with_capacity(ant_indices.len());
            let mut phases = Vec::with_capacity(ant_indices.len());
            let mut freqs = Vec::with_capacity(ant_indices.len());

            for (i, raw) in ant_data.iter().enumerate() {
                let result =
                    acquire_l1c_single(raw, sampling_freq, center_freq, search_band, prn);

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
                phases.push(result.codephase_frac);
                freqs.push(result.frequency);
            }

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
                phase_median,
                phase_mad,
                freq_median,
                freq_mad,
            }
        })
        .collect();

    // Sort by SV label to maintain ordering
    results.sort_by_key(|r| r.sv.clone());

    L1CAllAcquisitionOutput { results }
}

#[cfg(test)]
mod tests {
    use super::*;
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
        let code = l1c_code_resampled(samples_per_code, 1, 2.0);
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
}
