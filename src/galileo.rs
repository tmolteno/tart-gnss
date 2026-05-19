// Copyright (c) 2026 Tim Molteno <tim@elec.ac.nz>
// SPDX-License-Identifier: GPL-3.0

//! Galileo E1 Open Service pilot (E1-C) signal acquisition via FFT-based
//! circular cross-correlation.
//!
//! The E1-C pilot channel uses 4092-chip primary memory codes (4 ms period
//! at 1.023 Mcps), with a 25-chip secondary code overlay.

#[path = "galileo_codes.rs"]
mod galileo_codes;

use crate::observation::Observation;
use galileo_codes::{GALILEO_E1_C_CODES, GALILEO_E1_CHIPS, GALILEO_E1_NUM_SATS};
use num_complex::Complex;
use rayon::prelude::*;
use rustfft::FftPlanner;
use serde::Serialize;
use std::f64::consts::TAU;

/// Result of a Galileo E1 signal acquisition attempt for a single satellite.
#[derive(Debug, Clone, Serialize)]
pub struct GalileoAcquisitionResult {
    pub prn: usize,
    /// Peak correlation magnitude (normalised).
    pub signal_strength: f64,
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
}

/// Collection of acquisition results for all Galileo SVs, grouped by PRN.
#[derive(Debug, Clone, Serialize)]
pub struct GalileoAllAcquisitionOutput {
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
pub fn acquire_galileo_single(
    x: &[f64],
    sampling_freq: f64,
    center_freq: f64,
    search_band: f64,
    prn: usize,
) -> GalileoAcquisitionResult {
    let sampling_period = 1.0 / sampling_freq;
    // Galileo E1 code period is 4 ms
    let samples_per_code_period = (sampling_freq * 0.004) as usize;

    let epochs_available = (x.len() as f64 / samples_per_code_period as f64).floor();
    // Use at least 1 full code period
    let epochs = epochs_available.max(1.0);
    let total_samples = (epochs * samples_per_code_period as f64) as usize;

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
    let fft = planner.plan_fft_forward(total_samples);
    let ifft = planner.plan_fft_inverse(total_samples);

    // --- Local code replica ------------------------------------------------
    let code = e1c_code_resampled(samples_per_code_period as f64, prn, epochs);
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
    let phasepoints: Vec<f64> = (0..total_samples)
        .map(|i| phase_const * i as f64)
        .collect();

    // --- Per-frequency-bin correlation -------------------------------------
    let signal_len = total_samples.min(x.len());
    let mut best_peak: f32 = f32::NEG_INFINITY;
    let mut best_freq_idx: usize = 0;
    let mut best_codephase: usize = 0;

    let signal_f32: Vec<f32> = x[..signal_len].iter().map(|&v| v as f32).collect();

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

        let scale = 1.0 / (signal_len as f32).sqrt();
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

    GalileoAcquisitionResult {
        prn,
        signal_strength: best_peak as f64,
        codephase_frac,
        frequency,
    }
}

// ---------------------------------------------------------------------------
// All-SV search
// ---------------------------------------------------------------------------

/// Search for all 50 Galileo E1-C SVs across selected antennas.
///
/// If `ant_filter` is `Some(idx)`, only that antenna is used; otherwise all
/// antennas are searched.  PRN processing is parallelised via rayon.
///
/// For each PRN, per-antenna signal strengths, code-phase offsets, and
/// Doppler frequency offsets are collected.
pub fn acquire_all_galileo(
    obs: &Observation,
    center_freq: f64,
    search_band: f64,
    ant_filter: Option<usize>,
) -> GalileoAllAcquisitionOutput {
    let sampling_freq = obs.get_sampling_rate();
    let n_ant = obs.config.num_antenna();
    let num_samples_per_code_period = (sampling_freq * 0.004) as usize;
    // Use 8 ms of data (2 code periods) per antenna
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
    let mut results: Vec<GalileoPrnResult> = (1..=GALILEO_E1_NUM_SATS)
        .into_par_iter()
        .map(|prn| {
            let mut strengths = Vec::with_capacity(ant_indices.len());
            let mut phases = Vec::with_capacity(ant_indices.len());
            let mut freqs = Vec::with_capacity(ant_indices.len());

            for (i, raw) in ant_data.iter().enumerate() {
                let result = acquire_galileo_single(
                    raw,
                    sampling_freq,
                    center_freq,
                    search_band,
                    prn,
                );

                eprintln!(
                    "  galileo PRN {:2} ant {:2}: strength={:.3}  phase={:.6}  freq={:.1} Hz",
                    prn,
                    ant_indices[i],
                    result.signal_strength,
                    result.codephase_frac,
                    result.frequency
                );

                strengths.push(result.signal_strength);
                phases.push(result.codephase_frac);
                freqs.push(result.frequency);
            }

            GalileoPrnResult {
                sv: format!("GSAT{prn:02}"),
                strengths,
                phases,
                freqs,
            }
        })
        .collect();

    results.sort_by_key(|r| r.sv.clone());

    GalileoAllAcquisitionOutput { results }
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
