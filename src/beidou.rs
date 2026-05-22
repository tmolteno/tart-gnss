// Copyright (c) 2026 Tim Molteno <tim@elec.ac.nz>
// SPDX-License-Identifier: GPL-3.0

//! BeiDou B1C pilot (B1C_pilot) signal acquisition via FFT-based
//! circular cross-correlation.
//!
//! The B1C pilot channel uses 10230-chip Weil primary codes (10 ms period
//! at 1.023 Mcps), with BOC(1,1) modulation applied for acquisition.

#[path = "beidou_codes.rs"]
mod beidou_codes;

use crate::observation::Observation;
use beidou_codes::{
    b1c_code_resampled, BEIDOU_B1C_CODE_PERIOD, BEIDOU_B1C_NUM_SATS,
};
use num_complex::Complex;
use rayon::prelude::*;
use rustfft::FftPlanner;
use serde::Serialize;
use std::f64::consts::TAU;
use std::sync::atomic::{AtomicUsize, Ordering};

/// Result of a BeiDou B1C signal acquisition attempt for a single satellite.
#[derive(Debug, Clone, Serialize)]
pub struct BeiDouAcquisitionResult {
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

/// Per-SV BeiDou acquisition result with per-antenna measurements.
#[derive(Debug, Clone, Serialize)]
pub struct BeiDouPrnResult {
    /// Constellation label, e.g. "BEIDOU03".
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

/// Collection of acquisition results for all BeiDou B1C SVs, grouped by PRN.
#[derive(Debug, Clone, Serialize)]
pub struct BeiDouAllAcquisitionOutput {
    pub results: Vec<BeiDouPrnResult>,
}

// ---------------------------------------------------------------------------
// FFT-based acquisition (single PRN)
// ---------------------------------------------------------------------------

/// Perform parallel code-phase search (FFT circular cross-correlation) for a
/// single BeiDou B1C PRN over a frequency search band.
///
/// Uses the BOC(1,1)-modulated pilot code as the local replica.
///
/// `signal_f32` and `phasepoints` are pre-computed by the caller to avoid
/// redundant allocations across PRN/antenna iterations.
pub fn acquire_beidou_single(
    signal_f32: &[f32],
    phasepoints: &[f64],
    sampling_freq: f64,
    center_freq: f64,
    search_band: f64,
    prn: usize,
    samples_per_code_period: usize,
) -> BeiDouAcquisitionResult {
    let num_samples = signal_f32.len();

    let epochs_available = (num_samples as f64 / samples_per_code_period as f64).floor();
    let effective_epochs = epochs_available.max(1.0);

    // --- Frequency bins ----------------------------------------------------
    let freq_bin_size: f64 = 300.0;
    let n_freq_bins = (2.0 * search_band / freq_bin_size).round() as usize + 1;
    let fc: Vec<f64> = (0..n_freq_bins)
        .map(|i| {
            center_freq - search_band
                + (i as f64) * (2.0 * search_band) / (n_freq_bins as f64 - 1.0)
        })
        .collect();

    // --- FFT planner (thread-local, reused across calls) -------------------
    thread_local! {
        static PLANNER: std::cell::RefCell<FftPlanner<f32>> =
            std::cell::RefCell::new(FftPlanner::new());
    }

    // --- Local code replica ------------------------------------------------
    let code = b1c_code_resampled(samples_per_code_period as f64, prn, effective_epochs);
    let mut code_complex: Vec<Complex<f32>> = code
        .iter()
        .map(|&v| Complex::new(v as f32, 0.0))
        .collect();

    PLANNER.with(|p| {
        let mut planner = p.borrow_mut();
        let fft = planner.plan_fft_forward(num_samples);
        fft.process(&mut code_complex);
    });
    for c in &mut code_complex {
        *c = c.conj();
    }

    // --- Per-frequency-bin correlation (pre-allocated buffers) -------------
    let mut best_peak: f32 = f32::NEG_INFINITY;
    let mut best_freq_idx: usize = 0;
    let mut best_codephase: usize = 0;
    let mut best_corr: Vec<f32> = Vec::new();

    let mut iq_buf: Vec<Complex<f32>> = vec![Complex::new(0.0, 0.0); num_samples];
    let mut corr_buf: Vec<f32> = vec![0.0; num_samples];

    for (fi, &freq) in fc.iter().enumerate() {
        for (idx, (&p, &s)) in phasepoints.iter().zip(signal_f32.iter()).enumerate() {
            let phase = (p * freq) as f32;
            iq_buf[idx] = Complex::new(phase.cos(), phase.sin()) * Complex::new(s, 0.0);
        }

        PLANNER.with(|p| {
            let mut planner = p.borrow_mut();
            let fft = planner.plan_fft_forward(num_samples);
            fft.process(&mut iq_buf);
        });

        for (iq_elem, &code_elem) in iq_buf.iter_mut().zip(code_complex.iter()) {
            *iq_elem *= code_elem;
        }

        PLANNER.with(|p| {
            let mut planner = p.borrow_mut();
            let ifft = planner.plan_fft_inverse(num_samples);
            ifft.process(&mut iq_buf);
        });

        let scale = 1.0 / (num_samples as f32).sqrt();
        for (idx, c) in iq_buf.iter().enumerate() {
            corr_buf[idx] = c.norm() * scale;
        }

        let (peak_idx, &peak_val) = corr_buf
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap();

        if peak_val > best_peak {
            best_peak = peak_val;
            best_freq_idx = fi;
            best_codephase = peak_idx;
            best_corr = corr_buf.clone();
        }
    }

    // --- Find second peak (>1 chip from main peak) -------------------------
    let chip_width = (sampling_freq / 1.023e6).ceil() as usize;
    let main_codephase = best_codephase % samples_per_code_period;
    let mut second_peak: f32 = f32::NEG_INFINITY;

    for (idx, &val) in best_corr.iter().enumerate() {
        let idx_cp = idx % samples_per_code_period;
        let dist = if idx_cp >= main_codephase {
            (idx_cp - main_codephase).min(samples_per_code_period - idx_cp + main_codephase)
        } else {
            (main_codephase - idx_cp).min(samples_per_code_period - main_codephase + idx_cp)
        };
        if dist > chip_width && val > second_peak {
            second_peak = val;
        }
    }

    if second_peak <= 0.0 {
        second_peak = 1e-6;
    }

    let codephase_in_samples = best_codephase % samples_per_code_period;
    let codephase_frac = codephase_in_samples as f64 / samples_per_code_period as f64;
    let frequency = fc[best_freq_idx] - center_freq;

    BeiDouAcquisitionResult {
        prn,
        signal_strength: best_peak as f64,
        second_peak: second_peak as f64,
        codephase_frac,
        frequency,
    }
}

// ---------------------------------------------------------------------------
// All-SV search
// ---------------------------------------------------------------------------

/// Search for all 63 BeiDou B1C SVs across selected antennas.
///
/// If `ant_filter` is `Some(idx)`, only that antenna is used; otherwise all
/// antennas are searched.  PRN processing is parallelised via rayon.
///
/// For each PRN, per-antenna signal strengths, code-phase offsets, and
/// Doppler frequency offsets are collected.
pub fn acquire_all_beidou(
    obs: &Observation,
    center_freq: f64,
    search_band: f64,
    ant_filter: Option<usize>,
    prn_filter: Option<&[usize]>,
    debug: bool,
    cn0: bool,
) -> BeiDouAllAcquisitionOutput {
    let sampling_freq = obs.get_sampling_rate();
    let n_ant = obs.config.num_antenna();
    let num_samples_per_code_period = (sampling_freq * BEIDOU_B1C_CODE_PERIOD) as usize;
    // Use 20 ms of data (2 code periods) per antenna
    let num_samples = 2 * num_samples_per_code_period;

    // Which antennas to search
    let ant_indices: Vec<usize> = if let Some(ant) = ant_filter {
        vec![ant]
    } else {
        (0..n_ant).collect()
    };

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
        (1..=BEIDOU_B1C_NUM_SATS).collect()
    };
    let total = prn_list.len();
    let counter = AtomicUsize::new(0);

    let mut results: Vec<BeiDouPrnResult> = prn_list
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
                    acquire_beidou_single(signal_f32, phasepoints, sampling_freq, center_freq, search_band, prn, num_samples_per_code_period);

                if debug {
                    eprintln!(
                        "  beidou PRN {:2} ant {:2}: strength={:.3}  phase={:.6}  freq={:.1} Hz",
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
            eprintln!("  beidou [{n}/{total}]");
            if debug {
                eprintln!("  beidou [{n}/{total}] PRN {prn:02} complete");
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

            BeiDouPrnResult {
                sv: format!("BEIDOU{prn:02}"),
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

    BeiDouAllAcquisitionOutput { results }
}

#[cfg(test)]
mod tests {
    use super::*;
    use beidou_codes::{generate_b1c_pilot_code, BEIDOU_B1C_CHIPS};

    #[test]
    fn test_generate_b1c_pilot_code_length() {
        let code = generate_b1c_pilot_code(1);
        assert_eq!(code.len(), BEIDOU_B1C_CHIPS * 2);
    }

    #[test]
    fn test_generate_b1c_pilot_code_bipolar() {
        let code = generate_b1c_pilot_code(5);
        for &v in &code {
            assert!(v == 1.0 || v == -1.0, "unexpected value {v}");
        }
    }

    #[test]
    fn test_b1c_code_resampled_length() {
        let samples_per_code = (BEIDOU_B1C_CHIPS * 2) as f64;
        let code = b1c_code_resampled(samples_per_code, 1, 2.0);
        assert_eq!(code.len(), BEIDOU_B1C_CHIPS * 4); // 2 periods of BOC code
    }

    #[test]
    fn test_all_prns_generate() {
        for prn in 1..=BEIDOU_B1C_NUM_SATS {
            let code = generate_b1c_pilot_code(prn);
            assert_eq!(
                code.len(),
                BEIDOU_B1C_CHIPS * 2,
                "PRN {prn} wrong pilot code length"
            );
        }
    }
}
