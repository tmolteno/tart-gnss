// Copyright (c) 2026 Tim Molteno <tim@elec.ac.nz>
// SPDX-License-Identifier: GPL-3.0

//! SBAS (Satellite-Based Augmentation System) L1 C/A signal acquisition
//! via FFT-based circular cross-correlation.
//!
//! SBAS satellites transmit GPS L1 C/A codes with PRN numbers 120–138
//! (NMEA SV IDs 33–51).  The gold-code chip generation is reused from
//! the GPS acquisition module.

use num_complex::Complex;
use rayon::prelude::*;
use rustfft::FftPlanner;
use std::f64::consts::TAU;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::acquisition::{self, generate_ca_code_from_delay, gold_code_from_ca};
use crate::observation::Observation;

/// Number of SBAS PRNs (120..=138).
pub const SBAS_NUM_SATS: usize = 19;

/// SBAS L1 intermediate frequency in Hz (same as GPS L1 at the TART IF).
pub const SBAS_IF: f64 = 4.092e6;

/// Default search bandwidth in Hz.
pub const SBAS_SEARCH_BAND: f64 = 6000.0;

/// First SBAS PRN.
const SBAS_PRN_BASE: usize = 120;

// ---------------------------------------------------------------------------
// SBAS G2 delay table (PRN 120–138)
// ---------------------------------------------------------------------------

/// Code-phase delay table (G2 shift) for SBAS PRNs 120–138.
///
/// Values from IS-GPS-200 / GNSS-SDR.
const CODE_DELAY_TABLE: [usize; SBAS_NUM_SATS] = [
    145, 175, 52, 21, 237, 235, 886, 657, 634, 762, 355, 1012, 176, 603, 130, 359, 595, 68, 386,
];

// ---------------------------------------------------------------------------
// Output types
// ---------------------------------------------------------------------------

/// Result of an SBAS signal acquisition attempt for a single antenna.
#[derive(Debug, Clone, serde::Serialize)]
pub struct SbasAcquisitionResult {
    pub prn: usize,
    /// Peak correlation magnitude (normalised).
    pub signal_strength: f64,
    /// Code-phase offset in fractions of a millisecond [0, 1).
    pub codephase_frac: f64,
    /// Doppler frequency offset in Hz (relative to centre frequency).
    pub frequency: f64,
}

/// Per-SV SBAS acquisition result with per-antenna measurements.
#[derive(Debug, Clone, serde::Serialize)]
pub struct SbasPrnResult {
    /// Constellation label, e.g. "SBAS120".
    pub sv: String,
    /// Per-antenna signal strengths.
    pub strengths: Vec<f64>,
    /// Per-antenna code-phase offsets (fraction of a millisecond).
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

/// Collection of SBAS acquisition results for all PRNs.
#[derive(Debug, Clone, serde::Serialize)]
pub struct SbasAllAcquisitionOutput {
    pub results: Vec<SbasPrnResult>,
}

// ---------------------------------------------------------------------------
// SBAS code generation
// ---------------------------------------------------------------------------

/// Look up the G2 shift for an SBAS PRN (120-based, 0-indexed into the table).
fn sbas_g2shift(prn: usize) -> usize {
    CODE_DELAY_TABLE[prn - SBAS_PRN_BASE]
}

/// Generate the 1023-chip C/A gold code for an SBAS PRN.
pub fn generate_sbas_ca_code(prn: usize) -> [f64; acquisition::CA_CHIPS] {
    assert!(
        (SBAS_PRN_BASE..SBAS_PRN_BASE + SBAS_NUM_SATS).contains(&prn),
        "SBAS PRN must be in {SBAS_PRN_BASE}..={}, got {prn}",
        SBAS_PRN_BASE + SBAS_NUM_SATS - 1
    );
    generate_ca_code_from_delay(sbas_g2shift(prn))
}

/// Resample the SBAS C/A code to `samples_per_code` samples per period.
pub fn sbas_gold_code(samples_per_code: f64, prn: usize, epochs: f64) -> Vec<f64> {
    let ca = generate_sbas_ca_code(prn);
    gold_code_from_ca(samples_per_code, &ca, epochs)
}

// ---------------------------------------------------------------------------
// FFT-based acquisition (single PRN)
// ---------------------------------------------------------------------------

/// Perform parallel code-phase search (FFT circular cross-correlation) for a
/// single SBAS PRN over a frequency search band.
pub fn acquire_sbas(
    x: &[f64],
    sampling_freq: f64,
    center_freq: f64,
    search_band: f64,
    prn: usize,
) -> SbasAcquisitionResult {
    let sampling_period = 1.0 / sampling_freq;
    let samples_per_ms = sampling_freq / 1000.0;
    let samples_per_chunk = samples_per_ms as usize;

    let epochs_available = (x.len() as f64 / samples_per_ms).floor();
    let total_samples = (epochs_available * samples_per_chunk as f64) as usize;

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
    let code = sbas_gold_code(samples_per_ms, prn, epochs_available);
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

    let codephase_in_samples = best_codephase % samples_per_chunk;
    let codephase_frac = codephase_in_samples as f64 / samples_per_ms;
    let frequency = fc[best_freq_idx] - center_freq;

    SbasAcquisitionResult {
        prn,
        signal_strength: best_peak as f64,
        codephase_frac,
        frequency,
    }
}

// ---------------------------------------------------------------------------
// All-PRN search
// ---------------------------------------------------------------------------

/// Search for all SBAS L1 C/A PRNs across selected antennas.
///
/// If `ant_filter` is `Some(idx)`, only that antenna is used; otherwise all
/// antennas are searched.  PRN processing is parallelised via rayon.
pub fn acquire_all_sbas(
    obs: &Observation,
    center_freq: f64,
    search_band: f64,
    ant_filter: Option<usize>,
    debug: bool,
) -> SbasAllAcquisitionOutput {
    let sampling_freq = obs.get_sampling_rate();
    let n_ant = obs.config.num_antenna();
    let samples_per_ms = sampling_freq / 1000.0;
    let num_samples_per_ms = samples_per_ms as usize;
    let num_samples = 2 * num_samples_per_ms;

    let ant_indices: Vec<usize> = if let Some(ant) = ant_filter {
        vec![ant]
    } else {
        (0..n_ant).collect()
    };

    // Pre-extract and de-mean all antenna data.
    let ant_data: Vec<Vec<f64>> = ant_indices
        .iter()
        .map(|&ant_idx| {
            let bipolar = obs.get_antenna(ant_idx);
            let mean = bipolar.iter().sum::<f64>() / bipolar.len() as f64;
            let raw: Vec<f64> = bipolar.iter().map(|&v| v - mean).collect();
            raw[..num_samples.min(raw.len())].to_vec()
        })
        .collect();

    let prn_range: Vec<usize> = (SBAS_PRN_BASE..SBAS_PRN_BASE + SBAS_NUM_SATS).collect();
    let total = SBAS_NUM_SATS;
    let counter = AtomicUsize::new(0);

    let mut results: Vec<SbasPrnResult> = prn_range
        .into_par_iter()
        .map(|prn| {
            let mut strengths = Vec::with_capacity(ant_indices.len());
            let mut phases = Vec::with_capacity(ant_indices.len());
            let mut freqs = Vec::with_capacity(ant_indices.len());

            for (i, raw) in ant_data.iter().enumerate() {
                let result = acquire_sbas(raw, sampling_freq, center_freq, search_band, prn);

                if debug {
                    eprintln!(
                        "  sbas PRN {:3} ant {:2}: strength={:.3}  phase={:.6}  freq={:.1} Hz",
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
            eprintln!("  sbas [{n}/{total}]");
            if debug {
                eprintln!("  sbas [{n}/{total}] PRN {prn:03} complete");
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

            SbasPrnResult {
                sv: format!("SBAS{prn:03}"),
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

    results.sort_by_key(|r| r.sv.clone());

    SbasAllAcquisitionOutput { results }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_sbas_ca_code_length() {
        let ca = generate_sbas_ca_code(120);
        assert_eq!(ca.len(), acquisition::CA_CHIPS);
    }

    #[test]
    fn test_generate_sbas_ca_code_bipolar() {
        let ca = generate_sbas_ca_code(125);
        for &v in &ca {
            assert!(v == 1.0 || v == -1.0, "unexpected value {v}");
        }
    }

    #[test]
    fn test_sbas_gold_code_length() {
        let code = sbas_gold_code(1023.0, 120, 2.0);
        assert_eq!(code.len(), 2046);
    }

    #[test]
    #[should_panic]
    fn test_sbas_prn_out_of_range_low() {
        generate_sbas_ca_code(1); // GPS range, not SBAS
    }

    #[test]
    #[should_panic]
    fn test_sbas_prn_out_of_range_high() {
        generate_sbas_ca_code(200); // beyond SBAS range
    }
}
