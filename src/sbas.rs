// Copyright (c) 2026 Tim Molteno <tim@elec.ac.nz>
// SPDX-License-Identifier: GPL-3.0

//! SBAS (Satellite-Based Augmentation System) L1 C/A signal acquisition
//! via FFT-based circular cross-correlation.
//!
//! SBAS satellites transmit GPS L1 C/A codes with PRN numbers 120–158
//! (NMEA SV IDs 33–71).  The gold-code chip generation is reused from
//! the GPS acquisition module.

use num_complex::Complex;
use rayon::prelude::*;
use rustfft::FftPlanner;
use std::f64::consts::TAU;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::acquisition::{self, generate_ca_code_from_delay, gold_code_from_ca};
use crate::observation::Observation;

/// Number of SBAS PRNs (120..=158).
pub const SBAS_NUM_SATS: usize = 39;

/// SBAS L1 intermediate frequency in Hz (same as GPS L1 at the TART IF).
pub const SBAS_IF: f64 = 4.092e6;

/// Default search bandwidth in Hz.
pub const SBAS_SEARCH_BAND: f64 = 6000.0;

/// First SBAS PRN.
const SBAS_PRN_BASE: usize = 120;

// ---------------------------------------------------------------------------
// SBAS G2 delay table (PRN 120–158)
// ---------------------------------------------------------------------------

/// Code-phase delay table (G2 shift) for SBAS PRNs 120–158.
///
/// Values from L1 C/A PRN Code Assignments (Jan 2026).
const CODE_DELAY_TABLE: [usize; SBAS_NUM_SATS] = [
    145, 175,  52,  21, 237, 235, 886, 657, 634, 762,
    355, 1012, 176, 603, 130, 359, 595,  68, 386, 797,
    456, 499,  883, 307, 127, 211, 121, 118, 163, 628,
    853, 484,  289, 811, 202, 1021, 463, 568, 904,
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
    /// Second-highest correlation magnitude at the same frequency,
    /// at least 1 chip away from the peak (for ACR C/N0 estimation).
    pub second_peak: f64,
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
    signal_f32: &[f32],
    phasepoints: &[f64],
    sampling_freq: f64,
    center_freq: f64,
    search_band: f64,
    prn: usize,
    samples_per_chunk: usize,
) -> SbasAcquisitionResult {
    let num_samples = signal_f32.len();
    let samples_per_ms = samples_per_chunk as f64;
    let epochs_available = num_samples as f64 / samples_per_ms;

    // --- Frequency bins ----------------------------------------------------
    let freq_bin_size: f64 = 300.0;
    let n_freq_bins = (2.0 * search_band / freq_bin_size).round() as usize + 1;
    let fc: Vec<f64> = (0..n_freq_bins)
        .map(|i| {
            center_freq - search_band
                + (i as f64) * (2.0 * search_band) / (n_freq_bins as f64 - 1.0)
        })
        .collect();

    // --- FFT planner (thread-local) ---------------------------------------
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

    // --- Pre-allocated working buffers ------------------------------------
    let mut iq_buf: Vec<Complex<f32>> = vec![Complex::new(0.0, 0.0); num_samples];
    let mut corr_buf: Vec<f32> = vec![0.0; num_samples];

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

    // --- Per-frequency-bin correlation -------------------------------------
    let mut best_peak: f32 = f32::NEG_INFINITY;
    let mut best_freq_idx: usize = 0;
    let mut best_codephase: usize = 0;
    let mut best_corr: Vec<f32> = Vec::new();

    for (fi, &freq) in fc.iter().enumerate() {
        for (idx, (&p, &s)) in phasepoints.iter().zip(signal_f32.iter()).enumerate() {
            let phase = (p * freq) as f32;
            iq_buf[idx] = Complex::new(phase.cos(), phase.sin()) * Complex::new(s, 0.0);
        }

        fft.process(&mut iq_buf);

        for (iq_elem, &code_elem) in iq_buf.iter_mut().zip(code_complex.iter()) {
            *iq_elem *= code_elem;
        }

        ifft.process(&mut iq_buf);

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
    let main_codephase = best_codephase % samples_per_chunk;
    let mut second_peak: f32 = f32::NEG_INFINITY;

    for (idx, &val) in best_corr.iter().enumerate() {
        let idx_cp = idx % samples_per_chunk;
        let dist = if idx_cp >= main_codephase {
            (idx_cp - main_codephase).min(samples_per_chunk - idx_cp + main_codephase)
        } else {
            (main_codephase - idx_cp).min(samples_per_chunk - main_codephase + idx_cp)
        };
        if dist > chip_width && val > second_peak {
            second_peak = val;
        }
    }

    if second_peak <= 0.0 {
        second_peak = 1e-6;
    }

    let codephase_in_samples = best_codephase % samples_per_chunk;
    let codephase_frac = codephase_in_samples as f64 / samples_per_ms;
    let frequency = fc[best_freq_idx] - center_freq;

    SbasAcquisitionResult {
        prn,
        signal_strength: best_peak as f64,
        second_peak: second_peak as f64,
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
    prn_filter: Option<&[usize]>,
    debug: bool,
    cn0: bool,
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

    // Pre-extract, de-mean, convert to f32, and pre-compute phasepoints.
    let phase_const = TAU / sampling_freq;
    let ant_data: Vec<(Vec<f32>, Vec<f64>)> = ant_indices
        .iter()
        .map(|&ant_idx| {
            let bipolar = obs.get_antenna(ant_idx);
            let mean = bipolar.iter().sum::<f64>() / bipolar.len() as f64;
            let raw: Vec<f64> = bipolar.iter().map(|&v| v - mean).collect();
            let raw = &raw[..num_samples.min(raw.len())];
            let signal_f32: Vec<f32> = raw.iter().map(|&v| v as f32).collect();
            let phasepoints: Vec<f64> = (0..signal_f32.len())
                .map(|i| phase_const * i as f64)
                .collect();
            (signal_f32, phasepoints)
        })
        .collect();

    let prn_range: Vec<usize> = if let Some(filter) = prn_filter {
        filter.to_vec()
    } else {
        (SBAS_PRN_BASE..SBAS_PRN_BASE + SBAS_NUM_SATS).collect()
    };
    let total = prn_range.len();
    let counter = AtomicUsize::new(0);

    let mut results: Vec<SbasPrnResult> = prn_range
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
                let result = acquire_sbas(
                    signal_f32,
                    phasepoints,
                    sampling_freq,
                    center_freq,
                    search_band,
                    prn,
                    num_samples_per_ms,
                );

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
                cn0_acr,
                phase_median,
                phase_mad,
                freq_median,
                freq_mad,
            }
        })
        .collect();

    results.sort_by(|a, b| a.sv.cmp(&b.sv));

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

    #[test]
    fn test_all_sbas_prns_generate() {
        for prn in SBAS_PRN_BASE..SBAS_PRN_BASE + SBAS_NUM_SATS {
            let ca = generate_sbas_ca_code(prn);
            assert_eq!(
                ca.len(),
                acquisition::CA_CHIPS,
                "PRN {prn} wrong code length"
            );
            for &v in &ca {
                assert!(v == 1.0 || v == -1.0, "PRN {prn}: unexpected value {v}");
            }
        }
    }

    /// Convert 10 chips (±1) to an octal string (IS-GPS-200 convention).
    ///
    /// The 10-bit sequence is zero-padded to 12 bits (2 leading zeros),
    /// then split into 4 groups of 3 bits, each producing an octal digit.
    /// +1→0, -1→1 (C/A convention: 0 is in-phase).
    fn chips10_to_octal(chips: &[f64]) -> String {
        assert_eq!(chips.len(), 10);
        let bits: Vec<u8> = chips.iter().map(|&v| if v == 1.0 { 0 } else { 1 }).collect();
        // 10 bits → pad to 12 with two leading zeros, then 4 groups of 3
        let d0 = bits[0];
        let d1 = 4 * bits[1] + 2 * bits[2] + bits[3];
        let d2 = 4 * bits[4] + 2 * bits[5] + bits[6];
        let d3 = 4 * bits[7] + 2 * bits[8] + bits[9];
        format!("{d0}{d1}{d2}{d3}")
    }

    /// First 10 chips (octal) for each SBAS PRN from PRN assignment doc.
    /// Index 0 = PRN 120.
    ///
    /// Note: our C/A code generator uses -(G1 * G2) which produces XNOR
    /// (inverted relative to the IS-GPS-200 XOR definition).  The reference
    /// values below are therefore the *inverted* spec values, matching the
    /// "Initial G2 Setting" column rather than the "First 10 Chips" column.
    #[rustfmt::skip]
    const REF_FIRST10: [&str; 39] = [
        "1106", "1241", "0267", "0232", "1617", "1076", "1764", "0717",
        "1532", "1250", "0341", "0551", "0520", "1731", "0706", "1216",
        "0740", "1007", "0450", "0305", "1653", "1411", "1644", "1312",
        "1060", "1560", "0035", "0355", "0335", "1254", "1041", "0142",
        "1641", "1504", "0751", "1774", "0107", "1153", "1542",
    ];

    #[test]
    fn test_sbas_first10_octal() {
        for i in 0..SBAS_NUM_SATS {
            let prn = SBAS_PRN_BASE + i;
            let ca = generate_sbas_ca_code(prn);
            let octal = chips10_to_octal(&ca[0..10]);
            assert_eq!(
                octal, REF_FIRST10[i],
                "PRN {prn}: first 10 chips mismatch"
            );
        }
    }
}
