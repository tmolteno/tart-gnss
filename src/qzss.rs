// Copyright (c) 2026 Tim Molteno <tim@elec.ac.nz>
// SPDX-License-Identifier: GPL-3.0

//! QZSS (Quasi-Zenith Satellite System) L1 C/A signal acquisition
//! via FFT-based circular cross-correlation.
//!
//! QZSS satellites transmit GPS-compatible L1 C/A codes with PRN numbers
//! 184–206.  The gold-code chip generation is reused from the GPS
//! acquisition module.

use num_complex::Complex;
use rayon::prelude::*;
use rustfft::FftPlanner;
use std::f64::consts::TAU;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::acquisition::{self, generate_ca_code_from_delay, gold_code_from_ca};
use crate::correlate::correlate_code;
use crate::observation::Observation;

/// Number of QZSS PRNs (184..=206).
pub const QZSS_NUM_SATS: usize = 23;

/// QZSS L1 intermediate frequency in Hz (same as GPS L1 at the TART IF).
pub const QZSS_IF: f64 = 4.092e6;

/// Default search bandwidth in Hz.
pub const QZSS_SEARCH_BAND: f64 = 6000.0;

/// First QZSS PRN.
const QZSS_PRN_BASE: usize = 184;

// ---------------------------------------------------------------------------
// QZSS G2 delay table (PRN 184–206)
// ---------------------------------------------------------------------------

/// Code-phase delay table (G2 shift) for QZSS PRNs 184–206.
///
/// Values from L1 C/A PRN Code Assignments (Jan 2026).
#[rustfmt::skip]
const CODE_DELAY_TABLE: [usize; QZSS_NUM_SATS] = [
    476, 193, 109, 445, 291,  87, 399, 292, 901, 339,
    208, 711, 189, 263, 537, 663, 942, 173, 900,  30,
    500, 935, 556,
];

// ---------------------------------------------------------------------------
// Output types
// ---------------------------------------------------------------------------

/// Result of a QZSS signal acquisition attempt for a single antenna.
#[derive(Debug, Clone, serde::Serialize)]
pub struct QzssAcquisitionResult {
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

/// Per-SV QZSS acquisition result with per-antenna measurements.
#[derive(Debug, Clone, serde::Serialize)]
pub struct QzssPrnResult {
    /// Constellation label, e.g. "QZSS184".
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

/// Collection of QZSS acquisition results for all PRNs.
#[derive(Debug, Clone, serde::Serialize)]
pub struct QzssAllAcquisitionOutput {
    pub antenna_numbers: Vec<usize>,
    pub results: Vec<QzssPrnResult>,
}

// ---------------------------------------------------------------------------
// QZSS code generation
// ---------------------------------------------------------------------------

/// Look up the G2 shift for a QZSS PRN (184-based, 0-indexed into the table).
fn qzss_g2shift(prn: usize) -> usize {
    CODE_DELAY_TABLE[prn - QZSS_PRN_BASE]
}

/// Generate the 1023-chip C/A gold code for a QZSS PRN.
pub fn generate_qzss_ca_code(prn: usize) -> [f64; acquisition::CA_CHIPS] {
    assert!(
        (QZSS_PRN_BASE..QZSS_PRN_BASE + QZSS_NUM_SATS).contains(&prn),
        "QZSS PRN must be in {QZSS_PRN_BASE}..={}, got {prn}",
        QZSS_PRN_BASE + QZSS_NUM_SATS - 1
    );
    generate_ca_code_from_delay(qzss_g2shift(prn))
}

/// Resample the QZSS C/A code to `samples_per_code` samples per period.
pub fn qzss_gold_code(samples_per_code: f64, prn: usize, epochs: f64) -> Vec<f64> {
    let ca = generate_qzss_ca_code(prn);
    gold_code_from_ca(samples_per_code, &ca, epochs)
}

// ---------------------------------------------------------------------------
// FFT-based acquisition (single PRN)
// ---------------------------------------------------------------------------

/// Perform parallel code-phase search (FFT circular cross-correlation) for a
/// single QZSS PRN over a frequency search band.
///
/// `signal_f32` and `phasepoints` are pre-computed by the caller to avoid
/// redundant work when searching multiple PRNs against the same antenna data.
pub fn acquire_qzss(
    signal_f32: &[f32],
    phasepoints: &[f64],
    sampling_freq: f64,
    center_freq: f64,
    search_band: f64,
    prn: usize,
    samples_per_chunk: usize,
) -> QzssAcquisitionResult {
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
    let code = qzss_gold_code(samples_per_ms, prn, epochs_available);
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

    QzssAcquisitionResult {
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

/// Search for all QZSS L1 C/A PRNs across selected antennas.
///
/// If `ant_filter` is `Some(antennas)`, only those antennas are used; otherwise all
/// antennas are searched.  PRN processing is parallelised via rayon.
pub fn acquire_all_qzss(
    obs: &Observation,
    center_freq: f64,
    search_band: f64,
    ant_filter: Option<Vec<usize>>,
    prn_filter: Option<&[usize]>,
    debug: bool,
    cn0: bool,
) -> QzssAllAcquisitionOutput {
    let sampling_freq = obs.get_sampling_rate();
    let n_ant = obs.config.num_antenna();
    let samples_per_ms = sampling_freq / 1000.0;
    let num_samples_per_ms = samples_per_ms as usize;
    let num_samples = 2 * num_samples_per_ms;

    let ant_indices: Vec<usize> = ant_filter.unwrap_or_else(|| (0..n_ant).collect());

    // Pre-extract, de-mean, convert to f32, and pre-compute phasepoints.
    let phase_const = TAU / sampling_freq;
    let ant_data: Vec<(Vec<f32>, Vec<f64>)> = ant_indices
        .iter()
        .map(|&ant_idx| {
            let bipolar = obs.get_antenna(ant_idx);
            let mean = bipolar.iter().sum::<f64>() / bipolar.len() as f64;
            let raw: Vec<f64> = bipolar.iter().map(|&v| v - mean).collect();
            let n = num_samples.min(raw.len());
            let signal_f32: Vec<f32> = raw[..n].iter().map(|&v| v as f32).collect();
            let phasepoints: Vec<f64> = (0..n).map(|i| phase_const * i as f64).collect();
            (signal_f32, phasepoints)
        })
        .collect();

    let prn_range: Vec<usize> = if let Some(filter) = prn_filter {
        filter.to_vec()
    } else {
        (QZSS_PRN_BASE..QZSS_PRN_BASE + QZSS_NUM_SATS).collect()
    };
    let total = prn_range.len();
    let counter = AtomicUsize::new(0);

    let mut results: Vec<QzssPrnResult> = prn_range
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
                let result = acquire_qzss(
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
                        "  qzss PRN {:3} ant {:2}: strength={:.3}  phase={:.6}  freq={:.1} Hz",
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
            eprintln!("  qzss [{n}/{total}]");
            if debug {
                eprintln!("  qzss [{n}/{total}] PRN {prn:03} complete");
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

            QzssPrnResult {
                sv: format!("QZSS{prn:03}"),
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

    QzssAllAcquisitionOutput { antenna_numbers: ant_indices, results }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::{expected_recovered_sample, phasepoints, synth_signal};

    #[test]
    fn test_generate_qzss_ca_code_length() {
        let ca = generate_qzss_ca_code(184);
        assert_eq!(ca.len(), acquisition::CA_CHIPS);
    }

    #[test]
    fn test_generate_qzss_ca_code_bipolar() {
        let ca = generate_qzss_ca_code(190);
        for &v in &ca {
            assert!(v == 1.0 || v == -1.0, "unexpected value {v}");
        }
    }

    #[test]
    fn test_qzss_gold_code_length() {
        let code = qzss_gold_code(1023.0, 184, 2.0);
        assert_eq!(code.len(), 2046);
    }

    #[test]
    #[should_panic]
    fn test_qzss_prn_out_of_range_low() {
        generate_qzss_ca_code(100);
    }

    #[test]
    #[should_panic]
    fn test_qzss_prn_out_of_range_high() {
        generate_qzss_ca_code(300);
    }

    #[test]
    fn test_all_qzss_prns_generate() {
        for prn in QZSS_PRN_BASE..QZSS_PRN_BASE + QZSS_NUM_SATS {
            let ca = generate_qzss_ca_code(prn);
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

    /// First 10 chips (octal) for each QZSS PRN from PRN assignment doc.
    /// Index 0 = PRN 184.
    ///
    /// Note: our C/A code generator uses -(G1 * G2) which produces XNOR
    /// (inverted relative to the IS-GPS-200 XOR definition).  The reference
    /// values below are therefore the *inverted* spec values, matching the
    /// "Initial G2 Setting" column rather than the "First 10 Chips" column.
    #[rustfmt::skip]
    const REF_FIRST10: [&str; 23] = [
        "1003", "1454", "1665", "0471", "1750", "0307", "0272", "0764",
        "1422", "1050", "1607", "1747", "1305", "0540", "1363", "0727",
        "0147", "1206", "1045", "0476", "0604", "1757", "1330",
    ];

    #[test]
    fn test_qzss_first10_octal() {
        for i in 0..QZSS_NUM_SATS {
            let prn = QZSS_PRN_BASE + i;
            let ca = generate_qzss_ca_code(prn);
            let octal = chips10_to_octal(&ca[0..10]);
            assert_eq!(
                octal, REF_FIRST10[i],
                "PRN {prn}: first 10 chips mismatch"
            );
        }
    }

    #[test]
    fn test_acquire_qzss_recovers_delay() {
        let fs = 1.023e6;
        let period = (fs / 1000.0) as usize; // 1023
        let code = qzss_gold_code(period as f64, 184, 2.0);

        let delay = 250usize;
        let sig = synth_signal(&code, period, delay, 0.0, fs, 0.05, 51);
        let r = acquire_qzss(
            &sig, &phasepoints(fs, sig.len()), fs, 0.0, 3000.0, 184, period,
        );
        let recovered = (r.codephase_frac * period as f64).round() as usize;
        assert_eq!(recovered, expected_recovered_sample(period, delay));
        assert!(r.codephase_frac >= 0.0 && r.codephase_frac < 1.0);
    }

    #[test]
    fn test_acquire_qzss_noise_contrast() {
        let fs = 1.023e6;
        let period = (fs / 1000.0) as usize;
        let code = qzss_gold_code(period as f64, 184, 2.0);
        let n = code.len();
        let injected = synth_signal(&code, period, 300, 0.0, fs, 0.05, 52);
        let r_inj = acquire_qzss(
            &injected, &phasepoints(fs, n), fs, 0.0, 3000.0, 184, period,
        );
        let pure_noise: Vec<f32> = (0..n).map(|i| ((i as f64 * 3.9).sin() * 0.05) as f32).collect();
        let r_noise = acquire_qzss(
            &pure_noise, &phasepoints(fs, n), fs, 0.0, 3000.0, 184, period,
        );
        assert!(r_inj.signal_strength > 10.0 * r_noise.signal_strength);
    }
}
