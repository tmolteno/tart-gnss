// Copyright (c) 2026 Tim Molteno <tim@elec.ac.nz>
// SPDX-License-Identifier: GPL-3.0

//! Shared FFT-based parallel code-phase search (circular cross-correlation)
//! used by all GNSS signal acquisition modules.
//!
//! Each frequency bin is processed independently and in parallel via rayon,
//! replacing the sequential per-bin loop that was previously duplicated in
//! every constellation module.
//!
//! The search is dual-hypothesis: it correlates against two local replicas —
//! the two-period code `[c, c]` and the same code with its second period
//! negated, `[c, -c]`.  A sign flip between the two integrated periods
//! (pilot secondary codes such as the Galileo E1-C 25-chip sequence or the
//! GPS L1C overlay code, or navigation data bits on C/A signals) only ever
//! changes the *relative* sign of the two periods, so one of the two
//! hypotheses always matches and the stronger correlation wins.  The
//! wiped-signal forward FFT is shared between the hypotheses; only the
//! conjugate multiply and inverse FFT are repeated.

use num_complex::Complex;
use rayon::prelude::*;
use rustfft::Fft;
use std::sync::Arc;

/// Result of a parallel code-phase search over multiple frequency bins.
#[derive(Debug, Clone)]
pub(crate) struct CorrelationPeak {
    /// Peak correlation magnitude (normalised).
    pub best_peak: f32,
    /// Second-highest correlation magnitude at the same frequency,
    /// at least 1 chip away from the main peak (for ACR C/N0 estimation).
    pub second_peak: f32,
    /// Index of the best frequency bin (into the `fc` slice passed by the caller).
    pub best_freq_idx: usize,
    /// Sample index of the correlation peak.
    pub best_codephase: usize,
}

/// Build the pre-FFT'd, conjugated spectra of the two replica hypotheses
/// for a two-period `code`: `[c, c]` and `[c, -c]` (second period negated).
///
/// The negated hypothesis matches a sign flip between the two integrated
/// periods (pilot secondary code or navigation data bit).  For a signal
/// shorter than two periods the negation is a no-op and the two spectra are
/// identical.
pub(crate) fn code_spectra(
    code: &[f64],
    period: usize,
    fft: &Arc<dyn Fft<f32>>,
) -> (Vec<Complex<f32>>, Vec<Complex<f32>>) {
    let mut alt = code.to_vec();
    let split = period.min(alt.len());
    for c in &mut alt[split..] {
        *c = -*c;
    }
    (code_fft_conj(fft, code), code_fft_conj(fft, &alt))
}

fn code_fft_conj(fft: &Arc<dyn Fft<f32>>, code: &[f64]) -> Vec<Complex<f32>> {
    let mut cc: Vec<Complex<f32>> =
        code.iter().map(|&v| Complex::new(v as f32, 0.0)).collect();
    fft.process(&mut cc);
    for c in &mut cc {
        *c = c.conj();
    }
    cc
}

/// Perform parallel code-phase search over the given frequency bins.
///
/// Each frequency bin is processed independently and in parallel via rayon.
/// Per-bin working buffers (`iq` and `corr`) are allocated inside the parallel
/// closure; peak memory is bounded by `num_threads × 2 × num_samples × 8 B`.
///
/// # Arguments
/// * `signal_f32` — de-meaned antenna samples (real-valued, f32)
/// * `phasepoints` — pre-computed phase accumulators for each sample (`TAU * i / fs`)
/// * `code_complex` — pre-FFT'd and conjugated local code replica `[c, c]`
///   (length == num_samples)
/// * `code_complex_alt` — pre-FFT'd and conjugated replica `[c, -c]`
///   (second period negated), the secondary-code/data-bit flip hypothesis
/// * `fc` — frequency bin centre values in Hz
/// * `num_samples` — number of samples in the signal
/// * `code_period_samples` — length of one full code period in samples
///   (e.g. samples_per_ms for GPS C/A, samples_per_code_period for Galileo/BeiDou)
/// * `sampling_freq` — sample rate in Hz (used for chip-width calculation)
/// * `fft` — pre-planned forward FFT, shared across all bins
/// * `ifft` — pre-planned inverse FFT, shared across all bins
#[allow(clippy::too_many_arguments)]
pub(crate) fn correlate_code(
    signal_f32: &[f32],
    phasepoints: &[f64],
    code_complex: &[Complex<f32>],
    code_complex_alt: &[Complex<f32>],
    fc: &[f64],
    num_samples: usize,
    code_period_samples: usize,
    sampling_freq: f64,
    fft: Arc<dyn Fft<f32>>,
    ifft: Arc<dyn Fft<f32>>,
) -> CorrelationPeak {
    assert!(!fc.is_empty(), "frequency bin list must be non-empty");

    // Process each frequency bin in parallel via rayon
    let results: Vec<(usize, usize, f32, Vec<f32>)> = fc
        .par_iter()
        .enumerate()
        .map(|(fi, &freq)| {
            // 1. IQ carrier wipe-off (element-wise trig + multiply)
            let mut iq = vec![Complex::new(0.0f32, 0.0); num_samples];
            for (idx, (&p, &s)) in phasepoints.iter().zip(signal_f32.iter()).enumerate() {
                let phase = (p * freq) as f32;
                iq[idx] = Complex::new(phase.cos(), phase.sin()) * Complex::new(s, 0.0);
            }

            // 2. Forward FFT — shared by both replica hypotheses
            fft.process(&mut iq);

            // 3-6. Correlate against each replica hypothesis and keep the
            // stronger one.  The sign flip between the two integrated
            // periods selects either `[c, c]` (no flip) or `[c, -c]` (flip).
            let scale = 1.0 / (num_samples as f32).sqrt();
            let mut best: Option<(usize, f32, Vec<f32>)> = None;
            for replica in [code_complex, code_complex_alt] {
                let mut iq = iq.clone();
                for (iq_elem, &code_elem) in iq.iter_mut().zip(replica.iter()) {
                    *iq_elem *= code_elem;
                }
                ifft.process(&mut iq);
                let mut corr = vec![0.0f32; num_samples];
                for (idx, c) in iq.iter().enumerate() {
                    corr[idx] = c.norm() * scale;
                }
                let (peak_idx, &peak_val) = corr
                    .iter()
                    .enumerate()
                    .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                    .unwrap();
                if best.as_ref().map_or(true, |(_, v, _)| peak_val > *v) {
                    best = Some((peak_idx, peak_val, corr));
                }
            }
            let (peak_idx, peak_val, corr) = best.expect("two hypotheses always present");

            (fi, peak_idx, peak_val, corr)
        })
        .collect();

    // Find global maximum across all frequency bins
    let (best_freq_idx, best_codephase, best_peak, best_corr) = results
        .into_iter()
        .max_by(|(_, _, a, _), (_, _, b, _)| a.partial_cmp(b).unwrap())
        .expect("frequency bin list must be non-empty");

    // --- Second peak search (>1 chip away from main peak) ------------------
    let chip_width = (sampling_freq / 1.023e6).ceil() as usize;
    let main_codephase = best_codephase % code_period_samples;
    let mut second_peak: f32 = f32::NEG_INFINITY;

    for (idx, &val) in best_corr.iter().enumerate() {
        let idx_cp = idx % code_period_samples;
        let dist = if idx_cp >= main_codephase {
            (idx_cp - main_codephase).min(code_period_samples - idx_cp + main_codephase)
        } else {
            (main_codephase - idx_cp).min(code_period_samples - main_codephase + idx_cp)
        };
        if dist > chip_width && val > second_peak {
            second_peak = val;
        }
    }

    if second_peak <= 0.0 {
        second_peak = 1e-6;
    }

    CorrelationPeak {
        best_peak,
        second_peak,
        best_freq_idx,
        best_codephase,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::acquisition::gold_code;
    use crate::testutil::{expected_recovered_sample, phasepoints, synth_signal};
    use rustfft::FftPlanner;

    /// Build the pre-FFT'd, conjugated code replica exactly as the
    /// acquisition modules do.
    fn code_fft_conj(fft: Arc<dyn Fft<f32>>, code: &[f64]) -> Vec<Complex<f32>> {
        let mut cc: Vec<Complex<f32>> =
            code.iter().map(|&v| Complex::new(v as f32, 0.0)).collect();
        fft.process(&mut cc);
        for c in &mut cc {
            *c = c.conj();
        }
        cc
    }

    /// The flip-hypothesis replica `[c, -c]`: second period negated.
    fn code_alt_fft_conj(fft: Arc<dyn Fft<f32>>, code: &[f64], period: usize) -> Vec<Complex<f32>> {
        let mut alt = code.to_vec();
        for c in &mut alt[period..] {
            *c = -*c;
        }
        code_fft_conj(fft, &alt)
    }

    #[test]
    fn test_correlate_code_recovers_delay() {
        let period = 256usize;
        let code = gold_code(period as f64, 1, period); // 256-sample GPS C/A code
        let fs = 1.023e6;

        let mut planner = FftPlanner::<f32>::new();
        let fft = planner.plan_fft_forward(period);
        let ifft = planner.plan_fft_inverse(period);

        let cc = code_fft_conj(fft.clone(), &code);

        for delay in [0usize, 50, 128, 200] {
            let sig = synth_signal(&code, period, delay, 0.0, fs, 0.05, 42);
            // phasepoints + fc=[0.0] means carrier wipe is identity.
            let peak = correlate_code(
                &sig, &phasepoints(fs, period), &cc, &cc, &[0.0],
                period, period, fs, fft.clone(), ifft.clone(),
            );
            let recovered = peak.best_codephase % period;
            assert_eq!(
                recovered, expected_recovered_sample(period, delay),
                "delay {delay} not recovered correctly"
            );
            assert_eq!(peak.best_freq_idx, 0);
            assert!(peak.best_peak > 0.0);
        }
    }

    #[test]
    fn test_correlate_code_recovers_flipped_second_period() {
        // A sign flip between the two integrated periods (secondary code or
        // data bit, at the code epoch) must be recovered by the `[c, -c]`
        // hypothesis.
        let period = 256usize;
        let code = gold_code(period as f64, 1, 2 * period); // two periods
        let fs = 1.023e6;
        let delay = 50usize;

        let mut planner = FftPlanner::<f32>::new();
        let fft = planner.plan_fft_forward(2 * period);
        let ifft = planner.plan_fft_inverse(2 * period);

        let cc = code_fft_conj(fft.clone(), &code);
        let cc_alt = code_alt_fft_conj(fft.clone(), &code, period);

        // High C/N0 so the recovery is deterministic; flip at the epoch.
        let sig = crate::testutil::synth_signal_cn0(
            &code, period, delay, 0.0, fs, 60.0, 42, true,
        );

        let peak = correlate_code(
            &sig, &phasepoints(fs, 2 * period), &cc, &cc_alt, &[0.0],
            2 * period, period, fs, fft, ifft,
        );
        assert_eq!(
            peak.best_codephase % period,
            expected_recovered_sample(period, delay),
            "flipped signal delay not recovered"
        );
        assert!(peak.best_peak > 0.0);
    }

    #[test]
    fn test_correlate_code_noise_only_bounded() {
        let period = 256usize;
        let code = gold_code(period as f64, 1, period);
        let fs = 1.023e6;
        let mut planner = FftPlanner::<f32>::new();
        let fft = planner.plan_fft_forward(period);
        let ifft = planner.plan_fft_inverse(period);
        let cc = code_fft_conj(fft.clone(), &code);

        let sig = synth_signal(&code, period, 0, 0.0, fs, 0.0, 7);
        let peak = correlate_code(
            &sig, &phasepoints(fs, period), &cc, &cc, &[0.0],
            period, period, fs, fft, ifft,
        );
        assert!(peak.best_peak > 0.0);
        assert!(peak.best_peak.is_finite());
        assert!(peak.second_peak > 0.0);
        assert!(peak.best_codephase < period);
    }
}
