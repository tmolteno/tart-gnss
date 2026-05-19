//! GPS L1 C/A code generation and signal acquisition via FFT-based
//! circular cross-correlation.
//!
//! Ported from `tart/tart/operation/acquisition.py`.

use num_complex::Complex;
use rustfft::FftPlanner;
use std::f64::consts::TAU;

/// Result of a GPS signal acquisition attempt for a single antenna.
#[derive(Debug, Clone, serde::Serialize)]
pub struct AcquisitionResult {
    pub prn: usize,
    /// Peak correlation magnitude (normalised).
    pub signal_strength: f64,
    /// Code-phase offset in fractions of a millisecond [0, 1).
    pub codephase_frac: f64,
    /// Doppler frequency offset in Hz (relative to centre frequency).
    pub frequency: f64,
}

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
const CA_CHIPS: usize = 1023;

/// Generate the 1023-chip GPS C/A gold code for a given PRN (1-based).
///
/// Returns a `[f64; CA_CHIPS]` array of ±1 values.
pub fn generate_ca_code(prn: usize) -> [f64; CA_CHIPS] {
    assert!(prn >= 1 && prn <= 38, "PRN must be in 1..=38, got {prn}");
    let g2shift = CODE_DELAY_TABLE[prn - 1];

    // --- G1 ----------------------------------------------------------------
    let mut g1 = [0.0f64; CA_CHIPS];
    let mut lfsr = [-1.0f64; 10]; // initialised to all -1 (bipolar)
    for i in 0..CA_CHIPS {
        g1[i] = lfsr[9];
        let save_bit = lfsr[2] * lfsr[9];
        // shift right: lfsr[0..9] → lfsr[1..10], insert save_bit at 0
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

/// Resample the C/A code to `samples_per_code` samples per period,
/// repeating for `epochs` full periods.
pub fn gold_code(samples_per_code: f64, prn: usize, epochs: f64) -> Vec<f64> {
    let ca = generate_ca_code(prn);
    let samples_per_chip = samples_per_code / CA_CHIPS as f64;
    let num_samples = (samples_per_code * epochs).floor() as usize;

    (0..num_samples)
        .map(|n| {
            let idx = ((n as f64 / samples_per_chip).floor() as usize) % CA_CHIPS;
            ca[idx]
        })
        .collect()
}

// ---------------------------------------------------------------------------
// FFT-based acquisition
// ---------------------------------------------------------------------------

/// Perform parallel code-phase search (FFT circular cross-correlation) for a
/// single GPS PRN over a frequency search band.
///
/// Returns `[PRN, signal_strength, codephase_frac, frequency]`.
pub fn acquire_full(
    x: &[f64],
    sampling_freq: f64,
    center_freq: f64,
    search_band: f64,
    prn: usize,
) -> AcquisitionResult {
    let sampling_period = 1.0 / sampling_freq;
    let samples_per_ms = sampling_freq / 1000.0;
    let samples_per_chunk = samples_per_ms as usize;

    let epochs_available = (x.len() as f64 / samples_per_ms).floor();
    let total_samples = (epochs_available * samples_per_chunk as f64) as usize;

    // --- Frequency bins ----------------------------------------------------
    let freq_bin_size: f64 = 300.0;
    let n_freq_bins =
        (2.0 * search_band / freq_bin_size).round() as usize + 1;
    let fc: Vec<f64> = (0..n_freq_bins)
        .map(|i| center_freq - search_band + (i as f64) * (2.0 * search_band) / (n_freq_bins as f64 - 1.0))
        .collect();

    // --- FFT planner -------------------------------------------------------
    let mut planner: FftPlanner<f32> = FftPlanner::new();
    let fft = planner.plan_fft_forward(total_samples);
    let ifft = planner.plan_fft_inverse(total_samples);

    // --- Local code replica ------------------------------------------------
    let code = gold_code(samples_per_ms, prn, epochs_available);
    let mut code_complex: Vec<Complex<f32>> = code
        .iter()
        .map(|&v| Complex::new(v as f32, 0.0))
        .collect();
    fft.process(&mut code_complex);
    // conjugate in-place
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
        // exp(j * phasepoints * freq) * signal
        let mut iq: Vec<Complex<f32>> = phasepoints
            .iter()
            .zip(signal_f32.iter())
            .map(|(&p, &s)| {
                let phase = (p * freq) as f32;
                Complex::new(phase.cos(), phase.sin()) * Complex::new(s, 0.0)
            })
            .collect();

        fft.process(&mut iq);

        // Multiply by conjugated code spectrum
        for (iq_elem, &code_elem) in iq.iter_mut().zip(code_complex.iter()) {
            *iq_elem *= code_elem;
        }

        ifft.process(&mut iq);

        // Magnitude / sqrt(N)
        let scale = 1.0 / (signal_len as f32).sqrt();
        let corr: Vec<f32> = iq.iter().map(|c| c.norm() * scale).collect();

        // Find peak in this frequency bin
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

    AcquisitionResult {
        prn,
        signal_strength: best_peak as f64,
        codephase_frac,
        frequency,
    }
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
