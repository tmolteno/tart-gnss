// Copyright (c) 2026 Tim Molteno <tim@elec.ac.nz>
// SPDX-License-Identifier: GPL-3.0

//! Shared helpers for unit tests. Compiled only under `cargo test`
//! (registered as `#[cfg(test)] mod testutil;` in main.rs).

use std::f64::consts::TAU;

/// Phase-accumulator points used by the acquisition pipeline: `TAU * i / fs`.
pub fn phasepoints(fs: f64, n: usize) -> Vec<f64> {
    let pc = TAU / fs;
    (0..n).map(|i| pc * i as f64).collect()
}

/// Build a synthetic signal from a known `code` (length == `num_samples`),
/// delayed circularly by `delay` samples relative to a `code_period`, modulated
/// onto a real carrier at `f0` Hz, plus weak deterministic noise.
///
/// Returns the de-meaned f32 signal the acquisition functions expect
/// (de-meaned because the code has zero mean).
pub fn synth_signal(
    code: &[f64],
    code_period: usize,
    delay: usize,
    f0: f64,
    fs: f64,
    noise_amp: f64,
    seed: u64,
) -> Vec<f32> {
    let n = code.len();
    let mut rng = XorShift(seed);
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let cp = (i + delay) % code_period;
        let carrier = (TAU * f0 * i as f64 / fs).cos();
        let noise = noise_amp * (2.0 * rng.next_f64() - 1.0);
        out.push((code[cp] * carrier + noise) as f32);
    }
    out
}

/// Expected recovered code-phase sample (mod one period) given the FFT
/// correlation convention `signal[i] = code[(i+delay) % period]` →
/// `best_codephase % period == (period - delay) % period`.
pub fn expected_recovered_sample(period: usize, delay: usize) -> usize {
    (period - delay) % period
}

/// Build a synthetic signal like [`synth_signal`], but with white Gaussian
/// noise whose level yields a known C/N0 (dB-Hz).
///
/// Convention (Ma et al. 2024, eq. 7): C/N0 = A²·fs/(4σ²) with unit noise
/// variance σ² = 1, so the carrier amplitude is A = 2·sqrt(C/N0/fs).
/// The signal is `code` (two code periods of `period` samples), circularly
/// delayed by `delay` samples, on a real carrier at `f0` Hz.  With `flip`,
/// the sign flips at the code epoch between the two periods — the physical
/// secondary-code/data-bit flip, which lands at absolute sample
/// `period + (period - delay) % period` in the window.
#[allow(clippy::too_many_arguments)]
pub fn synth_signal_cn0(
    code: &[f64],
    period: usize,
    delay: usize,
    f0: f64,
    fs: f64,
    cn0_db: f64,
    seed: u64,
    flip: bool,
) -> Vec<f32> {
    let a = 2.0 * (10f64.powf(cn0_db / 10.0) / fs).sqrt();
    let epoch = period + (period - delay) % period; // absolute flip position
    let mut rng = XorShift(seed);
    (0..code.len())
        .map(|i| {
            let cp = (i + delay) % period;
            let sign = if flip && i >= epoch { -1.0 } else { 1.0 };
            let carrier = a * (TAU * f0 * i as f64 / fs).cos();
            (sign * code[cp] * carrier + rng.gaussian()) as f32
        })
        .collect()
}

/// Deterministic xorshift64 PRNG (independent of the crate's other PRNGs).
struct XorShift(u64);
impl XorShift {
    fn next_f64(&mut self) -> f64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        (x >> 11) as f64 / (1u64 << 53) as f64
    }
    /// Standard normal via Box-Muller.
    fn gaussian(&mut self) -> f64 {
        let u1 = self.next_f64().max(1e-12);
        let u2 = self.next_f64();
        (-2.0 * u1.ln()).sqrt() * (TAU * u2).cos()
    }
}
