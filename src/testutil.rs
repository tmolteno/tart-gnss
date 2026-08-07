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
}
