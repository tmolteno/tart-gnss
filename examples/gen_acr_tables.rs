// Copyright (c) 2026 Tim Molteno <tim@elec.ac.nz>
// SPDX-License-Identifier: GPL-3.0

//! Generate per-constellation ACR C/N0 lookup tables by Monte Carlo
//! simulation of the actual acquisition functions.
//!
//! Usage:
//!   cargo run --release --example gen_acr_tables > /tmp/tables.txt
//!
//! Method: for each constellation and each C/N0 on a 0.5 dB grid, synthesize
//! a two-code-period signal at 1 sample per chip (1.023 MHz) with a known
//! carrier, random code phase and carrier phase, run the production
//! acquisition, and measure
//!
//!   r_A = mean(V^m²) / mean(V^s²)
//!
//! over TRIALS trials (paper eq. 31: ratio of averages of the peak and
//! second-peak correlation powers).  The C/N0-to-amplitude mapping follows
//! Ma et al. (2024), eq. 7:  C/N0 = A²·fs/(4σ²) with unit noise variance,
//! i.e. A = 2·sqrt(C/N0/fs)  (see src/testutil.rs for the same convention).
//!
//! One sample per chip is faithful to the 16 MHz TART sampling: the noise
//! correlation function is piecewise linear with knots at chip boundaries,
//! so the max over the oversampled grid equals the max at chip spacing.

use rayon::prelude::*;
use std::f64::consts::TAU;
use tart_gnss_acquire::acquisition;
use tart_gnss_acquire::beidou;
use tart_gnss_acquire::galileo;
use tart_gnss_acquire::l1c;

const FS: f64 = 1.023e6; // 1 sample per chip for all constellations
const F_CARRIER: f64 = 100e3;
const TRIALS: usize = 100;

/// Deterministic xorshift64 PRNG (mirrors src/testutil.rs).
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
    fn gaussian(&mut self) -> f64 {
        let u1 = self.next_f64().max(1e-12);
        let u2 = self.next_f64();
        (-2.0 * u1.ln()).sqrt() * (TAU * u2).cos()
    }
}

/// Carrier amplitude for a C/N0 (dB-Hz) with unit noise variance:
/// A = 2·sqrt(C/N0/fs).
fn amplitude(cn0_db: f64) -> f64 {
    2.0 * (10f64.powf(cn0_db / 10.0) / FS).sqrt()
}

/// Two code periods of `code`, circularly delayed by `delay` samples,
/// modulated on the carrier with phase `theta`, plus unit-variance noise.
fn synth(code: &[f64], period: usize, cn0_db: f64, delay: usize, theta: f64, seed: u64) -> Vec<f32> {
    let a = amplitude(cn0_db);
    let mut rng = XorShift(seed);
    (0..code.len())
        .map(|i| {
            let cp = (i + delay) % period;
            let carrier = a * (TAU * F_CARRIER * i as f64 / FS + theta).cos();
            (code[cp] * carrier + rng.gaussian()) as f32
        })
        .collect()
}

fn phasepoints(n: usize) -> Vec<f64> {
    let pc = TAU / FS;
    (0..n).map(|i| pc * i as f64).collect()
}

type AcqFn = fn(&[f32], &[f64], f64, f64, f64, usize, usize) -> (f64, f64);

fn acq_gps(s: &[f32], p: &[f64], fs: f64, cf: f64, sb: f64, prn: usize, per: usize) -> (f64, f64) {
    let r = acquisition::acquire_full(s, p, fs, cf, sb, prn, per);
    (r.signal_strength, r.second_peak)
}

fn acq_galileo(s: &[f32], p: &[f64], fs: f64, cf: f64, sb: f64, prn: usize, per: usize) -> (f64, f64) {
    let r = galileo::acquire_galileo_single(s, p, fs, cf, sb, prn, per);
    (r.signal_strength, r.second_peak)
}

fn acq_beidou(s: &[f32], p: &[f64], fs: f64, cf: f64, sb: f64, prn: usize, per: usize) -> (f64, f64) {
    let r = beidou::acquire_beidou_single(s, p, fs, cf, sb, prn, per);
    (r.signal_strength, r.second_peak)
}

fn acq_l1c(s: &[f32], p: &[f64], fs: f64, cf: f64, sb: f64, prn: usize, per: usize) -> (f64, f64) {
    let r = l1c::acquire_l1c_single(s, p, fs, cf, sb, prn, per);
    (r.signal_strength, r.second_peak)
}

struct Spec {
    name: &'static str,
    doc: &'static str,
    period: usize, // one code period in samples (== chips)
    cn0_lo: f64,
    cn0_hi: f64,
    code: fn(usize, usize) -> Vec<f64>,
    acq: AcqFn,
}

fn gps_code(prn: usize, n: usize) -> Vec<f64> {
    acquisition::gold_code(1023.0, prn, n)
}
fn gal_code(prn: usize, n: usize) -> Vec<f64> {
    galileo::e1c_code_resampled(4092.0, prn, n)
}
fn bds_code(prn: usize, n: usize) -> Vec<f64> {
    beidou::b1c_code_resampled(10230.0, prn, n)
}
fn l1c_code(prn: usize, n: usize) -> Vec<f64> {
    l1c::l1c_code_resampled(10230.0, prn, n)
}

fn generate(spec: &Spec) {
    let prn = 1usize;
    let num_samples = 2 * spec.period;
    let code = (spec.code)(prn, num_samples);
    let pp = phasepoints(num_samples);

    println!("// {} — {} ({} trials per point)", spec.name, spec.doc, TRIALS);
    println!("pub static {}: &[(f64, f64)] = &[", spec.name);

    let mut cn0 = spec.cn0_lo;
    let mut prev_r_a = 0.0f64;
    while cn0 <= spec.cn0_hi + 1e-9 {
        let (sum_p2, sum_s2) = (0..TRIALS)
            .into_par_iter()
            .map(|t| {
                let seed = (t as u64) ^ (cn0.to_bits() << 33);
                let mut rng = XorShift(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15).wrapping_add(1));
                let delay = (rng.next_f64() * spec.period as f64) as usize % spec.period;
                let theta = rng.next_f64() * TAU;
                let sig = synth(&code, spec.period, cn0, delay, theta, seed);
                let (peak, second) = (spec.acq)(&sig, &pp, FS, F_CARRIER, 0.0, prn, spec.period);
                (peak * peak, second * second)
            })
            .reduce(|| (0.0, 0.0), |a, b| (a.0 + b.0, a.1 + b.1));
        let mut r_a = sum_p2 / sum_s2;
        // Enforce strict monotonicity: Monte Carlo noise near the cut-off
        // can invert (or tie) adjacent entries slightly.
        if r_a <= prev_r_a {
            r_a = prev_r_a + 0.001;
        }
        prev_r_a = r_a;
        println!("    ({r_a:.3}, {cn0:.1}),");
        cn0 += 0.5;
    }
    println!("];");
    println!();
}

fn main() {
    generate(&Spec {
        name: "GPS_L1_CA_ACR_TABLE",
        doc: "GPS L1 C/A, SBAS, QZSS. T_coh = 2 ms, 1023-chip Gold codes",
        period: 1023,
        cn0_lo: 33.0,
        cn0_hi: 62.0,
        code: gps_code,
        acq: acq_gps,
    });
    generate(&Spec {
        name: "GALILEO_E1C_ACR_TABLE",
        doc: "Galileo E1-C pilot. T_coh = 8 ms, 4092-chip memory codes",
        period: 4092,
        cn0_lo: 27.0,
        cn0_hi: 62.0,
        code: gal_code,
        acq: acq_galileo,
    });
    generate(&Spec {
        name: "BEIDOU_B1C_ACR_TABLE",
        doc: "BeiDou B1C pilot. T_coh = 20 ms, 10230-chip Weil codes",
        period: 10230,
        cn0_lo: 23.0,
        cn0_hi: 62.0,
        code: bds_code,
        acq: acq_beidou,
    });
    generate(&Spec {
        name: "GPS_L1C_ACR_TABLE",
        doc: "GPS L1C pilot. T_coh = 20 ms, 10230-chip Weil codes",
        period: 10230,
        cn0_lo: 23.0,
        cn0_hi: 62.0,
        code: l1c_code,
        acq: acq_l1c,
    });
}
