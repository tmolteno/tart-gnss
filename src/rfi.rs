// Copyright (c) 2026 Tim Molteno <tim@elec.ac.nz>
// SPDX-License-Identifier: GPL-3.0

//! Narrowband RFI report per antenna (`--rfi`).
//!
//! TART observations are one-bit quantized at 16.368 MHz, which makes them
//! vulnerable to strong narrowband interferers: a periodic signal that
//! dominates the sign decision (capture effect) both suppresses weak GNSS
//! signals and inflates the acquisition's correlation floor, pinning ACR
//! C/N0 estimates at the table floor (e.g. a ~2.4 MHz clock harmonic on
//! the Dunedin array).  This module scans each antenna's spectrum for
//! narrowband lines and reports an autocorrelation fingerprint so such
//! interferers can be identified before trusting C/N0 results.
//!
//! For each channel the report gives:
//! - the strongest narrowband spectral lines (frequency and dB above the
//!   median spectral floor), found as FFT local maxima with a 15 dB excess;
//! - the period in samples implied by the strongest line;
//! - the autocorrelation at lags 1, 3 and 16 — the fingerprint used to
//!   recognise periodic square-wave-like interferers (a ~2.4 MHz harmonic
//!   shows acf(1) ≈ +0.4, acf(3) ≈ −0.7, acf(16) ≈ −0.5);
//! - dead channels (no samples) are flagged rather than reported.

use num_complex::Complex;
use rustfft::FftPlanner;
use serde::Serialize;

use crate::observation::Observation;

/// Minimum spectral excess (dB over the median floor) for a line to be
/// reported as RFI.  The largest local maximum of white ±1 noise over
/// ~2^19 bins sits ~13 dB above the median floor, so 15 dB leaves a safe
/// false-positive margin while still catching modest narrowband lines
/// (a tone's power concentrates in a single FFT bin, so even a line only a
/// few dB above the floor in a 5 kHz average is far above it per bin).
const MIN_EXCESS_DB: f64 = 15.0;
/// Lines closer than this are treated as one interferer (the strongest
/// wins).  A sampled square wave leaks over several kHz around its
/// fundamental, so the cluster must be wider than that lobe structure;
/// genuine harmonics (e.g. 2.4/4.9/7.3 MHz) are far apart.
const CLUSTER_HZ: f64 = 10_000.0;
/// Maximum lines reported per channel.
const MAX_LINES: usize = 3;

/// A narrowband spectral line found on a channel.
#[derive(Debug, Clone, Serialize)]
pub struct RfiLine {
    /// Line centre frequency (Hz).
    pub frequency: f64,
    /// Peak power above the median spectral floor (dB).
    pub excess_db: f64,
}

/// Per-channel RFI analysis.
#[derive(Debug, Clone, Serialize)]
pub struct RfiChannel {
    /// Antenna index.
    pub antenna: usize,
    /// True if the channel carries no samples (e.g. an unconnected input);
    /// all other fields are then empty.
    pub dead: bool,
    /// Autocorrelation at lag 1 (adjacent-sample correlation).
    pub acf_lag1: Option<f64>,
    /// Autocorrelation at lag 3.
    pub acf_lag3: Option<f64>,
    /// Autocorrelation at lag 16 (one C/A chip at 16.368 MHz).
    pub acf_lag16: Option<f64>,
    /// Period in samples implied by the strongest line (fs / f_line).
    pub period_samples: Option<f64>,
    /// Narrowband spectral lines, strongest first.
    pub lines: Vec<RfiLine>,
}

/// JSON report produced by `--rfi`.
#[derive(Debug, Clone, Serialize)]
pub struct RfiReport {
    /// Antenna indices covered, in order (aligned with `channels`).
    pub antenna_numbers: Vec<usize>,
    /// Per-antenna rows, ordered by antenna index.
    pub channels: Vec<RfiChannel>,
}

fn autocorrelation(z: &[f64], lag: usize) -> f64 {
    let n = z.len();
    (0..n - lag).map(|i| z[i] * z[i + lag]).sum::<f64>() / (n - lag) as f64
}

/// Detect narrowband spectral lines in a de-meaned ±1 channel: FFT local
/// maxima whose power exceeds the median floor by `MIN_EXCESS_DB`, clustered
/// into interferers and truncated to `MAX_LINES`.
fn detect_lines(z: &[f64], fs: f64) -> Vec<RfiLine> {
    let n = z.len();
    if n < 8 {
        return Vec::new();
    }
    let mut planner = FftPlanner::<f64>::new();
    let fft = planner.plan_fft_forward(n);
    let mut buf: Vec<Complex<f64>> = z.iter().map(|&v| Complex::new(v, 0.0)).collect();
    fft.process(&mut buf);

    let half = n / 2;
    let bin_hz = fs / n as f64;
    let powers: Vec<f64> = (0..half).map(|k| buf[k].norm_sqr()).collect();

    // Robust floor: median bin power (narrow lines barely move it).
    let mut sorted = powers.clone();
    sorted.sort_by(|a, b| a.total_cmp(b));
    let floor = sorted[half / 2];

    let mut lines: Vec<RfiLine> = Vec::new();
    for k in 1..half - 1 {
        if powers[k] > powers[k - 1] && powers[k] >= powers[k + 1] {
            let excess_db = 10.0 * (powers[k] / floor).log10();
            if excess_db > MIN_EXCESS_DB {
                lines.push(RfiLine {
                    frequency: k as f64 * bin_hz,
                    excess_db,
                });
            }
        }
    }

    // Cluster: keep the strongest line of each interferer.
    lines.sort_by(|a, b| b.excess_db.total_cmp(&a.excess_db));
    let mut clustered: Vec<RfiLine> = Vec::new();
    for line in lines {
        if clustered
            .iter()
            .all(|c| (c.frequency - line.frequency).abs() > CLUSTER_HZ)
        {
            clustered.push(line);
            if clustered.len() == MAX_LINES {
                break;
            }
        }
    }
    clustered
}

/// Analyse one channel's bipolar ±1 samples (`Observation::get_antenna`
/// output), de-meaning exactly as the acquisition pipeline does.
pub fn analyze_channel(antenna: usize, bipolar: &[f64], fs: f64) -> RfiChannel {
    let mean = bipolar.iter().sum::<f64>() / bipolar.len() as f64;
    let z: Vec<f64> = bipolar.iter().map(|&v| v - mean).collect();
    let var = z.iter().map(|v| v * v).sum::<f64>() / z.len() as f64;
    if var < 1e-9 {
        return RfiChannel {
            antenna,
            dead: true,
            acf_lag1: None,
            acf_lag3: None,
            acf_lag16: None,
            period_samples: None,
            lines: Vec::new(),
        };
    }
    let lines = detect_lines(&z, fs);
    let period_samples = lines.first().map(|l| fs / l.frequency);
    RfiChannel {
        antenna,
        dead: false,
        acf_lag1: Some(autocorrelation(&z, 1)),
        acf_lag3: Some(autocorrelation(&z, 3)),
        acf_lag16: Some(autocorrelation(&z, 16)),
        period_samples,
        lines,
    }
}

/// Build the RFI report for every antenna of an observation.
pub fn run(obs: &Observation) -> RfiReport {
    let fs = obs.get_sampling_rate();
    let n_ant = obs.config.num_antenna();
    let channels: Vec<RfiChannel> = (0..n_ant)
        .map(|ant| {
            let bipolar = obs.get_antenna(ant);
            analyze_channel(ant, &bipolar, fs)
        })
        .collect();
    RfiReport {
        antenna_numbers: (0..n_ant).collect(),
        channels,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulator::{gaussian, XorShift64};
    use std::f64::consts::TAU;

    const FS: f64 = 16.368e6;
    const N: usize = 64 * 16_368; // 64 ms

    /// Square wave at ~2.43 MHz (the Dunedin array's clock harmonic).
    fn square_243(i: usize) -> f64 {
        let per = FS / 2.43e6;
        if (i as f64 % per) < per / 2.0 {
            1.0
        } else {
            -1.0
        }
    }

    /// One-bit quantize a continuous signal to bipolar ±1 (production path).
    fn bipolar(sig: &[f64]) -> Vec<f64> {
        crate::simulator::quantize_bipolar_demean(sig)
            .iter()
            .map(|&v| v as f64)
            .collect()
    }

    #[test]
    fn test_square_wave_detected() {
        // Strong 2.43 MHz square wave + noise, one-bit quantized.
        let mut rng = XorShift64::new(42);
        let sig: Vec<f64> = (0..N).map(|i| 2.0 * square_243(i) + gaussian(&mut rng)).collect();
        let ch = analyze_channel(0, &bipolar(&sig), FS);
        assert!(!ch.dead);
        let top = ch.lines.first().expect("expected a detected line");
        assert!(
            (top.frequency - 2.43e6).abs() < 1e3,
            "top line at {} Hz",
            top.frequency
        );
        assert!(top.excess_db > 15.0, "excess {}", top.excess_db);
        let p = ch.period_samples.unwrap();
        assert!((p - FS / 2.43e6).abs() < 0.2, "period {p}");
        // Square-wave fingerprint.
        assert!(ch.acf_lag1.unwrap() > 0.3, "acf1 {}", ch.acf_lag1.unwrap());
        assert!(ch.acf_lag3.unwrap() < -0.4, "acf3 {}", ch.acf_lag3.unwrap());
        assert!(ch.acf_lag16.unwrap() < 0.0, "acf16 {}", ch.acf_lag16.unwrap());
    }

    #[test]
    fn test_white_noise_no_lines() {
        let mut rng = XorShift64::new(7);
        let sig: Vec<f64> = (0..N).map(|_| gaussian(&mut rng)).collect();
        let ch = analyze_channel(0, &bipolar(&sig), FS);
        assert!(!ch.dead);
        assert!(
            ch.lines.is_empty(),
            "white noise flagged as RFI: {:?}",
            ch.lines
        );
    }

    #[test]
    fn test_line_excess_scales_with_amplitude() {
        // The FFT's coherent gain over 64 ms makes even small periodic
        // components detectable, so both weak and strong interferers are
        // reported — the excess quantifies severity.  Same noise seed for
        // both so the amplitude difference is isolated.
        let mut rng = XorShift64::new(9);
        let weak: Vec<f64> = (0..N).map(|i| 0.1 * square_243(i) + gaussian(&mut rng)).collect();
        let mut rng = XorShift64::new(9);
        let strong: Vec<f64> = (0..N).map(|i| 2.0 * square_243(i) + gaussian(&mut rng)).collect();
        let w = analyze_channel(0, &bipolar(&weak), FS);
        let s = analyze_channel(0, &bipolar(&strong), FS);
        let w_line = w.lines.first().expect("weak line not detected");
        let s_line = s.lines.first().expect("strong line not detected");
        assert!(
            s_line.excess_db > w_line.excess_db + 10.0,
            "strong {:.1} dB vs weak {:.1} dB",
            s_line.excess_db,
            w_line.excess_db
        );
    }

    #[test]
    fn test_gps_signal_not_flagged() {
        // A 40 dB-Hz GPS C/A signal is spread-spectrum (~2 MHz wide); its
        // per-bin power is far below the noise floor, so no narrowband line.
        let mut rng = XorShift64::new(11);
        let period = 16_368usize;
        let code = crate::acquisition::gold_code(period as f64, 1, 2 * period);
        let a = 2.0 * (10f64.powf(40.0 / 10.0) / FS).sqrt();
        let delay = 137usize;
        let sig: Vec<f64> = (0..N)
            .map(|i| {
                let cp = (i + delay) % period;
                code[cp] * (TAU * 100e3 * i as f64 / FS).cos() * a + gaussian(&mut rng)
            })
            .collect();
        let ch = analyze_channel(0, &bipolar(&sig), FS);
        assert!(ch.lines.is_empty(), "GPS signal flagged: {:?}", ch.lines);
    }

    #[test]
    fn test_dead_channel() {
        let dead = vec![0.0f64; N];
        let ch = analyze_channel(17, &dead, FS);
        assert!(ch.dead);
        assert!(ch.lines.is_empty());
        assert!(ch.acf_lag1.is_none());
    }
}
