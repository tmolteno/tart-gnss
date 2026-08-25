// Copyright (c) 2026 Tim Molteno <tim@elec.ac.nz>
// SPDX-License-Identifier: GPL-3.0

//! Minimal TART raw-data simulation.
//!
//! Synthesises per-antenna 0/1 samples from catalogue sources
//! (`{name, az, el, r, jy}`) with correct geometric phase delays, and
//! writes an HDF5 observation file readable by `Observation::from_hdf5`.

use serde::Deserialize;
use hdf5_writer::{DatasetBuilder, Hdf5Builder, Hdf5Writer, WriteOptions};

/// Speed of light (m/s), matching TART's `constants.V_LIGHT`.
pub const V_LIGHT: f64 = 2.99793e8;
/// Sources beyond this range are treated as plane waves (TART convention).
pub const PLANE_WAVE_RANGE: f64 = 1.0e4;
/// Fixed default RNG seed used when `--seed` is not supplied.
pub const DEFAULT_SEED: u64 = 0x5EED_2026;

/// A single catalogue source in local horizontal coordinates.
#[derive(Debug, Clone, Deserialize)]
pub struct Source {
    pub name: String,
    /// Azimuth in degrees (0 = North, increasing toward East).
    pub az: f64,
    /// Elevation in degrees.
    pub el: f64,
    /// Range in metres.
    pub r: f64,
    /// Flux density in Janskys.
    pub jy: f64,
}

/// Simulation parameters.
#[derive(Debug, Clone)]
pub struct SimConfig {
    pub sample_rate: f64,
    pub center_freq: f64,
    pub band: f64,
    pub samples: usize,
    pub gain: f64,
    pub seed: u64,
}

/// Parse a list of sources from a catalogue JSON file
/// (`[{"name": ..., "az": ..., "el": ..., "r": ..., "jy": ...}, ...]`).
pub fn parse_sources(json: &str) -> Result<Vec<Source>, serde_json::Error> {
    serde_json::from_str(json)
}

/// Parse a list of antenna positions `[[east, north, up], ...]` (metres).
/// The up (z) coordinate is forced to 0.
///
/// Returns an error (rather than panicking) if a position row has fewer than
/// two coordinates (`[east, north]`) or is empty.
pub fn parse_positions(json: &str) -> Result<Vec<[f64; 3]>, serde_json::Error> {
    use serde::de::Error as _;

    let raw: Vec<Vec<f64>> = serde_json::from_str(json)?;
    let mut out = Vec::with_capacity(raw.len());
    for v in raw {
        if v.len() < 2 {
            return Err(serde_json::Error::custom(
                "each antenna position must have at least 2 coordinates (east, north)",
            ));
        }
        out.push([v[0], v[1], 0.0]);
    }
    Ok(out)
}

/// Deterministic xorshift64 PRNG (samples layout + noise).
pub struct XorShift64(u64);
impl XorShift64 {
    pub fn new(seed: u64) -> Self {
        Self(seed | 1)
    }
    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }
    /// Uniform in [0,1).
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Standard-normal sample via Box-Muller from a seeded PRNG.
pub fn gaussian(rng: &mut XorShift64) -> f64 {
    let u1 = rng.next_f64().max(1e-12);
    let u2 = rng.next_f64();
    (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
}

/// Generate `n` antenna positions uniformly in a circle of radius
/// `diameter/2` in the East-North plane (z = 0), deterministically seeded.
pub fn random_positions(n: usize, diameter: f64, seed: u64) -> Vec<[f64; 3]> {
    let mut rng = XorShift64::new(seed);
    let r = diameter / 2.0;
    (0..n)
        .map(|_| {
            // Uniform in disc via sqrt on radius.
            let rr = r * rng.next_f64().sqrt();
            let theta = std::f64::consts::TAU * rng.next_f64();
            [rr * theta.cos(), rr * theta.sin(), 0.0]
        })
        .collect()
}

/// Unit direction vector toward a source: `[E, N, U]`.
/// Azimuth 0° = North, increasing East; elevation above the horizon.
pub fn source_direction(az_deg: f64, el_deg: f64) -> [f64; 3] {
    let az = az_deg.to_radians();
    let el = el_deg.to_radians();
    [
        az.sin() * el.cos(),
        az.cos() * el.cos(),
        el.sin(),
    ]
}

/// Geometric delay (seconds) from a source to an antenna position (ENU, m).
///
/// `Δ = (|r·ŝ − pos| − r)/c`. Negative means the wavefront reaches the antenna
/// before the array reference (origin). Sources beyond `PLANE_WAVE_RANGE` (1e4 m)
/// are treated as plane waves, matching TART.
pub fn geometric_delay(az_deg: f64, el_deg: f64, r: f64, pos: [f64; 3]) -> f64 {
    let rr = if r > PLANE_WAVE_RANGE { PLANE_WAVE_RANGE } else { r };
    let s = source_direction(az_deg, el_deg);
    let src = [rr * s[0], rr * s[1], rr * s[2]];
    let dx = src[0] - pos[0];
    let dy = src[1] - pos[1];
    let dz = src[2] - pos[2];
    let dist = (dx * dx + dy * dy + dz * dz).sqrt();
    (dist - rr) / V_LIGHT
}

/// Assign one distinct frequency (Hz) per source, evenly spread across
/// `[center - band/2, center + band/2]`. A single source sits at `center`.
pub fn assign_frequencies(center: f64, band: f64, n: usize) -> Vec<f64> {
    if n <= 1 {
        return vec![center];
    }
    (0..n)
        .map(|k| center + band * (k as f64 / (n as f64 - 1.0) - 0.5))
        .collect()
}

/// Tone amplitude proportional to the square root of the flux (Jy).
pub fn amplitude(jy: f64, gain: f64) -> f64 {
    gain * jy.sqrt()
}

/// Per-antenna Gaussian noise standard deviation from an SNR in dB:
/// `σ = sqrt( signal_power / 10^(snr/10) )`, `signal_power = Σ A²/2`.
pub fn noise_std(amplitudes: &[f64], snr_db: f64) -> f64 {
    let signal_power: f64 = amplitudes.iter().map(|a| a * a / 2.0).sum();
    (signal_power / 10f64.powf(snr_db / 10.0)).sqrt()
}

/// Per-antenna continuous-valued (pre-quantization) signals with a
/// per-antenna gain multiplier applied to the signal component.
///
/// For antenna `j`, `signals[j][n] = Σ_k amps[k]·gains[j]·cos(2π f_k (t_n + Δ_{k,j})) + noise`,
/// with `t_n = n / sample_rate`.  The noise (`sigma`) is computed from the
/// *unscaled* amplitudes, so gains[j] shifts that antenna's effective C/N0 by
/// `20·log10(gains[j])` dB — the lever used by the antenna-quality tests.
pub fn antenna_signals_with_gains(
    sources: &[Source],
    positions: &[[f64; 3]],
    cfg: &SimConfig,
    snr_db: f64,
    gains: &[f64],
) -> Vec<Vec<f64>> {
    let n_src = sources.len();
    let freqs = assign_frequencies(cfg.center_freq, cfg.band, n_src);
    let amps: Vec<f64> = sources.iter().map(|s| amplitude(s.jy, cfg.gain)).collect();

    // Delay per (source, antenna).
    let mut delays = vec![vec![0.0f64; positions.len()]; n_src];
    for (k, src) in sources.iter().enumerate() {
        for (j, pos) in positions.iter().enumerate() {
            delays[k][j] = geometric_delay(src.az, src.el, src.r, *pos);
        }
    }

    let sigma = noise_std(&amps, snr_db);
    let mut rng = XorShift64::new(cfg.seed);
    let two_pi = std::f64::consts::TAU;

    (0..positions.len())
        .map(|j| {
            let g = gains[j.min(gains.len() - 1)];
            let mut sig = Vec::with_capacity(cfg.samples);
            for i in 0..cfg.samples {
                let t = i as f64 / cfg.sample_rate;
                let mut v = 0.0f64;
                for k in 0..n_src {
                    v += amps[k] * g * (two_pi * freqs[k] * (t + delays[k][j])).cos();
                }
                v += gaussian(&mut rng) * sigma;
                sig.push(v);
            }
            sig
        })
        .collect()
}

/// Per-antenna continuous-valued (pre-quantization) signals, all gains 1.0.
///
/// For antenna `j`, `signals[j][n] = Σ_k A_k·cos(2π f_k (t_n + Δ_{k,j})) + noise`,
/// with `t_n = n / sample_rate`.
pub fn antenna_signals(
    sources: &[Source],
    positions: &[[f64; 3]],
    cfg: &SimConfig,
    snr_db: f64,
) -> Vec<Vec<f64>> {
    let ones = vec![1.0; positions.len()];
    antenna_signals_with_gains(sources, positions, cfg, snr_db, &ones)
}

/// One-bit quantize a continuous signal to 0/1 unipolar (NRZ; zero -> 1).
pub fn quantize(signal: &[f64]) -> Vec<u8> {
    signal.iter().map(|&v| if v >= 0.0 { 1 } else { 0 }).collect()
}

/// One-bit quantize a continuous signal to bipolar ±1 (0 -> -1, 1 -> +1,
/// matching `Observation::get_antenna`), then remove the mean exactly as
/// the acquisition pipeline does before correlating.  This is the data
/// path every TART observation takes (1-bit samples), so the ACR lookup
/// tables are calibrated against it (`examples/gen_acr_tables.rs`).
pub fn quantize_bipolar_demean(signal: &[f64]) -> Vec<f32> {
    let bipolar: Vec<f64> = signal
        .iter()
        .map(|&v| if v >= 0.0 { 1.0 } else { -1.0 })
        .collect();
    let mean = bipolar.iter().sum::<f64>() / bipolar.len() as f64;
    bipolar.iter().map(|&v| (v - mean) as f32).collect()
}

/// Produce raw 0/1 unipolar samples per antenna, scaling each antenna's
/// signal component by the corresponding entry in `gains`.
///
/// Primarily used by tests (e.g. the `--test-antennas` end-to-end ranking
/// test), so it is allowed to be dead code in production builds.
#[allow(dead_code)]
pub fn synthesize_with_gains(
    sources: &[Source],
    positions: &[[f64; 3]],
    cfg: &SimConfig,
    snr_db: f64,
    gains: &[f64],
) -> Vec<Vec<u8>> {
    antenna_signals_with_gains(sources, positions, cfg, snr_db, gains)
        .into_iter()
        .map(|s| quantize(&s))
        .collect()
}

/// Produce raw 0/1 unipolar samples per antenna (all gains 1.0).
pub fn synthesize(
    sources: &[Source],
    positions: &[[f64; 3]],
    cfg: &SimConfig,
    snr_db: f64,
) -> Vec<Vec<u8>> {
    antenna_signals(sources, positions, cfg, snr_db)
        .into_iter()
        .map(|s| quantize(&s))
        .collect()
}

/// Pack 0/1 bits MSB-first into bytes (matching `observation::unpack_bits` /
/// `numpy.packbits`). Trailing bits beyond a multiple of 8 are zero-padded.
pub fn pack_row(bits: &[u8]) -> Vec<u8> {
    bits.chunks(8)
        .map(|chunk| {
            let mut byte = 0u8;
            for (b, bit) in chunk.iter().enumerate() {
                if *bit == 1 {
                    byte |= 1 << (7 - b);
                }
            }
            byte
        })
        .collect()
}

/// Write a raw observation to HDF5 in the layout `Observation::from_hdf5`
/// expects: `config` (JSON string), `timestamp` (ISO-8601 string), and `data`
/// (2-D u8 of packed 1-bit samples, shape `[num_antenna, row_bytes]`).
pub fn write_observation(
    path: &str,
    num_antenna: usize,
    sampling_frequency: f64,
    timestamp: &str,
    data: &[Vec<u8>],
) -> Result<(), Box<dyn std::error::Error>> {
    let row_bytes = data.first().map(|r| r.len()).unwrap_or(0);
    let config_json =
        format!(r#"{{"num_antenna":{num_antenna},"sampling_frequency":{sampling_frequency}}}"#);
    let flat: Vec<u8> = data.iter().flatten().copied().collect();

    let config_ds = DatasetBuilder::fixed_string_data("config", [], &[config_json.as_str()])?;
    let ts_ds = DatasetBuilder::fixed_string_data("timestamp", [], &[timestamp])?;
    let data_ds = DatasetBuilder::typed_data::<u8>(
        "data",
        [num_antenna as u64, row_bytes as u64],
        &flat,
    )?;

    let plan = Hdf5Builder::new()
        .dataset(config_ds)
        .dataset(ts_ds)
        .dataset(data_ds)
        .into_plan()?;

    let file = std::fs::File::create(path)?;
    let writer = Hdf5Writer::new(file, WriteOptions::default());
    writer.finish(plan)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_sources() {
        let src = parse_sources(
            r#"[{"name":"sun","az":10.0,"el":30.0,"r":1e10,"jy":10000.0}]"#,
        )
        .unwrap();
        assert_eq!(src.len(), 1);
        assert_eq!(src[0].name, "sun");
        assert_eq!(src[0].az, 10.0);
        assert_eq!(src[0].jy, 10000.0);
    }

    #[test]
    fn test_parse_positions_forces_z_zero() {
        let p = parse_positions(r#"[[1.0,2.0,5.0],[3.0,4.0]]"#).unwrap();
        assert_eq!(p, vec![[1.0, 2.0, 0.0], [3.0, 4.0, 0.0]]);
    }

    #[test]
    fn test_parse_positions_rejects_short_row() {
        // A row with fewer than 2 coordinates must error, not panic.
        assert!(parse_positions(r#"[[1.0]]"#).is_err());
        assert!(parse_positions(r#"[[]]"#).is_err());
        assert!(parse_positions(r#"[[1.0,2.0],[5.0]]"#).is_err());
        // Valid rows still parse.
        assert!(parse_positions(r#"[[1.0,2.0]]"#).is_ok());
    }

    #[test]
    fn test_random_positions_deterministic_and_bounded() {
        let n = 24;
        let d = 3.0;
        let a = random_positions(n, d, 42);
        let b = random_positions(n, d, 42);
        assert_eq!(a, b);
        for pos in &a {
            let dist = (pos[0] * pos[0] + pos[1] * pos[1]).sqrt();
            assert!(dist <= d / 2.0 + 1e-9);
            assert_eq!(pos[2], 0.0);
        }
    }

    #[test]
    fn test_source_direction_north_horizontal() {
        let d = source_direction(0.0, 0.0);
        assert!(d[0].abs() < 1e-12, "east={}", d[0]);
        assert!((d[1] - 1.0).abs() < 1e-12, "north={}", d[1]);
        assert!(d[2].abs() < 1e-12, "up={}", d[2]);
    }

    #[test]
    fn test_source_direction_zenith() {
        let d = source_direction(30.0, 90.0);
        assert!((d[2] - 1.0).abs() < 1e-12, "up={}", d[2]);
        assert!(d[0].abs() < 1e-12 && d[1].abs() < 1e-12);
    }

    #[test]
    fn test_geometric_delay_plane_wave() {
        let d = geometric_delay(0.0, 90.0, 1e7, [0.0, 0.0, 0.0]);
        assert!(d.abs() < 1e-12);

        let d0 = geometric_delay(0.0, 45.0, 1e7, [0.0, 0.0, 0.0]);
        let d1 = geometric_delay(0.0, 45.0, 1e7, [0.0, 0.5, 0.0]);
        let diff = (d1 - d0) * V_LIGHT;
        let expect = -0.5 * 45.0f64.to_radians().cos();
        // Tolerance accounts for near-field curvature from the 1e4 m plane-wave cap.
        assert!((diff - expect).abs() < 1e-4, "diff={diff} expect={expect}");
    }

    #[test]
    fn test_geometric_delay_plane_cap() {
        let d_far = geometric_delay(0.0, 45.0, 1e10, [0.0, 0.5, 0.0]);
        let d_near = geometric_delay(0.0, 45.0, 1e4, [0.0, 0.5, 0.0]);
        assert!((d_far - d_near).abs() < 1e-12);
    }

    #[test]
    fn test_assign_frequencies_single() {
        assert_eq!(assign_frequencies(4.092e6, 2.0e6, 1), vec![4.092e6]);
    }

    #[test]
    fn test_assign_frequencies_spread() {
        let fs = assign_frequencies(4.0e6, 2.0e6, 3);
        assert_eq!(fs, vec![3.0e6, 4.0e6, 5.0e6]);
    }

    #[test]
    fn test_amplitude_sqrt_flux() {
        assert!((amplitude(4.0, 1.0) - 2.0).abs() < 1e-12);
        assert!((amplitude(4.0, 2.0) - 4.0).abs() < 1e-12);
        // 2x flux -> sqrt(2) amplitude.
        let r = amplitude(4.0, 1.0) / amplitude(2.0, 1.0);
        assert!((r - 2.0f64.sqrt()).abs() < 1e-9);
    }

    #[test]
    fn test_noise_std_snr() {
        assert!((noise_std(&[2.0], 0.0) - 2.0f64.sqrt()).abs() < 1e-9);
        assert!((noise_std(&[2.0], 10.0) - 0.2f64.sqrt()).abs() < 1e-9);
    }

    use std::f64::consts::TAU;

    fn make_config() -> SimConfig {
        SimConfig {
            sample_rate: 16.368e6,
            center_freq: 4.092e6,
            band: 2.0e6,
            samples: 8192,
            gain: 1.0,
            seed: 123,
        }
    }

    #[test]
    fn test_quantize() {
        assert_eq!(quantize(&[1.0, -0.5, 0.0, 3.0]), vec![1, 0, 1, 1]);
    }

    #[test]
    fn test_pack_row_msb_first() {
        assert_eq!(pack_row(&[1, 0, 1, 0, 1, 0, 1, 0]), vec![0b10101010]);
        let packed = pack_row(&[1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1]);
        assert_eq!(packed, vec![0b11110000, 0b01010000]);
    }

    #[test]
    fn test_phase_delay_recovery() {
        let src = Source { name: "t".into(), az: 30.0, el: 60.0, r: 1e7, jy: 1.0 };
        let pos0 = [0.0, 0.0, 0.0];
        let pos1 = [0.6, 0.2, 0.0];
        let cfg = SimConfig { samples: 8192, ..make_config() };
        let fs = cfg.sample_rate;
        let f = cfg.center_freq;

        let sigs = antenna_signals(std::slice::from_ref(&src), &[pos0, pos1], &cfg, 100.0);
        let d0 = geometric_delay(src.az, src.el, src.r, pos0);
        let d1 = geometric_delay(src.az, src.el, src.r, pos1);

        let dft_phase = |sig: &[f64]| -> f64 {
            let (mut re, mut im) = (0.0f64, 0.0f64);
            for (i, &sval) in sig.iter().enumerate() {
                let (s, c) = (TAU * f * i as f64 / fs).sin_cos();
                re += sval * c;
                im -= sval * s;
            }
            im.atan2(re)
        };
        let measured = (dft_phase(&sigs[1]) - dft_phase(&sigs[0])).rem_euclid(TAU);
        let expected = (TAU * f * (d1 - d0)).rem_euclid(TAU);
        assert!(
            (measured - expected).abs() < 1e-2,
            "measured {measured} expected {expected}"
        );
    }

    #[test]
    fn test_synthesize_dimensions_and_unipolar() {
        let src = Source { name: "a".into(), az: 10.0, el: 20.0, r: 5e3, jy: 100.0 };
        let cfg = make_config();
        let positions = random_positions(4, 3.0, 7);
        let data = synthesize(&[src], &positions, &cfg, 20.0);
        assert_eq!(data.len(), 4);
        for ant in &data {
            assert_eq!(ant.len(), cfg.samples);
            assert!(ant.iter().all(|&b| b == 0 || b == 1));
        }
    }

    #[test]
    fn test_hdf5_round_trip_via_observation() {
        let src = Source { name: "sun".into(), az: 0.0, el: 45.0, r: 1.5e11, jy: 1e4 };
        let positions = random_positions(3, 3.0, 5);
        let cfg = SimConfig { samples: 256, ..make_config() };
        let data = synthesize(&[src], &positions, &cfg, 20.0);

        let row_bytes = cfg.samples.div_ceil(8);
        let packed: Vec<Vec<u8>> = data.iter().map(|r| pack_row(r)).collect();
        assert_eq!(packed[0].len(), row_bytes);

        let path = std::env::temp_dir().join(format!("sim_roundtrip_{}.hdf", std::process::id()));
        let ts = "2026-08-08T00:00:00Z";
        write_observation(path.to_str().unwrap(), positions.len(), cfg.sample_rate, ts, &packed)
            .unwrap();

        let obs = crate::observation::Observation::from_hdf5(path.to_str().unwrap()).unwrap();
        assert_eq!(obs.config.num_antenna(), positions.len());
        assert_eq!(obs.get_sampling_rate(), cfg.sample_rate);
        assert_eq!(obs.data.len(), positions.len());
        for ant in &obs.data {
            assert_eq!(ant.len(), cfg.samples);
            assert!(ant.iter().all(|&b| b == 0 || b == 1));
        }
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_antenna_gains_scale_signal_amplitude() {
        // Two antennas at the same position, gains 1.0 and 2.0 → the second
        // antenna's signal amplitude at the carrier must be ~2x the first's.
        let src = Source { name: "sun".into(), az: 0.0, el: 45.0, r: 1e7, jy: 4.0 };
        let positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let cfg = SimConfig { samples: 4096, ..make_config() };
        let sigs = antenna_signals_with_gains(
            std::slice::from_ref(&src),
            &positions,
            &cfg,
            100.0, // essentially noise-free
            &[1.0, 2.0],
        );
        let f = cfg.center_freq;
        let fs = cfg.sample_rate;
        let dft_amp = |sig: &[f64]| -> f64 {
            let (mut re, mut im) = (0.0f64, 0.0f64);
            for (i, &sval) in sig.iter().enumerate() {
                let (s, c) = (TAU * f * i as f64 / fs).sin_cos();
                re += sval * c;
                im -= sval * s;
            }
            (re * re + im * im).sqrt()
        };
        let ratio = dft_amp(&sigs[1]) / dft_amp(&sigs[0]);
        assert!((ratio - 2.0).abs() < 1e-2, "gain ratio expected 2.0, got {ratio}");
    }

    #[test]
    fn test_synthesize_with_gains_unipolar() {
        let src = Source { name: "a".into(), az: 10.0, el: 20.0, r: 5e3, jy: 100.0 };
        let cfg = SimConfig { samples: 8192, ..make_config() };
        let positions = random_positions(3, 3.0, 7);
        let data = synthesize_with_gains(&[src], &positions, &cfg, 20.0, &[1.0, 0.5, 1.5]);
        assert_eq!(data.len(), 3);
        for ant in &data {
            assert_eq!(ant.len(), cfg.samples);
            assert!(ant.iter().all(|&b| b == 0 || b == 1));
        }
    }
}
