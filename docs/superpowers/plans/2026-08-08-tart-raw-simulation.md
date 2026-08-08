# TART Raw-Data Simulation (`--simulate`) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `--simulate` mode that synthesizes raw 0/1 TART observation samples from catalogue sources (`{name, az, el, r, jy}`) with correct geometric phase delays, and writes an HDF5 observation file loadable by `tart-gnss-acquire --file`.

**Architecture:** A new `src/simulator.rs` module holds all simulation logic (types, parsing, phase delays, tone synthesis, 1-bit quantization, packing, HDF5 writing) so it is unit-testable. `main.rs` gains a `--simulate` top-level branch reusing/extending `ParsedArgs`/`parse_args`. New dependency `hdf5-writer` (pure-Rust HDF5 encoder); `hdf5-reader` bumped to 0.9.x.

**Tech Stack:** Rust (edition 2021), `hdf5-writer`/`hdf5-reader` 0.9, `serde`/`serde_json`, `chrono`, `num` (f64 math only — no extra dep needed for the math).

## Global Constraints

- The simulated HDF5 file **must** be readable by the existing `Observation::from_hdf5` (layout: `config` JSON string, `timestamp` ISO-8601 string, `data` 2-D `u8` of packed 1-bit samples, shape `[N, row_bytes]`, row-major).
- `c = 2.99793e8` m/s. Sources with `r > 1e4` m are treated as plane waves (range `r` used in delay = `1e4`).
- Azimuth 0° = North, increasing toward East; `ŝ = [sin(az)·cos(el), cos(az)·cos(el), sin(el)]` (E,N,U).
- Per-antenna geometric delay `Δ = (|r·ŝ − pos| − r)/c`; tone phase `= 2π f (t + Δ)`.
- Amplitude `A = K·√jy` (default `K = 1.0`). SNR (dB) defines Gaussian noise std `σ = √( (Σ A²/2) / 10^(snr/10) )` per antenna.
- One-bit quantization: `1 if v >= 0 else 0` (NRZ, zero→+1). Packing is MSB-first within a byte (matches `observation::unpack_bits`).
- Random antenna generation uses a **fixed default seed** (overridable via `--seed`); positions uniformly in the EN plane within radius `diameter/2`, z = 0.
- Defaults: sample rate 16.368 MHz, center freq 4.092 MHz, band 2.0 MHz, samples 65 536, diameter 3.0 m. `--snr` has **no** default (must be given).
- Deterministic behavior required (fixed seed → identical output); tests must not depend on wall-clock.

---

### Task 1: Dependencies, module scaffold, and source/position parsing

**Files:**
- Modify: `Cargo.toml`
- Create: `src/simulator.rs`
- Modify: `src/main.rs` (register `mod simulator;`)

**Interfaces:**
- Produces:
  - `pub struct Source { pub name: String, pub az: f64, pub el: f64, pub r: f64, pub jy: f64 }` (`#[derive(Debug, Clone, Deserialize)]`)
  - `pub struct SimConfig { pub sample_rate: f64, pub center_freq: f64, pub band: f64, pub samples: usize, pub gain: f64, pub seed: u64 }`
  - `pub fn parse_sources(json: &str) -> Result<Vec<Source>, serde_json::Error>`
  - `pub fn parse_positions(json: &str) -> Result<Vec<[f64; 3]>, serde_json::Error>` (z forced to 0)

- [ ] **Step 1: Update dependencies in `Cargo.toml`**

Change:
```toml
hdf5-reader = "0.5"
```
to:
```toml
hdf5-reader = "0.9"
hdf5-writer = "0.9"
```

- [ ] **Step 2: Build to confirm the reader bump is API-compatible**

Run: `cargo build --release`
Expected: compiles with no changes needed in `observation.rs` (the 0.9 reader keeps `Hdf5File::dataset`, `read_string`, `shape`→`&[u64]`, `read_raw_bytes`).
If any method name changed, update `Observation::from_hdf5` to the 0.9 names.

- [ ] **Step 3: Create `src/simulator.rs`** (types, parsing, and the PRNG)

```rust
// Copyright (c) 2026 Tim Molteno <tim@elec.ac.nz>
// SPDX-License-Identifier: GPL-3.0

//! Minimal TART raw-data simulation.
//!
//! Synthesises per-antenna 0/1 samples from catalogue sources
//! (`{name, az, el, r, jy}`) with correct geometric phase delays, and
//! writes an HDF5 observation file readable by `Observation::from_hdf5`.

use serde::Deserialize;

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
pub fn parse_positions(json: &str) -> Result<Vec<[f64; 3]>, serde_json::Error> {
    let raw: Vec<Vec<f64>> = serde_json::from_str(json)?;
    Ok(raw
        .into_iter()
        .map(|v| [v[0], v[1], 0.0])
        .collect())
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
}
```

- [ ] **Step 4: Register the module in `src/main.rs`**

Add `mod simulator;` next to the other `mod` declarations (before the `#[cfg(test)] mod testutil;` line).

- [ ] **Step 5: Run tests**

Run: `cargo test --release simulator`
Expected: the three new tests pass; existing suite still compiles/passes.

- [ ] **Step 6: Commit**

```bash
git add Cargo.toml src/simulator.rs src/main.rs Cargo.lock
git commit -m "simulator: add deps, sources/positions parsing, and test scaffolding"
```

---

### Task 2: Phase delays, frequency assignment, amplitude, and noise

**Files:**
- Modify: `src/simulator.rs`

**Interfaces:**
- Produces (used by Task 3+):
  - `pub fn source_direction(az_deg: f64, el_deg: f64) -> [f64; 3]`
  - `pub fn geometric_delay(az_deg: f64, el_deg: f64, r: f64, pos: [f64; 3]) -> f64`
  - `pub fn assign_frequencies(center: f64, band: f64, n: usize) -> Vec<f64>`
  - `pub fn amplitude(jy: f64, gain: f64) -> f64`
  - `pub fn noise_std(amplitudes: &[f64], snr_db: f64) -> f64`

- [ ] **Step 1: Add the phase/frequency/amplitude/noise functions**

Append to `src/simulator.rs` (before the `#[cfg(test)]` module):

```rust
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
/// `Δ = (|r·ŝ − pos| − r)/c`. Negative means the wavefront reaches the antenna
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
```

- [ ] **Step 2: Add unit tests**

Append to the `#[cfg(test)] mod tests` in `src/simulator.rs`:

```rust
    #[test]
    fn test_source_direction_north_horizontal() {
        // az=0 (North), el=0 (horizon) -> pure North, zero up.
        let d = source_direction(0.0, 0.0);
        assert!((d[0]).abs() < 1e-12, "east={}", d[0]);
        assert!((d[1] - 1.0).abs() < 1e-12, "north={}", d[1]);
        assert!((d[2]).abs() < 1e-12, "up={}", d[2]);
    }

    #[test]
    fn test_source_direction_zenith() {
        let d = source_direction(30.0, 90.0);
        assert!((d[2] - 1.0).abs() < 1e-12, "up={}", d[2]);
        assert!(d[0].abs() < 1e-12 && d[1].abs() < 1e-12);
    }

    #[test]
    fn test_geometric_delay_plane_wave() {
        // Zenith source far away: zero horizontal position -> zero delay.
        let d = geometric_delay(0.0, 90.0, 1e7, [0.0, 0.0, 0.0]);
        assert!(d.abs() < 1e-12);

        // Antenna 0.5 m north of origin, source at az=0 (North), el=45.
        // Plane-wave path difference ≈ -0.5*cos(el); delay = diff/c.
        let d0 = geometric_delay(0.0, 45.0, 1e7, [0.0, 0.0, 0.0]);
        let d1 = geometric_delay(0.0, 45.0, 1e7, [0.0, 0.5, 0.0]);
        let diff = (d1 - d0) * V_LIGHT;
        let expect = -0.5 * 45.0f64.to_radians().cos();
        assert!((diff - expect).abs() < 1e-6, "diff={diff} expect={expect}");
    }

    #[test]
    fn test_geometric_delay_plane_cap() {
        // A very distant source (r=1e10) behaves like r=1e4 (plane wave):
        // the delay must NOT scale with the huge range.
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
        assert_eq!(fs, vec![3.0e6, 4.0e6, 5.0e6]); // edges at center ± band/2
    }

    #[test]
    fn test_amplitude_sqrt_flux() {
        assert!((amplitude(4.0, 1.0) - 2.0).abs() < 1e-12);
        assert!((amplitude(4.0, 2.0) - 4.0).abs() < 1e-12);
        // Doubling flux scales amplitude by sqrt(2).
        let r = amplitude(8.0, 1.0) / amplitude(2.0, 1.0);
        assert!((r - 2.0f64.sqrt()).abs() < 1e-9);
    }

    #[test]
    fn test_noise_std_snr() {
        // One source amplitude 2 -> signal power 2.
        // snr=0 dB -> sigma = sqrt(2/1) = sqrt(2).
        // snr=10 dB -> sigma = sqrt(2/10) = sqrt(0.2).
        assert!((noise_std(&[2.0], 0.0) - 2.0f64.sqrt()).abs() < 1e-9);
        assert!((noise_std(&[2.0], 10.0) - 0.2f64.sqrt()).abs() < 1e-9);
    }
```

- [ ] **Step 3: Run tests**

Run: `cargo test --release simulator`
Expected: all simulator tests pass.

- [ ] **Step 4: Commit**

```bash
git add src/simulator.rs
git commit -m "simulator: add phase delays, frequency assignment, amplitude, and noise model"
```

---

### Task 3: Signal synthesis, 1-bit quantization, and bit packing

**Files:**
- Modify: `src/simulator.rs`

**Interfaces:**
- Consumes: `Source`, `SimConfig`, `XorShift64`, `gaussian`, `source_direction`, `geometric_delay`, `assign_frequencies`, `amplitude`, `noise_std`.
- Produces (used by Task 4):
  - `pub fn antenna_signals(sources: &[Source], positions: &[[f64; 3]], cfg: &SimConfig, snr_db: f64) -> Vec<Vec<f64>>`
  - `pub fn quantize(signal: &[f64]) -> Vec<u8>`
  - `pub fn synthesize(sources: &[Source], positions: &[[f64; 3]], cfg: &SimConfig, snr_db: f64) -> Vec<Vec<u8>>`
  - `pub fn pack_row(bits: &[u8]) -> Vec<u8>`

- [ ] **Step 1: Add the synthesis/quantization/packing functions**

Append to `src/simulator.rs` (before the `#[cfg(test)]` module):

```rust
/// Per-antenna continuous-valued (pre-quantization) signals.
///
/// For antenna `j`, `signals[j][n] = Σ_k A_k·cos(2π f_k (t_n + Δ_{k,j})) + noise`,
/// with `t_n = n / sample_rate`.
pub fn antenna_signals(
    sources: &[Source],
    positions: &[[f64; 3]],
    cfg: &SimConfig,
    snr_db: f64,
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
            let mut sig = Vec::with_capacity(cfg.samples);
            for i in 0..cfg.samples {
                let t = i as f64 / cfg.sample_rate;
                let mut v = 0.0f64;
                for k in 0..n_src {
                    v += amps[k] * (two_pi * freqs[k] * (t + delays[k][j])).cos();
                }
                v += gaussian(&mut rng) * sigma;
                sig.push(v);
            }
            sig
        })
        .collect()
}

/// One-bit quantize a continuous signal to 0/1 unipolar (NRZ; zero -> 1).
pub fn quantize(signal: &[f64]) -> Vec<u8> {
    signal.iter().map(|&v| if v >= 0.0 { 1 } else { 0 }).collect()
}

/// Produce raw 0/1 unipolar samples per antenna.
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
```

- [ ] **Step 2: Add tests** (including the verified phase-recovery test)

Append to the `#[cfg(test)] mod tests` in `src/simulator.rs`:

```rust
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
        // bits 1,0,1,0,1,0,1,0 -> 0b10101010
        assert_eq!(pack_row(&[1, 0, 1, 0, 1, 0, 1, 0]), vec![0b10101010]);
        // 12 bits, last 4 zero-padded: 11110000 0101_0000
        let packed = pack_row(&[1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1]);
        assert_eq!(packed, vec![0b11110000, 0b01010000]);
    }

    #[test]
    fn test_phase_delay_recovery() {
        // Single source, two antennas; recover the geometric phase from the
        // pre-quantization tones via a single-bin DFT at the source frequency.
        let src = Source { name: "t".into(), az: 30.0, el: 60.0, r: 1e7, jy: 1.0 };
        let pos0 = [0.0, 0.0, 0.0];
        let pos1 = [0.6, 0.2, 0.0];
        let cfg = SimConfig { samples: 8192, ..make_config() };
        let fs = cfg.sample_rate;
        let n = cfg.samples;
        let f = cfg.center_freq;

        let sigs = antenna_signals(&[src.clone()], &[pos0, pos1], &cfg, 100.0);
        let d0 = geometric_delay(src.az, src.el, src.r, pos0);
        let d1 = geometric_delay(src.az, src.el, src.r, pos1);

        let dft_phase = |sig: &[f64]| -> f64 {
            let (mut re, mut im) = (0.0f64, 0.0f64);
            for i in 0..n {
                let (s, c) = (TAU * f * i as f64 / fs).sin_cos();
                re += sig[i] * c;
                im -= sig[i] * s;
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
```

Note: `snr` is **not** a field of `SimConfig`; it is a function argument to `antenna_signals`/`synthesize`. There is no `snr:` field in the `SimConfig { ... }` initializers.

- [ ] **Step 3: Run tests**

Run: `cargo test --release simulator`
Expected: all simulator tests pass (including `test_phase_delay_recovery`, verified to ~1e-7 in development).

- [ ] **Step 4: Commit**

```bash
git add src/simulator.rs
git commit -m "simulator: add signal synthesis, 1-bit quantization, and packing"
```

---

### Task 4: HDF5 observation writing + round-trip test

**Files:**
- Modify: `src/simulator.rs`

**Interfaces:**
- Consumes: `synthesize`, `pack_row`, `Source`, `SimConfig`, `random_positions`, `parse_positions`, `parse_sources`.
- Produces:
  - `pub fn write_observation(path: &str, num_antenna: usize, sampling_frequency: f64, timestamp: &str, data: &[Vec<u8>]) -> Result<(), Box<dyn std::error::Error>>`

- [ ] **Step 1: Add the HDF5 writer**

Append to `src/simulator.rs` (before the `#[cfg(test)]` module), and add at the top of the file:

```rust
use hdf5_writer::{DatasetBuilder, Hdf5Builder, Hdf5Writer, WriteOptions};
```

```rust
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
    let data_ds =
        DatasetBuilder::typed_data::<u8>("data", [num_antenna as u64, row_bytes as u64], &flat)?;

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
```

- [ ] **Step 2: Add a round-trip test**

Append to `#[cfg(test)] mod tests` in `src/simulator.rs`:

```rust
    #[test]
    fn test_hdf5_round_trip_via_observation() {
        // Build a small synthetic observation, write it, and read it back with
        // the crate's own Observation reader to prove the file is consumable.
        let src = Source { name: "sun".into(), az: 0.0, el: 45.0, r: 1.5e11, jy: 1e4 };
        let positions = random_positions(3, 3.0, 5);
        let cfg = SimConfig { samples: 256, ..make_config() };
        let data = synthesize(&[src], &positions, &cfg, 20.0);

        // Pack each antenna row and write.
        let row_bytes = (cfg.samples + 7) / 8;
        let packed: Vec<Vec<u8>> = data.iter().map(|r| pack_row(r)).collect();
        assert_eq!(packed[0].len(), row_bytes);

        let path = std::env::temp_dir().join(format!("sim_roundtrip_{}.hdf", std::process::id()));
        let ts = "2026-08-08T00:00:00Z";
        write_observation(path.to_str().unwrap(), positions.len(), cfg.sample_rate, ts, &packed)
            .unwrap();

        // Read back with the production reader.
        use crate::config::Config as _; // silence unused if needed
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
```

Note: `make_config` from Task 3 already exists in the test module; reuse it. `SimConfig` here uses `samples: 256`, so use struct-update syntax `SimConfig { samples: 256, ..make_config() }` (there is no `snr` field).

- [ ] **Step 3: Run tests**

Run: `cargo test --release simulator`
Expected: the round-trip test passes (verifying the written file loads via `Observation::from_hdf5`).

If the reader's `shape()` returns dims differently than the write shape, adjust `write_observation` so `[num_antenna, row_bytes]` matches what `from_hdf5` reads (`shape[0]=antennas`, `shape[1]=row_bytes`).

- [ ] **Step 4: Commit**

```bash
git add src/simulator.rs
git commit -m "simulator: write HDF5 observation and add round-trip test"
```

---

### Task 5: CLI integration (`--simulate` in `main.rs`)

**Files:**
- Modify: `src/main.rs`

**Interfaces:**
- Consumes: `simulator::{parse_sources, parse_positions, random_positions, SimConfig, synthesize, pack_row, write_observation, DEFAULT_SEED}`, `chrono::Utc`.
- Produces: extended `ParsedArgs` fields + `parse_args` flags + a `--simulate` branch in `main`.

- [ ] **Step 1: Add command-line flags to `ParsedArgs` and `parse_args`**

In `ParsedArgs` (in `src/main.rs`) add these fields:

```rust
    simulate: bool,
    sources: Option<String>,
    positions: Option<String>,
    n_generate: Option<usize>,
    diameter: Option<f64>,
    sim_seed: Option<u64>,
    sim_samples: Option<usize>,
    sample_rate: Option<f64>,
    center_freq: Option<f64>,
    band: Option<f64>,
    snr: Option<f64>,
```

In `parse_args`, add match arms:

```rust
            "--simulate" => p.simulate = true,
            "--sources" => {
                let v = flag_value(args, &mut i, "--sources")?;
                p.sources = Some(v.to_string());
            }
            "--positions" => {
                let v = flag_value(args, &mut i, "--positions")?;
                p.positions = Some(v.to_string());
            }
            "--N" => {
                let v = flag_value(args, &mut i, "--N")?;
                p.n_generate = Some(v.parse().map_err(|_| format!("invalid integer for --N: {v}"))?);
            }
            "--diameter" => {
                let v = flag_value(args, &mut i, "--diameter")?;
                p.diameter = Some(v.parse().map_err(|_| format!("invalid float for --diameter: {v}"))?);
            }
            "--seed" => {
                let v = flag_value(args, &mut i, "--seed")?;
                p.sim_seed = Some(v.parse().map_err(|_| format!("invalid integer for --seed: {v}"))?);
            }
            "--sample-rate" => {
                let v = flag_value(args, &mut i, "--sample-rate")?;
                p.sample_rate = Some(v.parse().map_err(|_| format!("invalid float for --sample-rate: {v}"))?);
            }
            "--center-freq" => {
                let v = flag_value(args, &mut i, "--center-freq")?;
                p.center_freq = Some(v.parse().map_err(|_| format!("invalid float for --center-freq: {v}"))?);
            }
            "--band" => {
                let v = flag_value(args, &mut i, "--band")?;
                p.band = Some(v.parse().map_err(|_| format!("invalid float for --band: {v}"))?);
            }
            "--snr" => {
                let v = flag_value(args, &mut i, "--snr")?;
                p.snr = Some(v.parse().map_err(|_| format!("invalid float for --snr: {v}"))?);
            }
```

Note: `--samples` is used by both benchmark and simulation contexts; keep the existing `--samples`-style value under `sim_samples`, or reuse an existing field if present. (If a prior task added a `samples` field, reuse it; otherwise add `sim_samples` as above.)

Also add these to the destructure in `main` (Step 3).

- [ ] **Step 2: Update the destructuring and the file-requirement check**

Destructure the new fields out of `ParsedArgs` (adding them to the existing `let ParsedArgs { ... } = parsed;`):

```rust
        simulate,
        sources,
        positions,
        n_generate,
        diameter,
        sim_seed,
        sim_samples,
        sample_rate,
        center_freq,
        band,
        snr,
```

Change the file-requirement so `--simulate` (like `--benchmark`) does not require `--file`:

```rust
    let file = if benchmark_flag || simulate {
        file
    } else {
        Some(file.unwrap_or_else(|| { ... }))
    };
```

- [ ] **Step 3: Add the `--simulate` branch**

Insert this right after the `if all_flag { ... }` block and before the `let file = ...` line:

```rust
    // --- Simulation mode -------------------------------------------------
    if simulate {
        use crate::simulator::{
            parse_positions, parse_sources, pack_row, random_positions, synthesize,
            write_observation, SimConfig, DEFAULT_SEED,
        };

        let src_path = sources.unwrap_or_else(|| {
            eprintln!("--simulate requires --sources <catalogue.json>");
            print_usage(&args[0]);
            std::process::exit(1);
        });
        let src_json = std::fs::read_to_string(&src_path)
            .unwrap_or_else(|e| {
                eprintln!("error reading sources file {src_path}: {e}");
                std::process::exit(1);
            });
        let src_list = parse_sources(&src_json).unwrap_or_else(|e| {
            eprintln!("error parsing sources file {src_path}: {e}");
            std::process::exit(1);
        });
        if src_list.is_empty() {
            eprintln!("--sources file contains no sources");
            std::process::exit(1);
        }

        let positions: Vec<[f64; 3]> = if let Some(pos_path) = positions {
            let pos_json = std::fs::read_to_string(&pos_path)
                .unwrap_or_else(|e| {
                    eprintln!("error reading positions file {pos_path}: {e}");
                    std::process::exit(1);
                });
            parse_positions(&pos_json).unwrap_or_else(|e| {
                eprintln!("error parsing positions file {pos_path}: {e}");
                std::process::exit(1);
            })
        } else {
            let n = n_generate.unwrap_or_else(|| {
                eprintln!("--simulate needs either --positions <file> or --N <int>");
                print_usage(&args[0]);
                std::process::exit(1);
            });
            random_positions(n, diameter.unwrap_or(3.0), sim_seed.unwrap_or(DEFAULT_SEED))
        };
        if positions.is_empty() {
            eprintln!("no antenna positions");
            std::process::exit(1);
        }

        let cfg = SimConfig {
            sample_rate: sample_rate.unwrap_or(16.368e6),
            center_freq: center_freq.unwrap_or(4.092e6),
            band: band.unwrap_or(2.0e6),
            samples: sim_samples.unwrap_or(65_536),
            gain: 1.0,
            seed: sim_seed.unwrap_or(DEFAULT_SEED),
        };
        let snr_db = snr.unwrap_or_else(|| {
            eprintln!("--simulate requires --snr <dB>");
            print_usage(&args[0]);
            std::process::exit(1);
        });

        eprintln!(
            "Simulating {} sources on {} antennas (fs={} Hz, fc={} Hz, band={} Hz, {} samples, snr={snr_db} dB)",
            src_list.len(),
            positions.len(),
            cfg.sample_rate,
            cfg.center_freq,
            cfg.band,
            cfg.samples,
        );
        for (k, s) in src_list.iter().enumerate() {
            let f = cfg.center_freq + cfg.band * (k as f64 / (src_list.len() as f64 - 1.0) - 0.5);
            eprintln!(
                "  {:>20}: az={:6.2} el={:6.2} r={:.3e} m jy={:.3e} f={:.3e} Hz",
                s.name, s.az, s.el, s.r, s.jy, f
            );
        }

        let raw = synthesize(&src_list, &positions, &cfg, snr_db);
        let packed: Vec<Vec<u8>> = raw.iter().map(|r| pack_row(r)).collect();

        let out_path = output_file.unwrap_or("simulation.hdf".to_string());
        let ts = chrono::Utc::now().to_rfc3339();
        write_observation(&out_path, positions.len(), cfg.sample_rate, &ts, &packed)
            .unwrap_or_else(|e| {
                eprintln!("error writing simulation file {out_path}: {e}");
                std::process::exit(1);
            });
        eprintln!("Wrote simulated observation to {out_path}");
        return;
    }
```

- [ ] **Step 4: Add `parse_args` unit tests for the new flags**

In the `#[cfg(test)] mod tests` in `src/main.rs`:

```rust
    #[test]
    fn test_parse_simulate_flags() {
        let p = parse_args(&a(&[
            "prog", "--simulate", "--sources", "cat.json",
            "--positions", "ants.json", "--diameter", "3.0",
            "--seed", "99", "--N", "12",
            "--sample-rate", "16.368e6", "--center-freq", "4.092e6",
            "--band", "2e6", "--snr", "20",
            "--output", "obs.hdf",
        ]))
        .unwrap();
        assert!(p.simulate);
        assert_eq!(p.sources.as_deref(), Some("cat.json"));
        assert_eq!(p.positions.as_deref(), Some("ants.json"));
        assert_eq!(p.n_generate, Some(12));
        assert_eq!(p.diameter, Some(3.0));
        assert_eq!(p.sim_seed, Some(99));
        assert_eq!(p.center_freq, Some(4.092e6));
        assert_eq!(p.snr, Some(20.0));
        assert_eq!(p.output_file.as_deref(), Some("obs.hdf"));
    }
```

- [ ] **Step 5: Build and run the full suite**

Run: `cargo build --release` then `cargo test --release`.
Expected: compiles cleanly; all tests pass.

- [ ] **Step 6: Smoke-test the CLI end to end**

Create a `data/sources_test.json` and `data/positions_test.json` (a couple of sources and antennas), then run:

```bash
cargo run --release -- --simulate --sources data/sources_test.json \
    --positions data/positions_test.json --snr 20 --output /tmp/sim_test.hdf
```

Then verify it loads with the real acquisition reader:

```bash
cargo run --release -- --file /tmp/sim_test.hdf --i 0 --j 1
```

Expected: `--simulate` prints the summary and writes the file; the `--i/--j` run loads it and prints a correlation.

- [ ] **Step 7: Commit**

```bash
git add src/main.rs data/sources_test.json data/positions_test.json
git commit -m "simulate: add --simulate CLI mode with HDF5 output"
```

---

### Task 6: Final verification

**Files:**
- None (verification only).

- [ ] **Step 1: Run the full test suite in release mode**

Run: `cargo test --release`
Expected: all tests pass (existing 98 + new simulator tests).

- [ ] **Step 2: Run clippy (best effort)**

Run: `cargo clippy --all-targets`
Expected: no *new* warnings introduced by this feature (pre-existing lints may remain; do not change behavior to silence them).

- [ ] **Step 3: Confirm a clean worktree and review the log**

Run: `git status` (clean) and `git log --oneline -12` to review the feature commits.

---

## Self-Review

**Spec coverage:**
- Catalogue source parsing (`{name, az, el, r, jy}`) → Task 1 ✓
- Correct phase delays (near-field + plane-wave cap) → Task 2 + `test_phase_delay_recovery` ✓
- Antenna positions (parse ENU, z=0; random from `--N` within `diameter/2`, fixed seed) → Tasks 1, 5 ✓
- Minimal coherent tones at distinct frequencies + `A=√jy` + SNR + 1-bit quantization → Tasks 2, 3 ✓
- HDF5 output loadable by `Observation::from_hdf5` → Task 4 round-trip ✓
- CLI surface (defaults 16.368 MHz / 4.092 MHz / ±1 MHz / 2¹⁶ samples / 3 m / fixed seed; required `--snr`) → Task 5 ✓
- Defaults all match spec ✓

**Placeholder scan:** No TBD/TODO; every code step is concrete. Numeric formulas (delay, amplitude, noise σ, frequency spread) verified in development spikes (phase recovery ~1e-7 rad; HDF5 round-trip confirmed).

**Type consistency:** `SimConfig` has no `snr` field (SNR is a function argument) — the plan flags the invalid illustrative fragment to avoid the error; `make_config()` helper shared across tests. `geometric_delay`, `assign_frequencies`, `amplitude`, `noise_std`, `antenna_signals`, `quantize`, `synthesize`, `pack_row`, `write_observation` signatures are consistent across tasks. Order of index args to `geometric_delay` is `(az, el, r, pos)` everywhere.
