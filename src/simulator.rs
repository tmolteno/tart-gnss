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
    Ok(raw.into_iter().map(|v| [v[0], v[1], 0.0]).collect())
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
