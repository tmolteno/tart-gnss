// Copyright (c) 2026 Tim Molteno <tim@elec.ac.nz>
// SPDX-License-Identifier: GPL-3.0

//! ACR (Acquisition Correlation Ratio) C/N0 estimation.
//!
//! Implements the method described in:
//!
//!   Ma Y, Li H, Zhou Z, Lu M (2024) "C/N0 estimation based on acquisition
//!   correlation ratio for short GNSS data."  GPS Solutions 28:143.
//!   https://doi.org/10.1007/s10291-024-01666-y
//!
//! The ACR method estimates C/N0 from the ratio of the global correlation
//! peak V^m to the largest off-peak correlation V^s at the same frequency.
//!
//! r_A = (V^m / V^s)² is the squared ratio of correlation magnitudes.
//! V^m is proportional to sqrt(I²+Q²) at the correlation peak.
//! V^s is proportional to sqrt(I²+Q²) at the largest off-peak (>1 chip).
//!
//! The lookup table below maps r_A → C/N0 (dB-Hz) using the third-order
//! hypothesis model for GPS L1 C/A (1023-chip Gold codes) with T_coh=1ms,
//! N_c=1.  Values computed from numerical integration of the Marcum Q
//! function and order-statistic CDF of non-central χ²(2, λ).

// ---------------------------------------------------------------------------
// Lookup table: r_A → C/N0 (dB-Hz)
// ---------------------------------------------------------------------------
// Pre-computed with T_coh=1ms, N_c=1, GPS L1 C/A code parameters.
// For r_A values below the first entry, C/N0 < 37 dB-Hz (cut-off region).
// For r_A values above the last entry, C/N0 saturates near ~62 dB-Hz.
static ACR_LOOKUP: &[(f64, f64)] = &[
    (1.020, 37.5),
    (1.050, 38.0),
    (1.090, 38.5),
    (1.140, 39.0),
    (1.200, 39.5),
    (1.270, 40.0),
    (1.350, 40.5),
    (1.440, 41.0),
    (1.550, 41.5),
    (1.670, 42.0),
    (1.810, 42.5),
    (1.970, 43.0),
    (2.140, 43.5),
    (2.340, 44.0),
    (2.560, 44.5),
    (2.810, 45.0),
    (3.100, 45.5),
    (3.420, 46.0),
    (3.790, 46.5),
    (4.220, 47.0),
    (4.700, 47.5),
    (5.250, 48.0),
    (5.880, 48.5),
    (6.600, 49.0),
    (7.420, 49.5),
    (8.360, 50.0),
    (9.430, 50.5),
    (10.660, 51.0),
    (12.060, 51.5),
    (13.670, 52.0),
    (15.500, 52.5),
    (17.600, 53.0),
    (20.000, 53.5),
    (22.800, 54.0),
    (26.000, 54.5),
    (29.700, 55.0),
    (34.000, 55.5),
    (38.900, 56.0),
    (44.600, 56.5),
    (51.200, 57.0),
    (58.800, 57.5),
    (67.600, 58.0),
    (77.800, 58.5),
    (89.500, 59.0),
    (103.000, 59.5),
    (118.500, 60.0),
    (136.000, 60.5),
    (156.000, 61.0),
    (179.000, 61.5),
    (200.000, 62.0),
];

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Estimate C/N0 (dB-Hz) from the acquisition correlation ratio
/// `r_a = (V^m / V^s)²` where V^m is the global correlation peak magnitude
/// and V^s is the largest off-peak (>1 chip away) correlation magnitude
/// at the same frequency.
///
/// Uses a pre-computed lookup table for T_coh = 1 ms, N_c = 1.
/// Linear interpolation is performed in log(r_A) space.
///
/// Returns `None` if r_a ≤ 1.0 (below theoretical cut-off).
pub fn estimate_cn0(r_a: f64) -> Option<f64> {
    if r_a <= 1.0 || !r_a.is_finite() {
        return None;
    }

    let table = ACR_LOOKUP;
    let first = table.first().unwrap();
    let last = table.last().unwrap();

    if r_a <= first.0 {
        return Some(first.1);
    }
    if r_a >= last.0 {
        return Some(last.1);
    }

    // Binary search for interpolation bracket
    let idx = table.partition_point(|&(r, _)| r < r_a);
    if idx == 0 {
        return Some(first.1);
    }
    if idx >= table.len() {
        return Some(last.1);
    }

    let (r_lo, cn0_lo) = table[idx - 1];
    let (r_hi, cn0_hi) = table[idx];

    // Linear interpolation in log(r_A) space for smoother results
    let log_ra = r_a.ln();
    let log_lo = r_lo.ln();
    let log_hi = r_hi.ln();
    let t = (log_ra - log_lo) / (log_hi - log_lo);
    let cn0 = cn0_lo + t * (cn0_hi - cn0_lo);

    Some(cn0.clamp(20.0, 70.0))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lookup_table_monotonic() {
        for w in ACR_LOOKUP.windows(2) {
            assert!(w[0].0 < w[1].0, "r_A not monotonic: {} -> {}", w[0].0, w[1].0);
            assert!(w[0].1 < w[1].1, "CN0 not monotonic: {} -> {}", w[0].1, w[1].1);
        }
    }

    #[test]
    fn test_estimate_cn0_below_cutoff() {
        assert!(estimate_cn0(0.5).is_none());
        assert!(estimate_cn0(1.0).is_none());
    }

    #[test]
    fn test_estimate_cn0_mid_range() {
        // r_A ≈ 3.1 corresponds to ~45.5 dB-Hz
        let cn0 = estimate_cn0(3.1);
        assert!(cn0.is_some());
        let c = cn0.unwrap();
        assert!(c > 40.0 && c < 50.0, "CN0 = {c} for r_a=3.1");
    }

    #[test]
    fn test_estimate_cn0_high() {
        let cn0 = estimate_cn0(50.0);
        assert!(cn0.is_some());
        let c = cn0.unwrap();
        assert!(c > 50.0 && c < 65.0, "CN0 = {c} for r_a=50");
    }

    #[test]
    fn test_estimate_cn0_saturation() {
        // Above max table value should clamp to max CN0
        let cn0 = estimate_cn0(1000.0);
        assert!(cn0.is_some());
        assert_eq!(cn0.unwrap(), 62.0);
    }
}
