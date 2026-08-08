// Copyright (c) 2026 Tim Molteno <tim@elec.ac.nz>
// SPDX-License-Identifier: GPL-3.0

//! Simple robust statistics helpers used across acquisition modules.

/// Median of a slice of `f64` values.
///
/// Returns the middle element after sorting a copy.
/// Panics if the slice is empty.
pub fn median(data: &[f64]) -> f64 {
    assert!(!data.is_empty(), "median of empty slice");
    let mut copy: Vec<f64> = data.to_vec();
    copy.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mid = copy.len() / 2;
    if copy.len().is_multiple_of(2) {
        (copy[mid - 1] + copy[mid]) / 2.0
    } else {
        copy[mid]
    }
}

/// Median Absolute Deviation (MAD) of a slice of `f64` values about a
/// given median.
///
/// Computes `median(|x_i - med|)`.  Panics if the slice is empty.
pub fn mad(data: &[f64], med: f64) -> f64 {
    assert!(!data.is_empty(), "MAD of empty slice");
    let abs_dev: Vec<f64> = data.iter().map(|&x| (x - med).abs()).collect();
    median(&abs_dev)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_median_odd() {
        assert_eq!(median(&[1.0, 3.0, 2.0]), 2.0);
    }

    #[test]
    fn test_median_even() {
        assert_eq!(median(&[1.0, 4.0, 2.0, 3.0]), 2.5);
    }

    #[test]
    fn test_median_single() {
        assert_eq!(median(&[7.0]), 7.0);
    }

    #[test]
    fn test_mad() {
        // data: 1, 2, 3, 4, 5 → deviations from med=3: 2,1,0,1,2 → median = 1
        assert_eq!(mad(&[1.0, 2.0, 3.0, 4.0, 5.0], 3.0), 1.0);
    }

    #[test]
    fn test_mad_all_equal() {
        assert_eq!(mad(&[2.0, 2.0, 2.0], 2.0), 0.0);
    }
}
