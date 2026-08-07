// Copyright (c) 2026 Tim Molteno <tim@elec.ac.nz>
// SPDX-License-Identifier: GPL-3.0

use serde::Deserialize;

/// Minimal configuration parsed from the TART JSON config stored in HDF5.
///
/// The JSON config contains many more fields; we only extract what is needed
/// for reading and interpreting observation data.
#[derive(Debug, Clone, Deserialize)]
pub struct Config {
    /// Number of antennas in the array
    #[serde(alias = "num_antenna")]
    num_antenna: usize,

    /// Sampling frequency in Hz
    #[serde(alias = "sampling_frequency")]
    sampling_frequency: f64,
}

impl Config {
    /// Parse a Config from a JSON string (as stored in the HDF5 config dataset).
    pub fn from_json(json: &str) -> Result<Self, serde_json::Error> {
        serde_json::from_str(json)
    }

    /// Number of antennas in the array.
    pub fn num_antenna(&self) -> usize {
        self.num_antenna
    }

    /// Sampling frequency in Hz.
    pub fn sampling_frequency(&self) -> f64 {
        self.sampling_frequency
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_valid() {
        let c = Config::from_json(
            r#"{"num_antenna":24,"sampling_frequency":20000000.0}"#,
        )
        .unwrap();
        assert_eq!(c.num_antenna(), 24);
        assert_eq!(c.sampling_frequency(), 20_000_000.0);
    }

    #[test]
    fn test_parse_ignores_extra_fields() {
        // Real TART configs carry many extra keys; they must be ignored.
        let c = Config::from_json(
            r#"{"num_antenna":4,"sampling_frequency":16e6,"instrument":"TART","latitude":-45.0}"#,
        )
        .unwrap();
        assert_eq!(c.num_antenna(), 4);
        assert_eq!(c.sampling_frequency(), 16e6);
    }

    #[test]
    fn test_parse_missing_required_field_fails() {
        assert!(Config::from_json(r#"{"sampling_frequency":16e6}"#).is_err());
        assert!(Config::from_json(r#"{"num_antenna":2}"#).is_err());
    }

    #[test]
    fn test_parse_malformed_json_fails() {
        assert!(Config::from_json("not json").is_err());
        assert!(Config::from_json("").is_err());
    }
}
