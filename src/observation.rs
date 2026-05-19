use chrono::{DateTime, Utc};
use hdf5_reader::{Hdf5File};

use crate::config::Config;

/// A TART radio-astronomy observation.
///
/// Holds timestamp, array configuration, and per-antenna unipolar (0/1) sample data.
///
/// This mirrors the Python `tart.operation.observation.Observation` class.
pub struct Observation {
    /// UTC timestamp of the observation.
    pub timestamp: DateTime<Utc>,

    /// Array configuration (number of antennas, sampling frequency, etc.).
    pub config: Config,

    /// Per-antenna unipolar data: `data[ant_idx][sample_idx]` is 0 or 1.
    pub data: Vec<Vec<u8>>,
}

impl Observation {
    // ---------------------------------------------------------------------------
    // Constructors
    // ---------------------------------------------------------------------------

    /// Create a new `Observation` directly from its components.
    pub fn new(timestamp: DateTime<Utc>, config: Config, data: Vec<Vec<u8>>) -> Self {
        Self {
            timestamp,
            config,
            data,
        }
    }

    /// Load an observation from an HDF5 file written by `Observation.to_hdf5()`.
    ///
    /// The file layout expected:
    ///   - `config`   – scalar dataset containing a JSON string
    ///   - `timestamp` – scalar dataset containing an ISO-8601 UTC string
    ///   - `data`      – 2-D dataset where each row is packed antenna data
    ///                   (as produced by `numpy.packbits`)
    pub fn from_hdf5(filename: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let h5f = Hdf5File::open(filename)?;

        // --- config -----------------------------------------------------------
        let config_ds = h5f.dataset("config")?;
        let config_json = config_ds.read_string()?;
        let config = Config::from_json(&config_json)?;

        // --- timestamp --------------------------------------------------------
        let ts_ds = h5f.dataset("timestamp")?;
        let ts_str = ts_ds.read_string()?;
        let timestamp: DateTime<Utc> = ts_str.parse()?;

        // --- data (packed bits) -----------------------------------------------
        let data_ds = h5f.dataset("data")?;
        let shape = data_ds.shape();
        let num_antennas = shape[0] as usize;
        let row_bytes = shape[1] as usize;

        let packed: Vec<u8> = data_ds.read_raw_bytes()?;

        // Each antenna's packed data is stored contiguously; chunk into rows.
        let packed_rows: Vec<&[u8]> = packed.chunks(row_bytes).collect();

        let data: Vec<Vec<u8>> = packed_rows
            .into_iter()
            .map(|packed_row| unpack_bits(packed_row))
            .collect();

        // Validate antenna count matches config.
        if data.len() != num_antennas {
            return Err(format!(
                "data antenna count ({}) != config num_antenna ({})",
                data.len(),
                config.num_antenna()
            )
            .into());
        }

        Ok(Self {
            timestamp,
            config,
            data,
        })
    }

    // ---------------------------------------------------------------------------
    // Accessors (mirroring the Python API)
    // ---------------------------------------------------------------------------

    /// Return bipolar (±1) sample data for a single antenna.
    ///
    /// Panics if `ant_num` is out of range.
    pub fn get_antenna(&self, ant_num: usize) -> Vec<f64> {
        assert!(
            ant_num < self.config.num_antenna(),
            "antenna index {} out of range (have {})",
            ant_num,
            self.config.num_antenna()
        );
        self.data[ant_num]
            .iter()
            .map(|&b| if b == 1 { 1.0 } else { -1.0 })
            .collect()
    }

    /// Return mean value (over time) for every antenna.
    ///
    /// Maps 0→0 and 1→1, averages, then scales by 2*x-1 so the final result
    /// is in [-1, 1].
    pub fn get_means(&self) -> Vec<f64> {
        self.data
            .iter()
            .map(|ant| {
                let sum: f64 = ant.iter().map(|&b| b as f64).sum();
                (sum / ant.len() as f64) * 2.0 - 1.0
            })
            .collect()
    }

    /// Zero-lag correlation between two antennas.
    ///
    /// Converts 0/1 samples to bipolar ±1, multiplies sample-by-sample, and
    /// returns the mean.  The result is in [-1, 1] where +1 means perfect
    /// positive correlation, -1 perfect anti-correlation, and 0 no correlation.
    ///
    /// Panics if either antenna index is out of range.
    pub fn correlate(&self, i: usize, j: usize) -> f64 {
        let n = self.config.num_antenna();
        assert!(i < n, "antenna index {i} out of range (have {n})");
        assert!(j < n, "antenna index {j} out of range (have {n})");

        let a = &self.data[i];
        let b = &self.data[j];
        assert_eq!(a.len(), b.len(), "antenna sample counts differ");

        let sum: f64 = a
            .iter()
            .zip(b.iter())
            .map(|(&ai, &bj)| {
                let ai_bipolar = if ai == 1 { 1.0 } else { -1.0 };
                let bj_bipolar = if bj == 1 { 1.0 } else { -1.0 };
                ai_bipolar * bj_bipolar
            })
            .sum();

        sum / a.len() as f64
    }

    /// Sampling rate in Hz.
    pub fn get_sampling_rate(&self) -> f64 {
        self.config.sampling_frequency()
    }
}

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------

/// Unpack a byte slice where each byte encodes 8 consecutive binary samples
/// (MSB-first, matching `numpy.packbits` default behaviour).
///
/// Returns a `Vec<u8>` of 0/1 values.
fn unpack_bits(packed: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(packed.len() * 8);
    for &byte in packed {
        for bit in (0..8).rev() {
            out.push((byte >> bit) & 1);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unpack_bits() {
        // 0b10101010 -> MSB-first: 1,0,1,0,1,0,1,0
        let packed = vec![0b10101010];
        let bits = unpack_bits(&packed);
        assert_eq!(bits, vec![1, 0, 1, 0, 1, 0, 1, 0]);
    }

    #[test]
    fn test_unpack_bits_all_ones() {
        let packed = vec![0xFF, 0x00];
        let bits = unpack_bits(&packed);
        assert_eq!(&bits[..8], &[1; 8]);
        assert_eq!(&bits[8..], &[0; 8]);
    }

    #[test]
    fn test_correlate_identical() {
        // Two identical antennas → correlation = 1.0
        let data = vec![vec![1, 0, 1, 0], vec![1, 0, 1, 0]];
        let cfg = Config::from_json(r#"{"num_antenna":2,"sampling_frequency":16e6}"#).unwrap();
        let obs = Observation::new(
            chrono::Utc::now(),
            cfg,
            data,
        );
        let c = obs.correlate(0, 1);
        assert!((c - 1.0).abs() < 1e-12, "expected 1.0, got {c}");
    }

    #[test]
    fn test_correlate_anti() {
        // Two inverted antennas → correlation = -1.0
        let data = vec![vec![1, 0, 1, 0], vec![0, 1, 0, 1]];
        let cfg = Config::from_json(r#"{"num_antenna":2,"sampling_frequency":16e6}"#).unwrap();
        let obs = Observation::new(
            chrono::Utc::now(),
            cfg,
            data,
        );
        let c = obs.correlate(0, 1);
        assert!((c + 1.0).abs() < 1e-12, "expected -1.0, got {c}");
    }
}
