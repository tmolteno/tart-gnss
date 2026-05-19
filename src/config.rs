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
