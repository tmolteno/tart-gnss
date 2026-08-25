// Copyright (c) 2026 Tim Molteno <tim@elec.ac.nz>
// SPDX-License-Identifier: GPL-3.0

//! Antenna-relative C/N0 quality testing (`--test-antennas`).
//!
//! Processes acquisition results (per-PRN ACR C/N0 estimates) to rank the
//! array's antennas by relative performance:
//!
//! 1. Find satellites "commonly visible with good signal strength" — PRNs whose
//!    ACR C/N0 estimate is ≥ `min_cn0` on **at least one** antenna.  A broken
//!    antenna never reaches the threshold itself, and need not produce an
//!    estimate at all: satellites other antennas see strongly still qualify.
//!    This is the reference set.
//! 2. For each antenna, its quality score is the median C/N0 across the
//!    reference set (robust to outlier satellites).  Antennas that fall below
//!    the threshold still contribute their (low) C/N0 values, so every antenna
//!    with at least one estimate is scored; one with none at all (a dead
//!    channel) is reported unscored (`null` median/offset/rank).
//! 3. Report each antenna's score relative to the best antenna
//!    (`offset_db = median − best`, so the best antenna has offset 0) and an
//!    overall rank.

use serde::Serialize;

use crate::CombinedOutput;

/// Per-antenna summary statistics.
#[derive(Debug, Clone, Serialize)]
pub struct AntennaStat {
    /// Antenna index.
    pub antenna: usize,
    /// Number of reference satellites for which this antenna met the C/N0
    /// threshold.  A broken or dead antenna scores 0 here even though it
    /// still appears in the report.
    pub n_sats: usize,
    /// Median C/N0 across reference satellites (dB-Hz), `null` if this
    /// antenna produced no ACR estimate at all (e.g. a dead channel).
    pub median_cn0_db_hz: Option<f64>,
    /// Median minus the best antenna's median (dB). Best antenna has 0.0;
    /// `null` for unscored antennas.
    pub offset_db: Option<f64>,
    /// Rank by median C/N0 (1 = best; ties broken by antenna index);
    /// `null` for unscored antennas.
    pub rank: Option<usize>,
}

/// Per-satellite C/N0 matrix (one entry per reference satellite).
#[derive(Debug, Clone, Serialize)]
pub struct SatelliteCn0 {
    /// Constellation label, e.g. "GPS03".
    pub sv: String,
    /// Per-antenna C/N0 estimates aligned with `antenna_numbers` (dB-Hz).
    /// `null` where an antenna produced no ACR estimate for this satellite.
    pub cn0_db_hz: Vec<Option<f64>>,
}

/// JSON report produced by `--test-antennas`.
#[derive(Debug, Clone, Serialize)]
pub struct AntennaTestOutput {
    /// Antenna indices the report covers, in order (aligned with the
    /// per-antenna fields of each satellite entry).
    pub antenna_numbers: Vec<usize>,
    /// C/N0 threshold used to select reference satellites (dB-Hz).
    pub min_cn0_db_hz: f64,
    /// Number of reference (commonly visible, strong) satellites used.
    pub n_reference_satellites: usize,
    /// Labels of the reference satellites.
    pub reference_satellites: Vec<String>,
    /// Per-antenna stat rows, ordered by antenna index.
    pub antennas: Vec<AntennaStat>,
    /// Per-satellite C/N0 matrix for the reference set.
    pub per_satellite: Vec<SatelliteCn0>,
}

/// Uniform read view over the six constellation output types.
trait AllOutput {
    fn antenna_numbers(&self) -> &[usize];
    fn rows(&self) -> Vec<(String, Option<Vec<Option<f64>>>)>;
}

macro_rules! impl_all_output {
    ($ty:ty) => {
        impl AllOutput for $ty {
            fn antenna_numbers(&self) -> &[usize] {
                &self.antenna_numbers
            }
            fn rows(&self) -> Vec<(String, Option<Vec<Option<f64>>>)> {
                self.results
                    .iter()
                    .map(|r| (r.sv.clone(), r.cn0_acr.clone()))
                    .collect()
            }
        }
    };
}

impl_all_output!(crate::acquisition::GpsAllAcquisitionOutput);
impl_all_output!(crate::galileo::GalileoAllAcquisitionOutput);
impl_all_output!(crate::beidou::BeiDouAllAcquisitionOutput);
impl_all_output!(crate::sbas::SbasAllAcquisitionOutput);
impl_all_output!(crate::l1c::L1CAllAcquisitionOutput);
impl_all_output!(crate::qzss::QzssAllAcquisitionOutput);

/// Build the antenna-relative C/N0 report from acquisition results.
pub fn run(output: &CombinedOutput, min_cn0: f64) -> AntennaTestOutput {
    let mut antenna_numbers: Vec<usize> = Vec::new();
    let mut candidates: Vec<(String, Vec<Option<f64>>)> = Vec::new();

    macro_rules! gather {
        ($field:ident) => {
            if let Some(all) = &output.$field {
                if antenna_numbers.is_empty() {
                    antenna_numbers = all.antenna_numbers().to_vec();
                }
                for (sv, cn0s) in all.rows() {
                    if let Some(cn0s) = cn0s {
                        // Rows are aligned with the output's antenna_numbers
                        // (None where an antenna produced no ACR estimate, e.g.
                        // a dead channel).  A satellite qualifies if ANY
                        // antenna reaches the threshold; antennas with no
                        // estimate do not disqualify it.
                        if cn0s.iter().any(|c| c.is_some_and(|c| c >= min_cn0)) {
                            candidates.push((sv, cn0s));
                        }
                    }
                }
            }
        };
    }
    gather!(gps);
    gather!(galileo);
    gather!(beidou);
    gather!(sbas);
    gather!(l1c);
    gather!(qzss);

    let n_ref = candidates.len();
    let reference_satellites: Vec<String> =
        candidates.iter().map(|(sv, _)| sv.clone()).collect();
    let per_satellite: Vec<SatelliteCn0> = candidates
        .iter()
        .map(|(sv, cn0s)| SatelliteCn0 {
            sv: sv.clone(),
            cn0_db_hz: cn0s.clone(),
        })
        .collect();

    if candidates.is_empty() {
        return AntennaTestOutput {
            antenna_numbers,
            min_cn0_db_hz: min_cn0,
            n_reference_satellites: 0,
            reference_satellites: Vec::new(),
            antennas: Vec::new(),
            per_satellite: Vec::new(),
        };
    }

    // Per-antenna C/N0 values across the reference set (column-major).
    // Antennas with no values at all (e.g. a dead channel) end up empty.
    let mut per_ant: Vec<(usize, Vec<f64>)> =
        antenna_numbers.iter().map(|&a| (a, Vec::new())).collect();
    for (_, cn0s) in &candidates {
        for (col, c) in cn0s.iter().enumerate() {
            if let Some(c) = c {
                per_ant[col].1.push(*c);
            }
        }
    }

    // Rank scored antennas by median C/N0 (descending), ties broken by
    // antenna index.  Unscored antennas (no ACR values at all) come last.
    let mut ranked: Vec<(usize, f64)> = per_ant
        .iter()
        .filter(|(_, vals)| !vals.is_empty())
        .map(|(a, vals)| (*a, crate::stats::median(vals)))
        .collect();
    ranked.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
    let best = ranked[0].1;

    let mut antennas: Vec<AntennaStat> = ranked
        .iter()
        .enumerate()
        .map(|(rank, (ant, med))| {
            let col = antenna_numbers.iter().position(|&a| a == *ant).unwrap();
            let n_sats = per_ant[col]
                .1
                .iter()
                .filter(|&&c| c >= min_cn0)
                .count();
            AntennaStat {
                antenna: *ant,
                n_sats,
                median_cn0_db_hz: Some(*med),
                offset_db: Some(*med - best),
                rank: Some(rank + 1),
            }
        })
        .collect();
    // Unscored antennas (all ACR values missing) still appear, with nulls.
    for (a, vals) in &per_ant {
        if vals.is_empty() {
            antennas.push(AntennaStat {
                antenna: *a,
                n_sats: 0,
                median_cn0_db_hz: None,
                offset_db: None,
                rank: None,
            });
        }
    }
    antennas.sort_by_key(|s| s.antenna); // stable output order

    AntennaTestOutput {
        antenna_numbers,
        min_cn0_db_hz: min_cn0,
        n_reference_satellites: n_ref,
        reference_satellites,
        antennas,
        per_satellite,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::acquisition::{
        acquire_all_gps, generate_ca_code, CA_CHIPS, GPS_IF, GPS_SEARCH_BAND,
    };
    use crate::config::Config;
    use crate::observation::Observation;
    use crate::simulator::quantize;
    use std::f64::consts::TAU;

    /// Build a GPS per-PRN result row directly (no acquisition needed).
    fn gps_row(sv: &str, cn0s: Option<Vec<Option<f64>>>) -> crate::acquisition::GpsPrnResult {
        crate::acquisition::GpsPrnResult {
            sv: sv.into(),
            strengths: Vec::new(),
            phases: Vec::new(),
            freqs: Vec::new(),
            cn0_acr: cn0s,
            phase_median: None,
            phase_mad: None,
            freq_median: None,
            freq_mad: None,
        }
    }

    /// Build a *Galileo* per-SV row so multi-constellation collection is tested.
    fn gal_row(sv: &str, cn0s: Option<Vec<Option<f64>>>) -> crate::galileo::GalileoPrnResult {
        crate::galileo::GalileoPrnResult {
            sv: sv.into(),
            strengths: Vec::new(),
            phases: Vec::new(),
            freqs: Vec::new(),
            cn0_acr: cn0s,
            phase_median: None,
            phase_mad: None,
            freq_median: None,
            freq_mad: None,
        }
    }

    fn output(
        gps: crate::acquisition::GpsAllAcquisitionOutput,
        galileo: crate::galileo::GalileoAllAcquisitionOutput,
    ) -> CombinedOutput {
        CombinedOutput {
            gps: Some(gps),
            galileo: Some(galileo),
            beidou: None,
            sbas: None,
            l1c: None,
            qzss: None,
        }
    }

    #[test]
    fn test_run_basic_ranking_and_offsets() {
        // Two reference satellites, three antennas.
        let gps = crate::acquisition::GpsAllAcquisitionOutput {
            antenna_numbers: vec![0, 1, 2],
            results: vec![
                gps_row("GPS01", Some(vec![Some(40.0), Some(41.0), Some(42.0)])),
                gps_row("GPS02", Some(vec![Some(43.0), Some(44.0), Some(45.0)])),
            ],
        };
        let gal = crate::galileo::GalileoAllAcquisitionOutput {
            antenna_numbers: vec![0, 1, 2],
            results: vec![],
        };
        let report = run(&output(gps, gal), 40.0);
        assert_eq!(report.antenna_numbers, vec![0, 1, 2]);
        assert_eq!(report.min_cn0_db_hz, 40.0);
        assert_eq!(report.n_reference_satellites, 2);
        assert_eq!(report.reference_satellites, vec!["GPS01", "GPS02"]);
        assert_eq!(report.per_satellite.len(), 2);
        assert_eq!(
            report.per_satellite[0].cn0_db_hz,
            vec![Some(40.0), Some(41.0), Some(42.0)]
        );
        // median per antenna: ant0 = 41.5, ant1 = 42.5, ant2 = 43.5 → best is ant2
        let stats = &report.antennas;
        assert_eq!(stats.len(), 3);
        assert_eq!(stats[0].antenna, 0);
        assert_eq!(stats[0].n_sats, 2);
        assert_eq!(stats[0].median_cn0_db_hz, Some(41.5));
        assert_eq!(stats[0].offset_db, Some(-2.0));
        assert_eq!(stats[0].rank, Some(3));
        assert_eq!(stats[2].antenna, 2);
        assert_eq!(stats[2].offset_db, Some(0.0), "best antenna offset 0");
        assert_eq!(stats[2].rank, Some(1));
    }

    #[test]
    fn test_run_threshold_filters_weak_satellites() {
        let gps = crate::acquisition::GpsAllAcquisitionOutput {
            antenna_numbers: vec![0, 1],
            results: vec![
                gps_row("GPS01", Some(vec![Some(40.0), Some(41.0)])),
                gps_row("GPS02", Some(vec![Some(45.0), Some(46.0)])),
            ],
        };
        let gal = crate::galileo::GalileoAllAcquisitionOutput {
            antenna_numbers: vec![0, 1],
            results: vec![],
        };
        let report = run(&output(gps, gal), 43.0);
        assert_eq!(report.n_reference_satellites, 1);
        assert_eq!(report.reference_satellites, vec!["GPS02"]);
        assert_eq!(report.per_satellite[0].sv, "GPS02");
        assert_eq!(report.antennas[0].n_sats, 1);
    }

    #[test]
    fn test_run_missing_antenna_nulls() {
        // Antenna 2 produced no ACR estimate for this satellite (and none at
        // all in this output): the satellite qualifies as soon as ≥ 1 antenna
        // reaches the threshold, its row keeps the None in place (aligned
        // with antenna_numbers), and the antenna with no values at all is
        // reported unscored (nulls).
        let gps = crate::acquisition::GpsAllAcquisitionOutput {
            antenna_numbers: vec![0, 1, 2],
            results: vec![gps_row("GPS01", Some(vec![Some(45.0), Some(46.0), None]))],
        };
        let gal = crate::galileo::GalileoAllAcquisitionOutput {
            antenna_numbers: vec![0, 1, 2],
            results: vec![],
        };
        let report = run(&output(gps, gal), 40.0);
        assert_eq!(report.n_reference_satellites, 1);
        assert_eq!(report.reference_satellites, vec!["GPS01"]);
        assert_eq!(
            report.per_satellite[0].cn0_db_hz,
            vec![Some(45.0), Some(46.0), None]
        );
        // antenna 2 has no values at all → unscored (nulls); antenna 1 best.
        let antenna_by_idx = |i: usize| report.antennas.iter().find(|a| a.antenna == i).unwrap();
        assert_eq!(antenna_by_idx(1).rank, Some(1));
        assert_eq!(antenna_by_idx(0).rank, Some(2));
        assert_eq!(antenna_by_idx(0).n_sats, 1);
        assert_eq!(antenna_by_idx(1).n_sats, 1);
        assert_eq!(antenna_by_idx(2).n_sats, 0);
        assert_eq!(antenna_by_idx(2).median_cn0_db_hz, None);
        assert_eq!(antenna_by_idx(2).offset_db, None);
        assert_eq!(antenna_by_idx(2).rank, None);
    }

    #[test]
    fn test_run_middle_antenna_null_alignment() {
        // A None in the middle of a row must not shift the values that
        // follow it: each entry stays bound to its antenna column.  (This is
        // the regression guard for the ACR rows being aligned with
        // antenna_numbers at the source.)
        let gps = crate::acquisition::GpsAllAcquisitionOutput {
            antenna_numbers: vec![0, 1, 2, 3],
            results: vec![gps_row(
                "GPS01",
                Some(vec![Some(45.0), None, Some(46.0), Some(47.0)]),
            )],
        };
        let gal = crate::galileo::GalileoAllAcquisitionOutput {
            antenna_numbers: vec![0, 1, 2, 3],
            results: vec![],
        };
        let report = run(&output(gps, gal), 40.0);
        assert_eq!(report.n_reference_satellites, 1);
        assert_eq!(
            report.per_satellite[0].cn0_db_hz,
            vec![Some(45.0), None, Some(46.0), Some(47.0)]
        );
        // ant1 is unscored (no values); ant3 best (47.0), ant2 next, ant0 last.
        let antenna_by_idx = |i: usize| report.antennas.iter().find(|a| a.antenna == i).unwrap();
        assert_eq!(antenna_by_idx(3).rank, Some(1));
        assert_eq!(antenna_by_idx(3).median_cn0_db_hz, Some(47.0));
        assert_eq!(antenna_by_idx(2).rank, Some(2));
        assert_eq!(antenna_by_idx(2).median_cn0_db_hz, Some(46.0));
        assert_eq!(antenna_by_idx(0).rank, Some(3));
        assert_eq!(antenna_by_idx(0).median_cn0_db_hz, Some(45.0));
        assert_eq!(antenna_by_idx(1).rank, None);
        assert_eq!(antenna_by_idx(1).median_cn0_db_hz, None);
    }

    #[test]
    fn test_report_json_serialization_shape() {
        // JSON contract of the --test-antennas report: per-satellite rows
        // keep one entry per antenna (null where ACR failed, so the matrix
        // stays aligned with antenna_numbers) and an antenna with no
        // estimates at all is unscored with null fields.
        let gps = crate::acquisition::GpsAllAcquisitionOutput {
            antenna_numbers: vec![0, 1, 2],
            results: vec![gps_row("GPS01", Some(vec![Some(45.0), None, Some(46.0)]))],
        };
        let gal = crate::galileo::GalileoAllAcquisitionOutput {
            antenna_numbers: vec![0, 1, 2],
            results: vec![],
        };
        let report = run(&output(gps, gal), 40.0);
        let json = serde_json::to_value(&report).expect("serialize antenna-test report");
        assert_eq!(json["n_reference_satellites"], serde_json::json!(1));
        assert_eq!(json["antenna_numbers"], serde_json::json!([0, 1, 2]));
        assert_eq!(json["min_cn0_db_hz"], serde_json::json!(40.0));
        assert_eq!(json["reference_satellites"], serde_json::json!(["GPS01"]));
        assert_eq!(json["per_satellite"][0]["sv"], serde_json::json!("GPS01"));
        assert_eq!(
            json["per_satellite"][0]["cn0_db_hz"],
            serde_json::json!([45.0, null, 46.0])
        );
        // antenna 1 (middle): no ACR estimates at all → unscored, all nulls.
        assert_eq!(json["antennas"][1]["antenna"], serde_json::json!(1));
        assert_eq!(json["antennas"][1]["n_sats"], serde_json::json!(0));
        assert!(json["antennas"][1]["median_cn0_db_hz"].is_null());
        assert!(json["antennas"][1]["offset_db"].is_null());
        assert!(json["antennas"][1]["rank"].is_null());
        // antennas 0 and 2 are scored with numbers (ant2 best, 46.0).
        assert!(json["antennas"][0]["median_cn0_db_hz"].is_number());
        assert_eq!(json["antennas"][0]["rank"], serde_json::json!(2));
        assert_eq!(json["antennas"][2]["rank"], serde_json::json!(1));
        // The raw document keeps the nulls (no field skipping).
        let doc = serde_json::to_string(&report).unwrap();
        assert!(doc.contains("\"median_cn0_db_hz\":null"), "{doc}");
    }

    #[test]
    fn test_run_missing_cn0_gives_empty_report() {
        let gps = crate::acquisition::GpsAllAcquisitionOutput {
            antenna_numbers: vec![0, 1],
            results: vec![gps_row("GPS01", None)],
        };
        let gal = crate::galileo::GalileoAllAcquisitionOutput {
            antenna_numbers: vec![0, 1],
            results: vec![],
        };
        let report = run(&output(gps, gal), 40.0);
        assert_eq!(report.n_reference_satellites, 0);
        assert!(report.antennas.is_empty());
    }

    #[test]
    fn test_run_no_outputs_at_all() {
        let report = run(&CombinedOutput {
            gps: None,
            galileo: None,
            beidou: None,
            sbas: None,
            l1c: None,
            qzss: None,
        }, 40.0);
        assert!(report.antenna_numbers.is_empty());
        assert_eq!(report.n_reference_satellites, 0);
    }

    #[test]
    fn test_run_single_antenna_carries_satellite() {
        // A broken antenna (30 dB-Hz) must not disqualify a satellite that
        // another antenna sees strongly (44 dB-Hz ≥ threshold).
        let gps = crate::acquisition::GpsAllAcquisitionOutput {
            antenna_numbers: vec![0, 1],
            results: vec![gps_row("GPS01", Some(vec![Some(44.0), Some(30.0)]))],
        };
        let gal = crate::galileo::GalileoAllAcquisitionOutput {
            antenna_numbers: vec![0, 1],
            results: vec![],
        };
        let report = run(&output(gps, gal), 44.0);
        assert_eq!(report.n_reference_satellites, 1);
        assert_eq!(report.reference_satellites, vec!["GPS01"]);
        // Both antennas are scored; the broken one shows its low median and
        // 0 satellites meeting the threshold.
        assert_eq!(report.antennas[0].antenna, 0);
        assert_eq!(report.antennas[0].n_sats, 1);
        assert_eq!(report.antennas[0].median_cn0_db_hz, Some(44.0));
        assert_eq!(report.antennas[0].offset_db, Some(0.0));
        assert_eq!(report.antennas[0].rank, Some(1));
        assert_eq!(report.antennas[1].antenna, 1);
        assert_eq!(report.antennas[1].n_sats, 0);
        assert_eq!(report.antennas[1].median_cn0_db_hz, Some(30.0));
        assert_eq!(report.antennas[1].offset_db, Some(-14.0));
        assert_eq!(report.antennas[1].rank, Some(2));
    }

    #[test]
    fn test_run_n_sats_varies_per_antenna() {
        // Both satellites carried by antenna 0 alone.
        let gps = crate::acquisition::GpsAllAcquisitionOutput {
            antenna_numbers: vec![0, 1],
            results: vec![
                gps_row("GPS01", Some(vec![Some(44.0), Some(30.0)])),
                gps_row("GPS02", Some(vec![Some(45.0), Some(31.0)])),
            ],
        };
        let gal = crate::galileo::GalileoAllAcquisitionOutput {
            antenna_numbers: vec![0, 1],
            results: vec![],
        };
        let report = run(&output(gps, gal), 44.0);
        assert_eq!(report.n_reference_satellites, 2);
        assert_eq!(report.antennas[0].n_sats, 2);
        assert_eq!(report.antennas[1].n_sats, 0);
        assert_eq!(report.antennas[0].median_cn0_db_hz, Some(44.5));
        assert_eq!(report.antennas[1].median_cn0_db_hz, Some(30.5));
    }

    #[test]
    fn test_run_every_antenna_below_threshold_excluded() {
        let gps = crate::acquisition::GpsAllAcquisitionOutput {
            antenna_numbers: vec![0, 1],
            results: vec![gps_row("GPS01", Some(vec![Some(40.0), Some(41.0)]))],
        };
        let gal = crate::galileo::GalileoAllAcquisitionOutput {
            antenna_numbers: vec![0, 1],
            results: vec![],
        };
        let report = run(&output(gps, gal), 44.0);
        assert_eq!(report.n_reference_satellites, 0);
        assert!(report.antennas.is_empty());
    }

    #[test]
    fn test_run_antenna_subset() {
        // antennas [1, 3] only (as if --ant 1,3 was given).
        let gps = crate::acquisition::GpsAllAcquisitionOutput {
            antenna_numbers: vec![1, 3],
            results: vec![
                gps_row("GPS01", Some(vec![Some(41.0), Some(42.0)])),
                gps_row("GPS02", Some(vec![Some(42.0), Some(40.0)])),
            ],
        };
        let gal = crate::galileo::GalileoAllAcquisitionOutput {
            antenna_numbers: vec![1, 3],
            results: vec![],
        };
        let report = run(&output(gps, gal), 40.0);
        assert_eq!(report.antenna_numbers, vec![1, 3]);
        assert_eq!(report.antennas.len(), 2);
        assert_eq!(report.antennas[0].antenna, 1);
        assert_eq!(report.antennas[1].antenna, 3);
        assert_eq!(report.reference_satellites, vec!["GPS01", "GPS02"]);
    }

    #[test]
    fn test_run_multiconstellation_collection() {
        let gps = crate::acquisition::GpsAllAcquisitionOutput {
            antenna_numbers: vec![0, 1],
            results: vec![gps_row("GPS10", Some(vec![Some(42.0), Some(43.0)]))],
        };
        let gal = crate::galileo::GalileoAllAcquisitionOutput {
            antenna_numbers: vec![0, 1],
            results: vec![gal_row("GSAT05", Some(vec![Some(44.0), Some(43.5)]))],
        };
        let report = run(&output(gps, gal), 40.0);
        assert_eq!(report.n_reference_satellites, 2);
        assert_eq!(
            report.reference_satellites,
            vec!["GPS10", "GSAT05"]
        );
        // ant0 median = (42+44)/2 = 43, ant1 = (43+43.5)/2 = 43.25 → ant1 best.
        assert_eq!(report.antennas[0].antenna, 0);
        assert_eq!(report.antennas[1].antenna, 1);
        assert_eq!(report.antennas[1].offset_db, Some(0.0));
        assert_eq!(report.antennas[0].offset_db, Some(-0.25));
    }

    // -----------------------------------------------------------------------
    // End-to-end: per-antenna gain offsets injected into synthesized GPS
    // signals, recovered by acquisition + antenna_test.
    // -----------------------------------------------------------------------

    const FS: f64 = 16.368e6;
    const SAMPLES_PER_MS: usize = 16_368;
    // Slightly offset from GPS_IF so the 16.368 MHz sampling grid does not
    // land exactly on carrier zeros (which would make 1-bit quantization
    // immune to noise). Still inside the ±6000 Hz search band.
    const CARRIER: f64 = GPS_IF - 2000.0;

    /// Small deterministic PRNG for the injected noise (kept local so this
    /// test module stays independent of the crate's private testutil).
    struct XorShift(u64);
    impl XorShift {
        fn new(seed: u64) -> Self {
            Self(seed | 1)
        }
        fn next_f64(&mut self) -> f64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            (x >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    fn upsample_ca(prn: usize) -> Vec<f64> {
        let code: [f64; CA_CHIPS] = generate_ca_code(prn);
        (0..SAMPLES_PER_MS).map(|i| code[i / 16]).collect()
    }

    /// Observed data for `n_ant` antennas, each receiving `gains[ant] ×` the
    /// GPS signals for `prns` (all at GPS IF, distinct delays), plus per-antenna
    /// uniform noise of amplitude `noise_levels[ant]`.
    fn observation_with_gains(
        prns: &[usize],
        n_ant: usize,
        gains: &[f64],
        noise_levels: &[f64],
    ) -> Observation {
        let delays: Vec<usize> = prns.iter().enumerate().map(|(k, _)| 137 * (k + 1)).collect();
        let codes: Vec<Vec<f64>> = prns.iter().map(|&p| upsample_ca(p)).collect();
        let mut data = Vec::with_capacity(n_ant);
        for j in 0..n_ant {
            let g = gains[j];
            let noise_amp = noise_levels[j.min(noise_levels.len() - 1)];
            let mut rng = XorShift::new(1000 + j as u64);
            let mut signal = Vec::with_capacity(2 * SAMPLES_PER_MS);
            for i in 0..2 * SAMPLES_PER_MS {
                let t = i as f64 / FS;
                let mut v = 0.0f64;
                for (k, code) in codes.iter().enumerate() {
                    // 1 ms code period; second ms repeats the first.
                    let cp = (i + delays[k]) % SAMPLES_PER_MS;
                    v += code[cp] * (TAU * CARRIER * t).cos();
                }
                v *= g;
                v += noise_amp * (2.0 * rng.next_f64() - 1.0);
                signal.push(v);
            }
            data.push(quantize(&signal));
        }
        let config = Config::from_json(&format!(
            r#"{{"num_antenna":{n_ant},"sampling_frequency":{FS}}}"#
        ))
        .expect("hard-coded config JSON should be valid");
        Observation::new(chrono::Utc::now(), config, data)
    }

    #[test]
    fn test_end_to_end_ranking_recovers_injected_gains() {
        // One-bit quantization (sign only) erases pure amplitude differences,
        // so the surviving lever for a "worse radio" is extra noise — antenna 2
        // gets 8x the noise of antennas 0/1. With these parameters the two
        // good antennas sit near 56 dB-Hz and the degraded one near 44 dB-Hz.
        let n_ant = 3;
        let prns: &[usize] = &[1, 7, 20];
        let obs = observation_with_gains(prns, n_ant, &[1.0; 3], &[2.0, 2.0, 16.0]);

        let acq = acquire_all_gps(&obs, GPS_IF, GPS_SEARCH_BAND, None, Some(prns), false, true);
        let report = run(
            &CombinedOutput {
                gps: Some(acq),
                galileo: None,
                beidou: None,
                sbas: None,
                l1c: None,
                qzss: None,
            },
            40.0,
        );

        assert_eq!(report.antenna_numbers, vec![0, 1, 2]);
        assert_eq!(report.reference_satellites.len(), prns.len());
        assert!(report.n_reference_satellites == prns.len(), "no reference satellites");
        let st = &report.antennas;
        // Best antenna must be 0 or 1 (quiet), never the noisy one.
        assert!(st[0].antenna != 2, "degraded antenna ranked best");
        // The noisy antenna must rank last with a clearly negative offset.
        let weak = st.iter().find(|s| s.antenna == 2).unwrap();
        assert_eq!(weak.rank, Some(n_ant));
        assert!(weak.offset_db.unwrap() < -5.0, "weak offset {:.2}", weak.offset_db.unwrap());
        // The two quiet antennas must be within 1 dB of each other, and both
        // well above the noisy antenna.
        let strong: Vec<&AntennaStat> = st.iter().filter(|s| s.antenna != 2).collect();
        let m0 = strong[0].median_cn0_db_hz.unwrap();
        let m1 = strong[1].median_cn0_db_hz.unwrap();
        assert!((m0 - m1).abs() < 1.0);
        let wm = weak.median_cn0_db_hz.unwrap();
        assert!(
            m0 - wm > 5.0,
            "median gap too small: strong {:.1} vs weak {:.1}",
            m0,
            wm
        );
    }
}
