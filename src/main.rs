// Copyright (c) 2026 Tim Molteno <tim@elec.ac.nz>
// SPDX-License-Identifier: GPL-3.0

mod acr;
mod acquisition;
mod beidou;
mod config;
mod correlate;
mod galileo;
mod l1c;
mod observation;
mod qzss;
mod sbas;
mod stats;

#[cfg(test)]
mod testutil;

use observation::Observation;
use serde::Serialize;

#[derive(Serialize)]
struct CombinedOutput {
    #[serde(skip_serializing_if = "Option::is_none")]
    gps: Option<acquisition::GpsAllAcquisitionOutput>,
    #[serde(skip_serializing_if = "Option::is_none")]
    galileo: Option<galileo::GalileoAllAcquisitionOutput>,
    #[serde(skip_serializing_if = "Option::is_none")]
    beidou: Option<beidou::BeiDouAllAcquisitionOutput>,
    #[serde(skip_serializing_if = "Option::is_none")]
    sbas: Option<sbas::SbasAllAcquisitionOutput>,
    #[serde(skip_serializing_if = "Option::is_none")]
    l1c: Option<l1c::L1CAllAcquisitionOutput>,
    #[serde(skip_serializing_if = "Option::is_none")]
    qzss: Option<qzss::QzssAllAcquisitionOutput>,
}

fn main() {
    let args: Vec<String> = std::env::args().collect();

    let mut file: Option<String> = None;
    let mut antenna_i: Option<usize> = None;
    let mut antenna_j: Option<usize> = None;
    let mut gps_flag = false;
    let mut galileo_flag = false;
    let mut beidou_flag = false;
    let mut sbas_flag = false;
    let mut l1c_flag = false;
    let mut qzss_flag = false;
    let mut all_flag = false;
    let mut ant_list: Option<Vec<usize>> = None;
    let mut filter_phase_mad: Option<f64> = None;
    let mut filter_freq_mad: Option<f64> = None;
    let mut output_file: Option<String> = None;
    let mut debug_flag = false;
    let mut cn0_flag = false;
    let mut benchmark_flag = false;
    let mut prn_filter: Option<Vec<usize>> = None;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--file" => {
                i += 1;
                file = Some(args[i].clone());
            }
            "--i" => {
                i += 1;
                antenna_i = Some(args[i].parse().expect("invalid integer for --i"));
            }
            "--j" => {
                i += 1;
                antenna_j = Some(args[i].parse().expect("invalid integer for --j"));
            }
            "--gps" => {
                gps_flag = true;
            }
            "--galileo" => {
                galileo_flag = true;
            }
            "--beidou" => {
                beidou_flag = true;
            }
            "--sbas" => {
                sbas_flag = true;
            }
            "--l1c" => {
                l1c_flag = true;
            }
            "--qzss" => {
                qzss_flag = true;
            }
            "--all" => {
                all_flag = true;
            }
            "--ant" => {
                i += 1;
                ant_list = Some(
                    args[i]
                        .split(',')
                        .map(|s| s.trim().parse().expect("invalid integer in --ant list"))
                        .collect(),
                );
            }
            "--filter-phase-mad" => {
                i += 1;
                filter_phase_mad = Some(args[i].parse().expect("invalid float for --filter-phase-mad"));
            }
            "--filter-freq-mad" => {
                i += 1;
                filter_freq_mad = Some(args[i].parse().expect("invalid float for --filter-freq-mad"));
            }
            "--output" => {
                i += 1;
                output_file = Some(args[i].clone());
            }
            "--debug" => {
                debug_flag = true;
            }
            "--cn0" => {
                cn0_flag = true;
            }
            "--benchmark" => {
                benchmark_flag = true;
            }
            "--prn" => {
                i += 1;
                prn_filter = Some(
                    args[i]
                        .split(',')
                        .map(|s| s.trim().parse().expect("invalid integer in --prn list"))
                        .collect(),
                );
            }
            "--version" => {
                println!("tart-gnss-acquire v{}", env!("CARGO_PKG_VERSION"));
                std::process::exit(0);
            }
            other => {
                eprintln!("unknown argument: {other}");
                std::process::exit(1);
            }
        }
        i += 1;
    }

    // --all implies all six acquisition modes
    if all_flag {
        gps_flag = true;
        galileo_flag = true;
        beidou_flag = true;
        sbas_flag = true;
        l1c_flag = true;
        qzss_flag = true;
    }

    let file = if benchmark_flag {
        file
    } else {
        Some(file.unwrap_or_else(|| {
            eprintln!(
                "usage: {} --file <observation.hdf> [--i <i> --j <j>] [--all] [--gps] [--galileo] [--beidou] [--sbas] [--l1c] [--qzss] [--cn0] [--prn <a,b,...>] [--ant <a,b,...>] [--filter-phase-mad <x>] [--filter-freq-mad <x>] [--output <path>] [--debug] [--benchmark]",
                args[0]
            );
            std::process::exit(1);
        }))
    };

    let obs = if let Some(ref f) = file {
        Observation::from_hdf5(f).unwrap_or_else(|e| {
            eprintln!("error loading HDF5 file: {e}");
            std::process::exit(1);
        })
    } else {
        // Benchmark mode without a data file — generate random noise.
        // 500k samples per antenna is enough for 31 ms at 16 MHz, covering
        // the longest code period (L1C/BeiDou = 20 ms with 2 epochs).
        eprintln!("--benchmark without --file: using random data");
        Observation::random(2, 16e6, 500_000)
    };

    let n_ant = obs.config.num_antenna();
    let sampling_freq = obs.get_sampling_rate();
    eprintln!("timestamp:   {}", obs.timestamp);
    eprintln!("antennas:    {n_ant}");
    eprintln!("sample rate: {sampling_freq} Hz");

    // --- Benchmark mode ---------------------------------------------------
    if benchmark_flag {
        use std::f64::consts::TAU;
        use std::time::Instant;

        #[derive(Serialize)]
        struct BenchResult {
            constellation: &'static str,
            prns: usize,
            setup_s: f64,
            search_s: f64,
            duration_s: f64,
        }

        /// Extract de-meaned f32 antenna data and pre-computed phasepoints
        /// for the given number of samples.  This is the shared setup work
        /// that all acquisition pipelines perform before the per-PRN search.
        fn extract_ant_data(
            obs: &Observation,
            num_samples: usize,
        ) -> Vec<(Vec<f32>, Vec<f64>)> {
            let sampling_freq = obs.get_sampling_rate();
            let n_ant = obs.config.num_antenna();
            let phase_const = TAU / sampling_freq;
            (0..n_ant)
                .map(|ant_idx| {
                    let bipolar = obs.get_antenna(ant_idx);
                    let mean = bipolar.iter().sum::<f64>() / bipolar.len() as f64;
                    let raw: Vec<f64> = bipolar.iter().map(|&v| v - mean).collect();
                    let len = num_samples.min(raw.len());
                    let signal_f32: Vec<f32> =
                        raw[..len].iter().map(|&v| v as f32).collect();
                    let phasepoints: Vec<f64> =
                        (0..len).map(|i| phase_const * i as f64).collect();
                    (signal_f32, phasepoints)
                })
                .collect()
        }

        let sampling_freq = obs.get_sampling_rate();
        let samples_per_ms = sampling_freq / 1000.0;
        let num_samples_per_ms = samples_per_ms as usize;

        // Galileo / BeiDou / L1C code periods (chip count / chip rate).
        const GAL_CODE_PERIOD: f64 = 4092.0 / 1.023e6; // 4 ms
        const BDS_CODE_PERIOD: f64 = 10230.0 / 1.023e6; // 10 ms
        const L1C_CODE_PERIOD: f64 = 10230.0 / 1.023e6; // 10 ms

        let gal_samples_per_period =
            (sampling_freq * GAL_CODE_PERIOD) as usize;
        let bds_samples_per_period =
            (sampling_freq * BDS_CODE_PERIOD) as usize;
        let l1c_samples_per_period =
            (sampling_freq * L1C_CODE_PERIOD) as usize;

        struct BenchEntry {
            name: &'static str,
            num_samples: usize,      // for data extraction (setup)
            prn: usize,
            center_freq: f64,
            search_band: f64,
            samples_per_chunk: usize, // passed to single-PRN function
        }

        let entries: [BenchEntry; 6] = [
            BenchEntry {
                name: "GPS",
                num_samples: 2 * num_samples_per_ms,
                prn: 1,
                center_freq: acquisition::GPS_IF,
                search_band: acquisition::GPS_SEARCH_BAND,
                samples_per_chunk: num_samples_per_ms,
            },
            BenchEntry {
                name: "Galileo",
                num_samples: 2 * gal_samples_per_period,
                prn: 1,
                center_freq: 4.092e6,
                search_band: 6000.0,
                samples_per_chunk: gal_samples_per_period,
            },
            BenchEntry {
                name: "BeiDou",
                num_samples: 2 * bds_samples_per_period,
                prn: 1,
                center_freq: 4.092e6,
                search_band: 6000.0,
                samples_per_chunk: bds_samples_per_period,
            },
            BenchEntry {
                name: "SBAS",
                num_samples: 2 * num_samples_per_ms,
                prn: 120,
                center_freq: sbas::SBAS_IF,
                search_band: sbas::SBAS_SEARCH_BAND,
                samples_per_chunk: num_samples_per_ms,
            },
            BenchEntry {
                name: "QZSS",
                num_samples: 2 * num_samples_per_ms,
                prn: 184,
                center_freq: qzss::QZSS_IF,
                search_band: qzss::QZSS_SEARCH_BAND,
                samples_per_chunk: num_samples_per_ms,
            },
            BenchEntry {
                name: "L1C",
                num_samples: 2 * l1c_samples_per_period,
                prn: 1,
                center_freq: l1c::L1C_IF,
                search_band: l1c::L1C_SEARCH_BAND,
                samples_per_chunk: l1c_samples_per_period,
            },
        ];

        let mut results: Vec<BenchResult> = Vec::with_capacity(6);
        let mut total_setup: f64 = 0.0;
        let mut total_search: f64 = 0.0;

        eprintln!("Running benchmark (one PRN per constellation, setup timed separately)...");

        for entry in &entries {
            eprintln!("  {} (PRN {})...", entry.name, entry.prn);

            // --- setup ----------------------------------------------------
            let setup_start = Instant::now();
            let ant_data = extract_ant_data(&obs, entry.num_samples);
            let setup_s = setup_start.elapsed().as_secs_f64();

            // --- single-PRN search (antenna 0 only) -----------------------
            let (signal_f32, phasepoints) = &ant_data[0];
            let search_start = Instant::now();
            match entry.name {
                "GPS" => {
                    acquisition::acquire_full(
                        signal_f32, phasepoints, sampling_freq,
                        entry.center_freq, entry.search_band,
                        entry.prn, entry.samples_per_chunk,
                    );
                }
                "Galileo" => {
                    galileo::acquire_galileo_single(
                        signal_f32, phasepoints, sampling_freq,
                        entry.center_freq, entry.search_band,
                        entry.prn, entry.samples_per_chunk,
                    );
                }
                "BeiDou" => {
                    beidou::acquire_beidou_single(
                        signal_f32, phasepoints, sampling_freq,
                        entry.center_freq, entry.search_band,
                        entry.prn, entry.samples_per_chunk,
                    );
                }
                "SBAS" => {
                    sbas::acquire_sbas(
                        signal_f32, phasepoints, sampling_freq,
                        entry.center_freq, entry.search_band,
                        entry.prn, entry.samples_per_chunk,
                    );
                }
                "QZSS" => {
                    qzss::acquire_qzss(
                        signal_f32, phasepoints, sampling_freq,
                        entry.center_freq, entry.search_band,
                        entry.prn, entry.samples_per_chunk,
                    );
                }
                "L1C" => {
                    l1c::acquire_l1c_single(
                        signal_f32, phasepoints, sampling_freq,
                        entry.center_freq, entry.search_band,
                        entry.prn, entry.samples_per_chunk,
                    );
                }
                _ => unreachable!(),
            }
            let search_s = search_start.elapsed().as_secs_f64();
            let duration_s = setup_s + search_s;

            eprintln!(
                "    {}: setup={setup_s:.4}s  search={search_s:.4}s  total={duration_s:.4}s",
                entry.name
            );
            total_setup += setup_s;
            total_search += search_s;
            results.push(BenchResult {
                constellation: entry.name,
                prns: 1,
                setup_s,
                search_s,
                duration_s,
            });
        }

        let total_prns: usize = results.iter().map(|r| r.prns).sum();
        let total_duration_s = total_setup + total_search;
        let acq_per_s = total_prns as f64 / total_search;
        eprintln!(
            "  Combined: startup={total_setup:.4}s  search={total_search:.4}s  total={total_duration_s:.4}s  acq/s={acq_per_s:.1}"
        );

        let bench_output = serde_json::to_string_pretty(&results).expect("JSON serialization failed");
        if let Some(path) = &output_file {
            std::fs::write(path, &bench_output).unwrap_or_else(|e| {
                eprintln!("error writing output file {path}: {e}");
                std::process::exit(1);
            });
            eprintln!("Wrote benchmark JSON to {path}");
        } else {
            println!("{bench_output}");
        }
        return;
    }

    let any_acq = gps_flag || galileo_flag || beidou_flag || sbas_flag || l1c_flag || qzss_flag;

    if any_acq {
        // Validate --ant indices once for all acquisition modes
        if let Some(ref ants) = ant_list {
            for &ant in ants {
                if ant >= n_ant {
                    eprintln!("antenna index {ant} out of range (have {n_ant})");
                    std::process::exit(1);
                }
            }
        }

        let mut output = CombinedOutput {
            gps: None,
            galileo: None,
            beidou: None,
            sbas: None,
            l1c: None,
            qzss: None,
        };

        // --- GPS -----------------------------------------------------------
        if gps_flag {
            let n_prns = prn_filter.as_ref().map_or(acquisition::GPS_NUM_SATS, |f| f.len());
            eprintln!(
                "Running GPS L1 C/A search ({} PRNs)...",
                n_prns
            );
            output.gps = Some(acquisition::acquire_all_gps(
                &obs,
                acquisition::GPS_IF,
                acquisition::GPS_SEARCH_BAND,
                ant_list.clone(),
                prn_filter.as_deref(),
                debug_flag,
                cn0_flag,
            ));
        }

        // --- Galileo -------------------------------------------------------
        if galileo_flag {
            let galileo_if = 4.092e6;
            let galileo_search_band = 6000.0;
            let n_prns = prn_filter.as_ref().map_or(50, |f| f.len());
            eprintln!("Running Galileo E1-C search ({} PRNs)...", n_prns);
            output.galileo = Some(galileo::acquire_all_galileo(
                &obs,
                galileo_if,
                galileo_search_band,
                ant_list.clone(),
                prn_filter.as_deref(),
                debug_flag,
                cn0_flag,
            ));
        }

        // --- BeiDou --------------------------------------------------------
        if beidou_flag {
            let beidou_if = 4.092e6;
            let beidou_search_band = 6000.0;
            let n_prns = prn_filter.as_ref().map_or(63, |f| f.len());
            eprintln!("Running BeiDou B1C search ({} PRNs)...", n_prns);
            output.beidou = Some(beidou::acquire_all_beidou(
                &obs,
                beidou_if,
                beidou_search_band,
                ant_list.clone(),
                prn_filter.as_deref(),
                debug_flag,
                cn0_flag,
            ));
        }

        // --- SBAS ----------------------------------------------------------
        if sbas_flag {
            let n_prns = prn_filter.as_ref().map_or(sbas::SBAS_NUM_SATS, |f| f.len());
            eprintln!(
                "Running SBAS L1 C/A search ({} PRNs)...",
                n_prns
            );
            output.sbas = Some(sbas::acquire_all_sbas(
                &obs,
                sbas::SBAS_IF,
                sbas::SBAS_SEARCH_BAND,
                ant_list.clone(),
                prn_filter.as_deref(),
                debug_flag,
                cn0_flag,
            ));
        }

        // --- GPS L1C -------------------------------------------------------
        if l1c_flag {
            let n_prns = prn_filter.as_ref().map_or(l1c::L1C_NUM_SATS, |f| f.len());
            eprintln!("Running GPS L1C search ({} PRNs)...", n_prns);
            output.l1c = Some(l1c::acquire_all_l1c(
                &obs,
                l1c::L1C_IF,
                l1c::L1C_SEARCH_BAND,
                ant_list.clone(),
                prn_filter.as_deref(),
                debug_flag,
                cn0_flag,
            ));
        }

        // --- QZSS ----------------------------------------------------------
        if qzss_flag {
            let n_prns = prn_filter.as_ref().map_or(qzss::QZSS_NUM_SATS, |f| f.len());
            eprintln!(
                "Running QZSS L1 C/A search ({} PRNs)...",
                n_prns
            );
            output.qzss = Some(qzss::acquire_all_qzss(
                &obs,
                qzss::QZSS_IF,
                qzss::QZSS_SEARCH_BAND,
                ant_list.clone(),
                prn_filter.as_deref(),
                debug_flag,
                cn0_flag,
            ));
        }

        // --- Apply MAD filters ---------------------------------------------
        if filter_phase_mad.is_some() || filter_freq_mad.is_some() {
            let mut filter_count = 0u64;

            if let Some(ref mut gps_out) = output.gps {
                let before = gps_out.results.len();
                if let Some(thresh) = filter_phase_mad {
                    gps_out.results.retain(|r| r.phase_mad.map_or(true, |m| m <= thresh));
                }
                if let Some(thresh) = filter_freq_mad {
                    gps_out.results.retain(|r| r.freq_mad.map_or(true, |m| m <= thresh));
                }
                filter_count += (before - gps_out.results.len()) as u64;
            }
            if let Some(ref mut gal_out) = output.galileo {
                let before = gal_out.results.len();
                if let Some(thresh) = filter_phase_mad {
                    gal_out.results.retain(|r| r.phase_mad.map_or(true, |m| m <= thresh));
                }
                if let Some(thresh) = filter_freq_mad {
                    gal_out.results.retain(|r| r.freq_mad.map_or(true, |m| m <= thresh));
                }
                filter_count += (before - gal_out.results.len()) as u64;
            }
            if let Some(ref mut bd_out) = output.beidou {
                let before = bd_out.results.len();
                if let Some(thresh) = filter_phase_mad {
                    bd_out.results.retain(|r| r.phase_mad.map_or(true, |m| m <= thresh));
                }
                if let Some(thresh) = filter_freq_mad {
                    bd_out.results.retain(|r| r.freq_mad.map_or(true, |m| m <= thresh));
                }
                filter_count += (before - bd_out.results.len()) as u64;
            }
            if let Some(ref mut sb_out) = output.sbas {
                let before = sb_out.results.len();
                if let Some(thresh) = filter_phase_mad {
                    sb_out.results.retain(|r| r.phase_mad.map_or(true, |m| m <= thresh));
                }
                if let Some(thresh) = filter_freq_mad {
                    sb_out.results.retain(|r| r.freq_mad.map_or(true, |m| m <= thresh));
                }
                filter_count += (before - sb_out.results.len()) as u64;
            }
            if let Some(ref mut l1c_out) = output.l1c {
                let before = l1c_out.results.len();
                if let Some(thresh) = filter_phase_mad {
                    l1c_out.results.retain(|r| r.phase_mad.map_or(true, |m| m <= thresh));
                }
                if let Some(thresh) = filter_freq_mad {
                    l1c_out.results.retain(|r| r.freq_mad.map_or(true, |m| m <= thresh));
                }
                filter_count += (before - l1c_out.results.len()) as u64;
            }
            if let Some(ref mut qz_out) = output.qzss {
                let before = qz_out.results.len();
                if let Some(thresh) = filter_phase_mad {
                    qz_out.results.retain(|r| r.phase_mad.map_or(true, |m| m <= thresh));
                }
                if let Some(thresh) = filter_freq_mad {
                    qz_out.results.retain(|r| r.freq_mad.map_or(true, |m| m <= thresh));
                }
                filter_count += (before - qz_out.results.len()) as u64;
            }

            if filter_count > 0 {
                eprintln!("Filtered out {filter_count} results (MAD thresholds)");
            }
        }

        let json = serde_json::to_string_pretty(&output).expect("JSON serialization failed");
        if let Some(path) = &output_file {
            std::fs::write(path, &json).unwrap_or_else(|e| {
                eprintln!("error writing output file {path}: {e}");
                std::process::exit(1);
            });
            eprintln!("Wrote JSON output to {path}");
        } else {
            println!("{json}");
        }
        return;
    }

    // --- Correlation mode --------------------------------------------------
    let i = antenna_i.unwrap_or_else(|| {
        eprintln!("missing --i <antenna_i> (or use --gps, --galileo, --beidou, --sbas, --l1c, --qzss, or --all)");
        std::process::exit(1);
    });
    let j = antenna_j.unwrap_or_else(|| {
        eprintln!("missing --j <antenna_j>");
        std::process::exit(1);
    });

    if i >= n_ant {
        eprintln!("antenna index {i} out of range (have {n_ant})");
        std::process::exit(1);
    }
    if j >= n_ant {
        eprintln!("antenna index {j} out of range (have {n_ant})");
        std::process::exit(1);
    }

    let correlation = obs.correlate(i, j);
    println!("antenna_i:   {i}");
    println!("antenna_j:   {j}");
    println!("correlation: {correlation:.6}");
}
