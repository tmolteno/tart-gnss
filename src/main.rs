// Copyright (c) 2026 Tim Molteno <tim@elec.ac.nz>
// SPDX-License-Identifier: GPL-3.0

use tart_gnss_acquire::{
    acquisition, antenna_test, beidou, galileo, l1c, observation, qzss, rfi, sbas, simulator,
    CombinedOutput,
};

use observation::Observation;
use serde::Serialize;

/// All CLI options parsed from the command line.
#[derive(Debug, Default)]
struct ParsedArgs {
    file: Option<String>,
    antenna_i: Option<usize>,
    antenna_j: Option<usize>,
    gps_flag: bool,
    galileo_flag: bool,
    beidou_flag: bool,
    sbas_flag: bool,
    l1c_flag: bool,
    qzss_flag: bool,
    all_flag: bool,
    ant_list: Option<Vec<usize>>,
    filter_phase_mad: Option<f64>,
    filter_freq_mad: Option<f64>,
    output_file: Option<String>,
    debug_flag: bool,
    cn0_flag: bool,
    benchmark_flag: bool,
    test_antennas_flag: bool,
    test_min_cn0: Option<f64>,
    prn_filter: Option<Vec<usize>>,
    version_flag: bool,
    rfi_flag: bool,
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
}

/// Read the value following a flag, advancing `i`, or return an error if absent.
fn flag_value<'a>(args: &'a [String], i: &mut usize, flag: &str) -> Result<&'a str, String> {
    if *i + 1 >= args.len() {
        return Err(format!("missing value for {flag}"));
    }
    *i += 1;
    Ok(args[*i].as_str())
}

/// Parse command-line args (excluding the program name in `args[0]`), matching
/// the original `main` loop. Returns an error string on bad input.
fn parse_args(args: &[String]) -> Result<ParsedArgs, String> {
    let mut p = ParsedArgs::default();
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--file" => {
                let v = flag_value(args, &mut i, "--file")?;
                p.file = Some(v.to_string());
            }
            "--i" => {
                let v = flag_value(args, &mut i, "--i")?;
                p.antenna_i = Some(v.parse().map_err(|_| format!("invalid integer for --i: {v}"))?);
            }
            "--j" => {
                let v = flag_value(args, &mut i, "--j")?;
                p.antenna_j = Some(v.parse().map_err(|_| format!("invalid integer for --j: {v}"))?);
            }
            "--gps" => p.gps_flag = true,
            "--galileo" => p.galileo_flag = true,
            "--beidou" => p.beidou_flag = true,
            "--sbas" => p.sbas_flag = true,
            "--l1c" => p.l1c_flag = true,
            "--qzss" => p.qzss_flag = true,
            "--all" => p.all_flag = true,
            "--ant" => {
                let v = flag_value(args, &mut i, "--ant")?;
                p.ant_list = Some(
                    v.split(',')
                        .map(|s| s.trim().parse().map_err(|_| format!("invalid integer in --ant: {s}")))
                        .collect::<Result<Vec<_>, _>>()?,
                );
            }
            "--filter-phase-mad" => {
                let v = flag_value(args, &mut i, "--filter-phase-mad")?;
                p.filter_phase_mad = Some(v.parse().map_err(|_| format!("invalid float for --filter-phase-mad: {v}"))?);
            }
            "--filter-freq-mad" => {
                let v = flag_value(args, &mut i, "--filter-freq-mad")?;
                p.filter_freq_mad = Some(v.parse().map_err(|_| format!("invalid float for --filter-freq-mad: {v}"))?);
            }
            "--output" => {
                let v = flag_value(args, &mut i, "--output")?;
                p.output_file = Some(v.to_string());
            }
            "--debug" => p.debug_flag = true,
            "--cn0" => p.cn0_flag = true,
            "--benchmark" => p.benchmark_flag = true,
            "--test-antennas" => p.test_antennas_flag = true,
            "--test-min-cn0" => {
                let v = flag_value(args, &mut i, "--test-min-cn0")?;
                p.test_min_cn0 = Some(
                    v.parse().map_err(|_| format!("invalid float for --test-min-cn0: {v}"))?,
                );
            }
            "--prn" => {
                let v = flag_value(args, &mut i, "--prn")?;
                p.prn_filter = Some(
                    v.split(',')
                        .map(|s| s.trim().parse().map_err(|_| format!("invalid integer in --prn: {s}")))
                        .collect::<Result<Vec<_>, _>>()?,
                );
            }
            "--version" => p.version_flag = true,
            "--rfi" => p.rfi_flag = true,
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
            other => return Err(format!("unknown argument: {other}")),
        }
        i += 1;
    }
    Ok(p)
}

/// Resolve flag implications (`--all` ⇒ six constellations; `--test-antennas`
/// ⇒ force `--cn0` and, unless a specific constellation flag is given, imply
/// `--all`). Returns an error string for incompatible combinations.
fn apply_mode_implications(p: &mut ParsedArgs) -> Result<(), String> {
    if p.all_flag {
        p.gps_flag = true;
        p.galileo_flag = true;
        p.beidou_flag = true;
        p.sbas_flag = true;
        p.l1c_flag = true;
        p.qzss_flag = true;
    }
    if p.test_antennas_flag {
        p.cn0_flag = true;
        if p.benchmark_flag {
            return Err("--test-antennas cannot be combined with --benchmark".to_string());
        }
        if !(p.gps_flag || p.galileo_flag || p.beidou_flag || p.sbas_flag || p.l1c_flag || p.qzss_flag) {
            p.gps_flag = true;
            p.galileo_flag = true;
            p.beidou_flag = true;
            p.sbas_flag = true;
            p.l1c_flag = true;
            p.qzss_flag = true;
        }
    }
    if p.rfi_flag {
        if p.benchmark_flag {
            return Err("--rfi cannot be combined with --benchmark".to_string());
        }
        if p.test_antennas_flag {
            return Err("--rfi cannot be combined with --test-antennas".to_string());
        }
        if p.antenna_i.is_some() || p.antenna_j.is_some() {
            return Err("--rfi cannot be combined with --i/--j".to_string());
        }
        if p.simulate {
            return Err("--rfi cannot be combined with --simulate".to_string());
        }
    }
    Ok(())
}

fn print_usage(prog: &str) {
    eprintln!(
        "usage: {prog} --file <observation.hdf> [--i <i> --j <j>] [--all] [--gps] [--galileo] [--beidou] [--sbas] [--l1c] [--qzss] [--cn0] [--test-antennas] [--test-min-cn0 <x>] [--rfi] [--prn <a,b,...>] [--ant <a,b,...>] [--filter-phase-mad <x>] [--filter-freq-mad <x>] [--output <path>] [--debug] [--benchmark]"
    );
}

/// Common accessors for MAD-based filtering across PRN result types.
trait MadFilterable {
    fn phase_mad(&self) -> Option<f64>;
    fn freq_mad(&self) -> Option<f64>;
}

macro_rules! impl_mad_filterable {
    ($ty:ty) => {
        impl MadFilterable for $ty {
            fn phase_mad(&self) -> Option<f64> { self.phase_mad }
            fn freq_mad(&self) -> Option<f64> { self.freq_mad }
        }
    };
}

impl_mad_filterable!(acquisition::GpsPrnResult);
impl_mad_filterable!(galileo::GalileoPrnResult);
impl_mad_filterable!(beidou::BeiDouPrnResult);
impl_mad_filterable!(sbas::SbasPrnResult);
impl_mad_filterable!(l1c::L1CPrnResult);
impl_mad_filterable!(qzss::QzssPrnResult);

/// Retain only results whose MADs are at or below the given thresholds.
/// A `None` MAD (single antenna) always passes. Returns the number removed.
fn apply_mad_filters<T: MadFilterable>(
    results: &mut Vec<T>,
    phase_thresh: Option<f64>,
    freq_thresh: Option<f64>,
) -> usize {
    let before = results.len();
    if let Some(t) = phase_thresh {
        results.retain(|r| r.phase_mad().map_or(true, |m| m <= t));
    }
    if let Some(t) = freq_thresh {
        results.retain(|r| r.freq_mad().map_or(true, |m| m <= t));
    }
    before - results.len()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut parsed = match parse_args(&args) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("{e}");
            print_usage(&args[0]);
            std::process::exit(1);
        }
    };

    if parsed.version_flag {
        println!("tart-gnss-acquire v{}", env!("CARGO_PKG_VERSION"));
        std::process::exit(0);
    }

    // --all / --test-antennas flag implications, before destructuring.
    if let Err(e) = apply_mode_implications(&mut parsed) {
        eprintln!("{e}");
        print_usage(&args[0]);
        std::process::exit(1);
    }

    let ParsedArgs {
        file,
        antenna_i,
        antenna_j,
        gps_flag,
        galileo_flag,
        beidou_flag,
        sbas_flag,
        l1c_flag,
        qzss_flag,
        all_flag: _,
        ant_list,
        filter_phase_mad,
        filter_freq_mad,
        output_file,
        debug_flag,
        cn0_flag,
        benchmark_flag,
        test_antennas_flag,
        test_min_cn0,
        prn_filter,
        version_flag: _,
        rfi_flag,
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
    } = parsed;

    // --- Simulation mode -------------------------------------------------
    if simulate {
        use simulator::{
            pack_row, parse_positions, parse_sources, random_positions, synthesize,
            write_observation, SimConfig, DEFAULT_SEED,
        };

        let src_path = sources.unwrap_or_else(|| {
            eprintln!("--simulate requires --sources <catalogue.json>");
            print_usage(&args[0]);
            std::process::exit(1);
        });
        let src_json = std::fs::read_to_string(&src_path).unwrap_or_else(|e| {
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
            let pos_json = std::fs::read_to_string(&pos_path).unwrap_or_else(|e| {
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
            let f = cfg.center_freq
                + cfg.band * (k as f64 / (src_list.len() as f64 - 1.0) - 0.5);
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

    let file = if benchmark_flag || simulate {
        file
    } else {
        Some(file.unwrap_or_else(|| {
            eprintln!("missing --file (or use --benchmark)");
            print_usage(&args[0]);
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

    // --- RFI report mode ---------------------------------------------------
    if rfi_flag {
        let report = rfi::run(&obs);
        for ch in &report.channels {
            if ch.dead {
                eprintln!("  antenna {:2}: dead channel (no samples)", ch.antenna);
                continue;
            }
            let acf = |k: Option<f64>| k.map(|v| format!("{v:+.3}")).unwrap_or_else(|| "  n/a".to_string());
            let period = ch
                .period_samples
                .map(|p| format!("{p:.2}"))
                .unwrap_or_else(|| "-".to_string());
            let lines: Vec<String> = ch
                .lines
                .iter()
                .map(|l| format!("{:.0} Hz ({:+.1} dB)", l.frequency, l.excess_db))
                .collect();
            eprintln!(
                "  antenna {:2}: acf(1)={} acf(3)={} acf(16)={}  period={} samples",
                ch.antenna,
                acf(ch.acf_lag1),
                acf(ch.acf_lag3),
                acf(ch.acf_lag16),
                period
            );
            if !lines.is_empty() {
                eprintln!("                lines: {}", lines.join("; "));
            }
        }
        let json = serde_json::to_string_pretty(&report).expect("JSON serialization failed");
        if let Some(path) = &output_file {
            std::fs::write(path, &json).unwrap_or_else(|e| {
                eprintln!("error writing output file {path}: {e}");
                std::process::exit(1);
            });
            eprintln!("Wrote RFI JSON to {path}");
        } else {
            println!("{json}");
        }
        return;
    }

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
            if let Some(o) = output.gps.as_mut() {
                filter_count += apply_mad_filters(&mut o.results, filter_phase_mad, filter_freq_mad) as u64;
            }
            if let Some(o) = output.galileo.as_mut() {
                filter_count += apply_mad_filters(&mut o.results, filter_phase_mad, filter_freq_mad) as u64;
            }
            if let Some(o) = output.beidou.as_mut() {
                filter_count += apply_mad_filters(&mut o.results, filter_phase_mad, filter_freq_mad) as u64;
            }
            if let Some(o) = output.sbas.as_mut() {
                filter_count += apply_mad_filters(&mut o.results, filter_phase_mad, filter_freq_mad) as u64;
            }
            if let Some(o) = output.l1c.as_mut() {
                filter_count += apply_mad_filters(&mut o.results, filter_phase_mad, filter_freq_mad) as u64;
            }
            if let Some(o) = output.qzss.as_mut() {
                filter_count += apply_mad_filters(&mut o.results, filter_phase_mad, filter_freq_mad) as u64;
            }
            if filter_count > 0 {
                eprintln!("Filtered out {filter_count} results (MAD thresholds)");
            }
        }

        // --- Antenna test mode ---------------------------------------------
        if test_antennas_flag {
            let min_cn0 = test_min_cn0.unwrap_or(44.0);
            let report = antenna_test::run(&output, min_cn0);
            if report.n_reference_satellites == 0 {
                eprintln!(
                    "warning: --test-antennas found no satellite with C/N0 ≥ {min_cn0:.1} dB-Hz on at least one antenna"
                );
            } else {
                eprintln!(
                    "Antenna relative C/N0 report: {} reference satellite(s) @ ≥ {min_cn0:.1} dB-Hz",
                    report.n_reference_satellites
                );
                for st in &report.antennas {
                    let med = st
                        .median_cn0_db_hz
                        .map(|m| format!("{m:5.1}"))
                        .unwrap_or_else(|| "  n/a".to_string());
                    let off = st
                        .offset_db
                        .map(|m| format!("{m:+.1}"))
                        .unwrap_or_else(|| " n/a".to_string());
                    let rk = st.rank.map(|r| r.to_string()).unwrap_or_else(|| "-".to_string());
                    eprintln!(
                        "  antenna {:2}: median C/N0 {} dB-Hz  offset {} dB  rank {}",
                        st.antenna, med, off, rk
                    );
                }
                eprintln!(
                    "  reference SVs: {}",
                    report.reference_satellites.join(", ")
                );
            }
            let json = serde_json::to_string_pretty(&report).expect("JSON serialization failed");
            if let Some(path) = &output_file {
                std::fs::write(path, &json).unwrap_or_else(|e| {
                    eprintln!("error writing output file {path}: {e}");
                    std::process::exit(1);
                });
                eprintln!("Wrote antenna-test JSON to {path}");
            } else {
                println!("{json}");
            }
            return;
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

#[cfg(test)]
mod tests {
    use super::*;

    fn a(list: &[&str]) -> Vec<String> {
        list.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn test_parse_all_flags() {
        let p = parse_args(&a(&[
            "tart-gnss-acquire", "--file", "obs.hdf", "--all", "--cn0", "--debug",
            "--prn", "1,5,12", "--ant", "0,2",
            "--filter-phase-mad", "0.002",
            "--filter-freq-mad", "100.0",
            "--output", "out.json",
        ]))
        .unwrap();
        assert_eq!(p.file.as_deref(), Some("obs.hdf"));
        assert_eq!(p.prn_filter.as_deref(), Some(&[1usize, 5, 12][..]));
        assert_eq!(p.ant_list.as_deref(), Some(&[0usize, 2][..]));
        assert_eq!(p.filter_phase_mad, Some(0.002));
        assert_eq!(p.filter_freq_mad, Some(100.0));
        assert_eq!(p.output_file.as_deref(), Some("out.json"));
        assert!(p.cn0_flag && p.debug_flag);
    }

    #[test]
    fn test_parse_test_antennas_flags() {
        let p = parse_args(&a(&[
            "prog", "--file", "obs.hdf", "--test-antennas",
            "--test-min-cn0", "42.5",
        ]))
        .unwrap();
        assert!(p.test_antennas_flag);
        assert_eq!(p.test_min_cn0, Some(42.5));
        // --test-antennas alone must imply cn0 + all constellations.
        let mut p2 = parse_args(&a(&["prog", "--file", "obs.hdf", "--test-antennas"])).unwrap();
        apply_mode_implications(&mut p2).unwrap();
        assert!(p2.cn0_flag);
        assert!(p2.gps_flag && p2.galileo_flag && p2.beidou_flag);
        assert!(p2.sbas_flag && p2.l1c_flag && p2.qzss_flag);
    }

    #[test]
    fn test_apply_mode_implications_respects_explicit_constellation() {
        let mut p = parse_args(&a(&[
            "prog", "--file", "obs.hdf", "--test-antennas", "--gps",
        ]))
        .unwrap();
        apply_mode_implications(&mut p).unwrap();
        assert!(p.cn0_flag);
        assert!(p.gps_flag);
        assert!(!p.galileo_flag && !p.qzss_flag, "must respect explicit --gps");
    }

    #[test]
    fn test_apply_mode_implications_benchmark_conflict() {
        let mut p = parse_args(&a(&["prog", "--test-antennas", "--benchmark"])).unwrap();
        assert!(apply_mode_implications(&mut p).is_err());
    }

    #[test]
    fn test_parse_all_sets_flags() {
        let p = parse_args(&a(&["prog", "--all", "--file", "x.hdf"])).unwrap();
        assert!(p.all_flag);
    }

    #[test]
    fn test_parse_unknown_arg_err() {
        assert!(parse_args(&a(&["prog", "--nope"])).is_err());
        assert!(parse_args(&a(&["prog", "--file"])).is_err()); // missing value
        assert!(parse_args(&a(&["prog", "--i", "bogus"])).is_err());
        assert!(parse_args(&a(&["prog", "--filter-freq-mad", "x"])).is_err());
    }

    #[test]
    fn test_parse_version() {
        let p = parse_args(&a(&["prog", "--version"])).unwrap();
        assert!(p.version_flag);
    }

    #[test]
    fn test_parse_rfi_flag() {
        let p = parse_args(&a(&["prog", "--file", "obs.hdf", "--rfi"])).unwrap();
        assert!(p.rfi_flag);
    }

    #[test]
    fn test_apply_mode_implications_rfi_conflicts() {
        let mut p = parse_args(&a(&["prog", "--rfi", "--i", "0", "--j", "1"])).unwrap();
        assert!(apply_mode_implications(&mut p).is_err());
        let mut p = parse_args(&a(&["prog", "--rfi", "--test-antennas"])).unwrap();
        assert!(apply_mode_implications(&mut p).is_err());
        let mut p = parse_args(&a(&["prog", "--rfi", "--benchmark"])).unwrap();
        assert!(apply_mode_implications(&mut p).is_err());
        let mut p = parse_args(&a(&["prog", "--rfi", "--file", "obs.hdf"])).unwrap();
        assert!(apply_mode_implications(&mut p).is_ok());
    }

    #[test]
    fn test_apply_mad_filters_threshold() {
        let mut results = vec![
            acquisition::GpsPrnResult { sv: "GPS01".into(), strengths: vec![], phases: vec![], freqs: vec![],
                cn0_acr: None, phase_median: Some(0.1), phase_mad: Some(0.001),
                freq_median: Some(0.0), freq_mad: Some(50.0) },
            acquisition::GpsPrnResult { sv: "GPS02".into(), strengths: vec![], phases: vec![], freqs: vec![],
                cn0_acr: None, phase_median: Some(0.2), phase_mad: Some(0.005),
                freq_median: Some(0.0), freq_mad: Some(20.0) },
            acquisition::GpsPrnResult { sv: "GPS03".into(), strengths: vec![], phases: vec![], freqs: vec![],
                cn0_acr: None, phase_median: None, phase_mad: None,
                freq_median: None, freq_mad: None },
        ];
        let removed = apply_mad_filters(&mut results, Some(0.002), None);
        assert_eq!(removed, 1);
        assert_eq!(results.len(), 2);
        assert!(results.iter().all(|r| r.sv != "GPS02"));
        assert!(results.iter().any(|r| r.sv == "GPS03"));
    }

    #[test]
    fn test_apply_mad_filters_freq() {
        let mut results = vec![
            acquisition::GpsPrnResult { sv: "GPS01".into(), strengths: vec![], phases: vec![], freqs: vec![],
                cn0_acr: None, phase_median: Some(0.0), phase_mad: Some(0.0),
                freq_median: Some(0.0), freq_mad: Some(100.0) },
            acquisition::GpsPrnResult { sv: "GPS02".into(), strengths: vec![], phases: vec![], freqs: vec![],
                cn0_acr: None, phase_median: Some(0.0), phase_mad: Some(0.0),
                freq_median: Some(0.0), freq_mad: Some(90.0) },
        ];
        let removed = apply_mad_filters(&mut results, None, Some(95.0));
        assert_eq!(removed, 1);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].sv, "GPS02");
    }

    #[test]
    fn test_apply_mad_filters_no_threshold_noop() {
        let mut results = vec![acquisition::GpsPrnResult { sv: "GPS01".into(), strengths: vec![], phases: vec![],
            freqs: vec![], cn0_acr: None, phase_median: Some(0.0), phase_mad: Some(1000.0),
            freq_median: Some(0.0), freq_mad: Some(1000.0) }];
        assert_eq!(apply_mad_filters(&mut results, None, None), 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].sv, "GPS01");
    }


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
        assert_eq!(p.sample_rate, Some(16.368e6));
        assert_eq!(p.center_freq, Some(4.092e6));
        assert_eq!(p.band, Some(2e6));
        assert_eq!(p.snr, Some(20.0));
        assert_eq!(p.output_file.as_deref(), Some("obs.hdf"));
    }
}
