mod acquisition;
mod config;
mod galileo;
mod observation;

use observation::Observation;
use serde::Serialize;

#[derive(Serialize)]
struct GpsAcquisitionOutput {
    prn: usize,
    /// Per-antenna signal strengths.
    strengths: Vec<f64>,
    /// Per-antenna code-phase offsets (fraction of a millisecond).
    phases: Vec<f64>,
    /// Per-antenna Doppler frequency offsets (Hz).
    freqs: Vec<f64>,
}

fn main() {
    let args: Vec<String> = std::env::args().collect();

    let mut file: Option<String> = None;
    let mut antenna_i: Option<usize> = None;
    let mut antenna_j: Option<usize> = None;
    let mut gps_prn: Option<usize> = None;
    let mut galileo_flag = false;
    let mut single_ant: Option<usize> = None;

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
                i += 1;
                gps_prn = Some(args[i].parse().expect("invalid integer for --gps"));
            }
            "--galileo" => {
                galileo_flag = true;
            }
            "--ant" => {
                i += 1;
                single_ant = Some(args[i].parse().expect("invalid integer for --ant"));
            }
            other => {
                eprintln!("unknown argument: {other}");
                std::process::exit(1);
            }
        }
        i += 1;
    }

    let file = file.unwrap_or_else(|| {
        eprintln!(
            "usage: {} --file <observation.hdf> [--i <i> --j <j>] [--gps <prn>] [--galileo]",
            args[0]
        );
        std::process::exit(1);
    });

    let obs = Observation::from_hdf5(&file).unwrap_or_else(|e| {
        eprintln!("error loading HDF5 file: {e}");
        std::process::exit(1);
    });

    let n_ant = obs.config.num_antenna();
    let sampling_freq = obs.get_sampling_rate();
    eprintln!("timestamp:   {}", obs.timestamp);
    eprintln!("antennas:    {n_ant}");
    eprintln!("sample rate: {sampling_freq} Hz");

    // --- GPS acquisition mode ----------------------------------------------
    if let Some(prn) = gps_prn {
        let num_samples_per_ms = (sampling_freq / 1000.0) as usize;
        // Use 2 ms of data per antenna as in check_sv_strength.py.
        let num_samples = 2 * num_samples_per_ms;

        let ant_range: Vec<usize> = if let Some(ant) = single_ant {
            if ant >= n_ant {
                eprintln!("antenna index {ant} out of range (have {n_ant})");
                std::process::exit(1);
            }
            vec![ant]
        } else {
            (0..n_ant).collect()
        };

        let mut strengths = Vec::with_capacity(ant_range.len());
        let mut phases = Vec::with_capacity(ant_range.len());
        let mut freqs = Vec::with_capacity(ant_range.len());

        for &ant_idx in &ant_range {
            let bipolar = obs.get_antenna(ant_idx);
            let mean = bipolar.iter().sum::<f64>() / bipolar.len() as f64;
            let raw: Vec<f64> = bipolar.iter().map(|&v| v - mean).collect();
            let chunk_len = num_samples.min(raw.len());

            let result = acquisition::acquire_full(
                &raw[..chunk_len],
                sampling_freq,
                4.092e6, // GPS L1 intermediate frequency
                6000.0,  // search bandwidth Hz
                prn,
            );

            eprintln!(
                "  ant {:2}: strength={:.3}  phase={:.6}  freq={:.1} Hz",
                ant_idx,
                result.signal_strength,
                result.codephase_frac,
                result.frequency
            );

            strengths.push(result.signal_strength);
            phases.push(result.codephase_frac);
            freqs.push(result.frequency);
        }

        let output = GpsAcquisitionOutput {
            prn,
            strengths,
            phases,
            freqs,
        };

        println!(
            "{}",
            serde_json::to_string_pretty(&output).expect("JSON serialization failed")
        );
        return;
    }

    // --- Galileo all-SV acquisition mode -----------------------------------
    if galileo_flag {
        // Galileo E1 carrier: 1575.42 MHz; typical IF: 4.092 MHz
        let galileo_if = 4.092e6;
        let galileo_search_band = 6000.0;

        if let Some(ant) = single_ant {
            if ant >= n_ant {
                eprintln!("antenna index {ant} out of range (have {n_ant})");
                std::process::exit(1);
            }
        }

        eprintln!("Running Galileo E1-C all-SV search (50 PRNs)...");
        let result =
            galileo::acquire_all_galileo(&obs, galileo_if, galileo_search_band, single_ant);

        println!(
            "{}",
            serde_json::to_string_pretty(&result).expect("JSON serialization failed")
        );
        return;
    }

    // --- Correlation mode --------------------------------------------------
    let i = antenna_i.unwrap_or_else(|| {
        eprintln!("missing --i <antenna_i> (or use --gps <prn> or --galileo)");
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
