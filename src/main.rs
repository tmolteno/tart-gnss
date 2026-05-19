mod acquisition;
mod config;
mod galileo;
mod observation;

use observation::Observation;

fn main() {
    let args: Vec<String> = std::env::args().collect();

    let mut file: Option<String> = None;
    let mut antenna_i: Option<usize> = None;
    let mut antenna_j: Option<usize> = None;
    let mut gps_flag = false;
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
                gps_flag = true;
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
            "usage: {} --file <observation.hdf> [--i <i> --j <j>] [--gps] [--galileo] [--ant <idx>]",
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

    // --- GPS all-PRN acquisition mode --------------------------------------
    if gps_flag {
        if let Some(ant) = single_ant {
            if ant >= n_ant {
                eprintln!("antenna index {ant} out of range (have {n_ant})");
                std::process::exit(1);
            }
        }

        eprintln!(
            "Running GPS L1 C/A all-PRN search ({} PRNs)...",
            acquisition::GPS_NUM_SATS
        );
        let result = acquisition::acquire_all_gps(
            &obs,
            acquisition::GPS_IF,
            acquisition::GPS_SEARCH_BAND,
            single_ant,
        );

        println!(
            "{}",
            serde_json::to_string_pretty(&result).expect("JSON serialization failed")
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
        eprintln!("missing --i <antenna_i> (or use --gps or --galileo)");
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
