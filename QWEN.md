# QWEN.md — tart-gnss-acquire

Rust CLI for GNSS signal acquisition (GPS L1 C/A, GPS L1C, Galileo E1-C,
BeiDou B1C, SBAS L1 C/A, QZSS L1 C/A) and antenna cross-correlation using
data from the TART radio telescope.  Ported from the Python
`tart.tart.operation` module.

## Build & test

```bash
cargo build              # debug build
cargo build --release    # optimised build
cargo test               # run all unit tests
cargo clippy             # lint (when available)
```

## Usage

```bash
# GPS L1 C/A all-PRN search
cargo run -- --file data/observation.hdf --gps

# Galileo E1-C all-SV pilot channel search
cargo run -- --file data/observation.hdf --galileo

# BeiDou B1C all-SV pilot channel search
cargo run -- --file data/observation.hdf --beidou

# SBAS L1 C/A all-PRN search
cargo run -- --file data/observation.hdf --sbas

# QZSS L1 C/A all-PRN search
cargo run -- --file data/observation.hdf --qzss

# GPS L1C all-SV pilot search
cargo run -- --file data/observation.hdf --l1c

# All constellations
cargo run -- --file data/observation.hdf --all

# Antenna relative C/N0 quality (implies --cn0 and --all unless a
# constellation flag is given)
cargo run -- --file data/observation.hdf --test-antennas

# Same, but only GPS, with a higher reference threshold
cargo run -- --file data/observation.hdf --test-antennas --gps --test-min-cn0 45

# Single PRN
cargo run -- --file data/observation.hdf --gps --prn 3

# Multiple PRNs
cargo run -- --file data/observation.hdf --galileo --prn 1,5,12

# Single antenna only
cargo run -- --file data/observation.hdf --gps --ant 0

# Filter by inter-antenna MAD
cargo run -- --file data/observation.hdf --all --filter-phase-mad 0.002 --filter-freq-mad 100.0

# Antenna cross-correlation
cargo run -- --file data/observation.hdf --i 0 --j 1
```

Arguments:
- `--file <path>`           — HDF5 observation file (required)
- `--gps`                   — GPS L1 C/A all-PRN (1–38) search
- `--galileo`               — Galileo E1-C all-SV pilot (1–50) search
- `--beidou`                — BeiDou B1C all-SV pilot (1–63) search
- `--sbas`                  — SBAS L1 C/A all-PRN (120–158) search
- `--qzss`                  — QZSS L1 C/A all-PRN (184–206) search
- `--l1c`                   — GPS L1C all-SV pilot (1–63) search
- `--all`                   — run all constellations
- `--prn <a,b,...>`         — restrict search to specific PRN numbers (comma-separated)
- `--ant <idx>`             — restrict acquisition to a single antenna
- `--cn0`                   — enable ACR C/N0 estimation per antenna
- `--output <path>`         — write JSON output to file instead of stdout
- `--debug`                 — print per-PRN diagnostics to stderr
- `--filter-phase-mad <x>`  — drop PRNs with phase MAD > x (multi-antenna only)
- `--filter-freq-mad <x>`   — drop PRNs with frequency MAD > x (multi-antenna only)
- `--i <idx>`               — first antenna index (correlation mode)
- `--j <idx>`               — second antenna index (correlation mode)
- `--test-antennas`         — rank antennas by relative C/N0 (implies `--cn0`, and `--all` unless a constellation flag is given)
- `--test-min-cn0 <x>`      — reference-satellite C/N0 threshold in dB-Hz (default 44.0; satellite qualifies if ≥ this on **at least one** antenna)

## Source layout

| File                  | Purpose                                              |
|-----------------------|------------------------------------------------------|
| `src/main.rs`         | CLI argument parsing, mode dispatch, combined output  |
| `src/antenna_test.rs` | antenna-relative C/N0 quality report (`--test-antennas`) |
| `src/config.rs`       | `Config` — deserialises TART JSON from HDF5           |
| `src/observation.rs`  | `Observation` — HDF5 reader, unpacking, correlation   |
| `src/acquisition.rs`  | GPS C/A code gen + FFT circular cross-correlation     |
| `src/galileo.rs`      | Galileo E1-C pilot acquisition (4092-chip codes)      |
| `src/galileo_codes.rs`| Galileo E1-C primary codes (50 PRNs, hex-encoded)     |
| `src/beidou.rs`       | BeiDou B1C pilot acquisition (10230-chip Weil codes)  |
| `src/beidou_codes.rs` | BeiDou B1C Legendre/Weil code generation              |
| `src/l1c.rs`          | GPS L1C pilot acquisition (10230-chip Weil codes)      |
| `src/l1c_codes.rs`    | GPS L1C Legendre/Weil code generation                  |
| `src/sbas.rs`         | SBAS L1 C/A acquisition (PRN 120–158, gold codes)      |
| `src/qzss.rs`         | QZSS L1 C/A acquisition (PRN 184–206, gold codes)      |
| `src/stats.rs`        | `median()` and `mad()` robust statistics               |

## Key dependencies

- `hdf5-reader` — reads HDF5 files (no HDF5 C lib required)
- `serde` / `serde_json` — config deserialisation, JSON output
- `chrono` — UTC timestamps
- `rustfft` / `num-complex` — FFT-based acquisition
- `rayon` — parallel PRN search for multi-SV modes

## Repository

<https://github.com/tmolteno/tart-gnss>
