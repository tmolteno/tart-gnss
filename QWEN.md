# QWEN.md — tart-gnss-acquire

Rust CLI for GNSS signal acquisition (GPS L1 C/A, Galileo E1-C, BeiDou B1C)
and antenna cross-correlation using data from the TART radio telescope.
Ported from the Python `tart.tart.operation` module.

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

# All three constellations
cargo run -- --file data/observation.hdf --all

# Single antenna only
cargo run -- --file data/observation.hdf --gps --ant 0

# Antenna cross-correlation
cargo run -- --file data/observation.hdf --i 0 --j 1
```

Arguments:
- `--file <path>` — HDF5 observation file (required)
- `--gps`         — GPS L1 C/A all-PRN (1–38) search
- `--galileo`     — Galileo E1-C all-SV pilot (1–50) search
- `--beidou`      — BeiDou B1C all-SV pilot (1–63) search
- `--all`         — run all three constellations
- `--ant <idx>`   — restrict acquisition to a single antenna
- `--i <idx>`     — first antenna index (correlation mode)
- `--j <idx>`     — second antenna index (correlation mode)

## Source layout

| File                  | Purpose                                              |
|-----------------------|------------------------------------------------------|
| `src/main.rs`         | CLI argument parsing, mode dispatch, combined output  |
| `src/config.rs`       | `Config` — deserialises TART JSON from HDF5           |
| `src/observation.rs`  | `Observation` — HDF5 reader, unpacking, correlation   |
| `src/acquisition.rs`  | GPS C/A code gen + FFT circular cross-correlation     |
| `src/galileo.rs`      | Galileo E1-C pilot acquisition (4092-chip codes)      |
| `src/galileo_codes.rs`| Galileo E1-C primary codes (50 PRNs, hex-encoded)     |
| `src/beidou.rs`       | BeiDou B1C pilot acquisition (10230-chip Weil codes)  |
| `src/beidou_codes.rs` | BeiDou B1C Legendre/Weil code generation              |

## Key dependencies

- `hdf5-reader` — reads HDF5 files (no HDF5 C lib required)
- `serde` / `serde_json` — config deserialisation, JSON output
- `chrono` — UTC timestamps
- `rustfft` / `num-complex` — FFT-based acquisition
- `rayon` — parallel PRN search for multi-SV modes

## Repository

<https://github.com/tmolteno/tart-gnss>
