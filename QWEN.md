# QWEN.md — gnss-tart

Rust CLI for GPS L1 C/A signal acquisition and antenna cross-correlation using
data from the TART radio telescope. Ported from the Python
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
# GPS acquisition for a single PRN
cargo run -- --file data/observation.hdf --gps 5

# Galileo E1-C all-SV pilot channel search
cargo run -- --file data/observation.hdf --galileo

# Antenna cross-correlation
cargo run -- --file data/observation.hdf --i 0 --j 1
```

Arguments:
- `--file <path>` — HDF5 observation file (required)
- `--gps <prn>`   — GPS PRN 1–38; runs acquisition mode
- `--galileo`     — run Galileo E1-C all-SV (PRN 1–50) pilot channel search
- `--ant <idx>`   — restrict GPS/Galileo acquisition to a single antenna
- `--i <idx>`     — first antenna index (correlation mode)
- `--j <idx>`     — second antenna index (correlation mode)

## Source layout

| File              | Purpose                                          |
|-------------------|--------------------------------------------------|
| `src/main.rs`     | CLI argument parsing, mode dispatch              |
| `src/config.rs`   | `Config` — deserialises TART JSON from HDF5      |
| `src/observation.rs` | `Observation` — HDF5 reader, unpacking, correlation |
| `src/acquisition.rs` | GPS C/A code gen + FFT circular cross-correlation |
| `src/galileo.rs`  | Galileo E1-C pilot acquisition (4092-chip codes) |
| `src/galileo_codes.rs` | Galileo E1-C primary codes (50 PRNs, hex-encoded) |

## Key dependencies

- `hdf5-reader` — reads HDF5 files (no HDF5 C lib required)
- `serde` / `serde_json` — config deserialisation, JSON output
- `chrono` — UTC timestamps
- `rustfft` / `num-complex` — FFT-based acquisition
- `rayon` — parallel PRN search for Galileo all-SV mode

## Repository

Part of the `elec-otago/molteno` monorepo (this crate lives at
`physics/gnss-tart/`).
