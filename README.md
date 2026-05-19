# tart-gnss-acquire

[![Crates.io](https://img.shields.io/crates/v/tart-gnss-acquire)](https://crates.io/crates/tart-gnss-acquire)

GNSS signal acquisition for the TART radio telescope — GPS L1 C/A, Galileo
E1-C, BeiDou B1C, and SBAS.

## Quick start

```bash
cargo install tart-gnss-acquire
tart-gnss-acquire --file observation.hdf --all
```

## Usage

```bash
# All four constellations
tart-gnss-acquire --file data/observation.hdf --all

# GPS only
tart-gnss-acquire --file data/observation.hdf --gps

# SBAS only
tart-gnss-acquire --file data/observation.hdf --sbas

# Single antenna
tart-gnss-acquire --file data/observation.hdf --galileo --ant 0

# Filter out results with high inter-antenna dispersion
tart-gnss-acquire --file data/observation.hdf --all \
    --filter-phase-mad 0.002 --filter-freq-mad 100.0

# Antenna cross-correlation
tart-gnss-acquire --file data/observation.hdf --i 0 --j 1
```

| Flag                      | Description                                         |
|---------------------------|-----------------------------------------------------|
| `--file PATH`             | HDF5 observation file (required)                    |
| `--gps`                   | GPS L1 C/A all-PRN search (1–38)                    |
| `--galileo`               | Galileo E1-C all-SV pilot search (1–50)             |
| `--beidou`                | BeiDou B1C all-SV pilot search (1–63)               |
| `--sbas`                  | SBAS L1 C/A all-PRN search (120–138)                |
| `--all`                   | Run GPS + Galileo + BeiDou + SBAS                   |
| `--ant IDX`               | Restrict acquisition to a single antenna            |
| `--filter-phase-mad X`    | Drop PRNs with phase MAD > X (multi-antenna only)   |
| `--filter-freq-mad X`     | Drop PRNs with frequency MAD > X (multi-antenna only) |
| `--i IDX`                 | First antenna index (correlation mode)              |
| `--j IDX`                 | Second antenna index (correlation mode)             |

## JSON output format

Acquisition results are written to stdout as JSON.  The top-level object
contains one key per requested constellation:

| Key        | Constellation | PRN range | SV label format |
|------------|---------------|-----------|-----------------|
| `gps`      | GPS L1 C/A    | 1–38      | `GPS01`         |
| `galileo`  | Galileo E1-C  | 1–50      | `GSAT01`        |
| `beidou`   | BeiDou B1C    | 1–63      | `BEIDOU01`      |
| `sbas`     | SBAS L1 C/A   | 120–138   | `SBAS120`       |

Each constellation value is an object with a `results` array.  Each element
in the array describes one satellite/PRN:

```json
{
  "gps": {
    "results": [
      {
        "sv": "GPS03",
        "strengths": [12.3, 11.8, 12.1],
        "phases": [0.001230, 0.001190, 0.001210],
        "freqs": [-1500.0, -1490.0, -1510.0],
        "phase_median": 0.00121,
        "phase_mad": 0.00002,
        "freq_median": -1500.0,
        "freq_mad": 10.0
      }
    ]
  }
}
```

### Per-PRN fields

| Field            | Type          | Description                                              |
|------------------|---------------|----------------------------------------------------------|
| `sv`             | string        | Constellation label (e.g. `GPS03`, `SBAS120`)           |
| `strengths`      | array[number] | Per-antenna peak correlation magnitudes                  |
| `phases`         | array[number] | Per-antenna code-phase offsets (fraction of code period) |
| `freqs`          | array[number] | Per-antenna Doppler frequency offsets (Hz)               |
| `phase_median`   | number?       | Median phase across antennas (only when >1 antenna)      |
| `phase_mad`      | number?       | Median absolute deviation of phase (multi-antenna only)  |
| `freq_median`    | number?       | Median frequency across antennas (multi-antenna only)    |
| `freq_mad`       | number?       | Median absolute deviation of frequency (multi-antenna only) |

The `_median` and `_mad` fields are omitted when acquisition is restricted
to a single antenna via `--ant`.  Use `--filter-phase-mad` /
`--filter-freq-mad` to suppress results with large inter-antenna spread.

## License

GPL-3.0

---

```
                ,odO
             .dPYqP9
            dY.  dP                  made with
           dY   dB
           Yb   Yb                  D E E P S E E K
           '9._.dP
            'VPY"                       v4
```

https://deepseek.ai
