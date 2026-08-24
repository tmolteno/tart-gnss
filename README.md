# tart-gnss-acquire

[![Crates.io](https://img.shields.io/crates/v/tart-gnss-acquire)](https://crates.io/crates/tart-gnss-acquire)

GNSS signal acquisition for the TART radio telescope — GPS L1 C/A, GPS L1C,
Galileo E1-C, BeiDou B1C, SBAS, and QZSS.

## Quick start

```bash
cargo install tart-gnss-acquire
tart-gnss-acquire --file observation.hdf --all
```

## Usage

```bash
# All constellations
tart-gnss-acquire --file data/observation.hdf --all

# GPS L1 C/A only
tart-gnss-acquire --file data/observation.hdf --gps

# GPS L1C only
tart-gnss-acquire --file data/observation.hdf --l1c

# SBAS only
tart-gnss-acquire --file data/observation.hdf --sbas

# QZSS only
tart-gnss-acquire --file data/observation.hdf --qzss

# Single PRN
tart-gnss-acquire --file data/observation.hdf --gps --prn 3

# Multiple PRNs
tart-gnss-acquire --file data/observation.hdf --galileo --prn 1,5,12

# Single antenna
tart-gnss-acquire --file data/observation.hdf --galileo --ant 0

# With C/N0 estimation
tart-gnss-acquire --file data/observation.hdf --gps --cn0

# Save output to file
tart-gnss-acquire --file data/observation.hdf --all --output results.json

# Filter out results with high inter-antenna dispersion
tart-gnss-acquire --file data/observation.hdf --all \
    --filter-phase-mad 0.002 --filter-freq-mad 100.0

# Antenna cross-correlation
tart-gnss-acquire --file data/observation.hdf --i 0 --j 1

# Rank antennas by relative C/N0 quality
tart-gnss-acquire --file data/observation.hdf --test-antennas

# Same, but restrict to GPS and raise the reference-satellite threshold
tart-gnss-acquire --file data/observation.hdf --test-antennas --gps \
    --test-min-cn0 45

# Benchmark (single PRN per constellation, timed breakdown)
tart-gnss-acquire --file data/observation.hdf --benchmark

# Print version
tart-gnss-acquire --version
```

| Flag                      | Description                                         |
|---------------------------|-----------------------------------------------------|
| `--file PATH`             | HDF5 observation file (required)                    |
| `--gps`                   | GPS L1 C/A all-PRN search (1–38)                    |
| `--galileo`               | Galileo E1-C all-SV pilot search (1–50)             |
| `--beidou`                | BeiDou B1C all-SV pilot search (1–63)               |
| `--sbas`                  | SBAS L1 C/A all-PRN search (120–158)                |
| `--qzss`                  | QZSS L1 C/A all-PRN search (184–206)                |
| `--l1c`                   | GPS L1C all-SV pilot search (1–63)                  |
| `--all`                   | Run all constellations                              |
| `--prn A,B,...`           | Restrict search to specific PRN numbers             |
| `--ant IDX`               | Restrict acquisition to a single antenna            |
| `--cn0`                   | Enable ACR C/N0 estimation per antenna              |
| `--output PATH`           | Write JSON output to file instead of stdout         |
| `--debug`                 | Print per-PRN, per-antenna diagnostics to stderr    |
| `--filter-phase-mad X`    | Drop PRNs with phase MAD > X (multi-antenna only)   |
| `--filter-freq-mad X`     | Drop PRNs with frequency MAD > X (multi-antenna only) |
| `--i IDX`                 | First antenna index (correlation mode)              |
| `--j IDX`                 | Second antenna index (correlation mode)             |
| `--benchmark`             | Single-PRN timing with startup/search breakdown     |
| `--test-antennas`         | Rank antennas by relative C/N0 quality (implies `--cn0` and, unless a constellation flag is given, `--all`) |
| `--test-min-cn0 X`        | Reference-satellite C/N0 threshold in dB-Hz (default 44.0) |
| `--version`               | Print version and exit                              |

## Simulating data

Instead of reading a recorded observation, you can generate a synthetic HDF5
observation file with `--simulate`.  This is useful for testing acquisition
without a telescope recording.

```bash
# Generate an observation from a source catalogue and random antenna layout
tart-gnss-acquire --simulate --sources catalogue.json --N 16 --snr 20 \
    --output simulation.hdf

# Use a fixed antenna layout from a positions file
tart-gnss-acquire --simulate --sources catalogue.json --positions positions.json \
    --snr 20 --output simulation.hdf

# Then acquire from the generated file
tart-gnss-acquire --file simulation.hdf --all
```

The generated file is written to `--output <path>` (default `simulation.hdf`).
The source catalogue is a JSON array of sources, each with `name`, `az`
(azimuth, degrees), `el` (elevation, degrees), `r` (distance, m), and `jy`
(flux density, Jy).  Antenna positions come from either a `--positions` JSON
file or `--N <int>` random positions spread over a `--diameter` (default 3.0 m)
disk.

| Flag                | Description                                            |
|---------------------|--------------------------------------------------------|
| `--simulate`        | Generate a synthetic HDF5 observation instead of acquiring |
| `--sources PATH`    | Source catalogue JSON file (required)                  |
| `--positions PATH`  | Antenna positions JSON file (alternative to `--N`)     |
| `--N INT`           | Number of random antenna positions to generate         |
| `--diameter X`      | Disk diameter for random positions (m, default 3.0)    |
| `--snr X`           | Signal-to-noise ratio in dB (required)                 |
| `--seed INT`        | RNG seed for reproducible random positions/noise       |
| `--sample-rate X`   | Sample rate in Hz (default 16.368e6)                   |
| `--center-freq X`   | Centre frequency in Hz (default 4.092e6)               |
| `--band X`          | Signal bandwidth in Hz (default 2.0e6)                 |

## JSON output format

Acquisition results are written to stdout as JSON.  The top-level object
contains one key per requested constellation:

| Key        | Constellation | PRN range | SV label format |
|------------|---------------|-----------|-----------------|
| `gps`      | GPS L1 C/A    | 1–38      | `GPS01`         |
| `galileo`  | Galileo E1-C  | 1–50      | `GSAT01`        |
| `beidou`   | BeiDou B1C    | 1–63      | `BEIDOU01`      |
| `sbas`     | SBAS L1 C/A   | 120–158   | `SBAS120`       |
| `l1c`      | GPS L1C       | 1–63      | `GPSL1C01`      |
| `qzss`     | QZSS L1 C/A   | 184–206   | `QZSS184`       |

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
| `sv`             | string        | Constellation label (e.g. `GPS03`, `QZSS195`)           |
| `strengths`      | array[number] | Per-antenna peak correlation magnitudes                  |
| `phases`         | array[number] | Per-antenna code-phase offsets (fraction of code period) |
| `freqs`          | array[number] | Per-antenna Doppler frequency offsets (Hz)               |
| `cn0_acr`        | array[number]?| Per-antenna ACR C/N0 estimates in dB-Hz (only with `--cn0`) |
| `phase_median`   | number?       | Median phase across antennas (only when >1 antenna)      |
| `phase_mad`      | number?       | Median absolute deviation of phase (multi-antenna only)  |
| `freq_median`    | number?       | Median frequency across antennas (multi-antenna only)    |
| `freq_mad`       | number?       | Median absolute deviation of frequency (multi-antenna only) |

The `_median` and `_mad` fields are omitted when acquisition is restricted
to a single antenna via `--ant`.  Use `--filter-phase-mad` /
`--filter-freq-mad` to suppress results with large inter-antenna spread.

### Antenna test output (with `--test-antennas`)

With `--test-antennas`, the acquisition JSON is replaced by an
`antenna_test` report estimating the *relative* C/N0 quality of each
radio/antenna.  It finds satellites commonly visible with good signal
strength (ACR C/N0 ≥ `--test-min-cn0` on **at least one** antenna; a broken
antenna never reaches the threshold itself, but satellites it can still see
are scored for all antennas), then scores each antenna as the median C/N0
over that reference set — every antenna is scored, including ones below the
threshold:

```json
{
  "antenna_numbers": [0, 1, 2],
  "min_cn0_db_hz": 44.0,
  "n_reference_satellites": 4,
  "reference_satellites": ["GPS03", "GPS17", "GSAT21", "SBAS124"],
  "antennas": [
    {"antenna": 0, "n_sats": 4, "median_cn0_db_hz": 46.1, "offset_db": 0.0, "rank": 1},
    {"antenna": 1, "n_sats": 3, "median_cn0_db_hz": 44.8, "offset_db": -1.3, "rank": 2},
    {"antenna": 2, "n_sats": 0, "median_cn0_db_hz": 38.2, "offset_db": -7.9, "rank": 3}
  ],
  "per_satellite": [
    {"sv": "GPS03",   "cn0_db_hz": [46.0, 44.2, 37.5]},
    {"sv": "SBAS124", "cn0_db_hz": [46.5, 45.1, 39.0]}
  ]
}
```

| Field                     | Type            | Description                                             |
|---------------------------|-----------------|---------------------------------------------------------|
| `antenna_numbers`         | array[number]   | Antenna indices covered, in order                        |
| `min_cn0_db_hz`           | number          | Reference-satellite C/N0 threshold (default 44.0)       |
| `n_reference_satellites`  | number          | Reference satellites meeting the threshold on at least one antenna |
| `reference_satellites`    | array[string]   | SV labels of the reference set                         |
| `antennas[].antenna`      | number          | Antenna index                                          |
| `antennas[].n_sats`       | number          | Reference satellites this antenna met the threshold on (0 for a broken antenna) |
| `antennas[].median_cn0_db_hz` | number     | Median C/N0 over the reference set (dB-Hz)             |
| `antennas[].offset_db`    | number          | Median minus best antenna's median (dB; best = 0, negative = worse) |
| `antennas[].rank`         | number          | Rank by median C/N0 (1 = best)                         |
| `per_satellite[].sv`      | string          | Satellite label                                       |
| `per_satellite[].cn0_db_hz` | array[number] | Per-antenna C/N0 aligned with `antenna_numbers`        |

If no satellite meets the threshold on at least one antenna,
`n_reference_satellites` is 0, `reference_satellites`/`per_satellite`/
`antennas` are empty, and a warning is printed to stderr.  A short per-antenna
summary is also printed to stderr while the JSON report goes to stdout
(or `--output <path>`).

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
