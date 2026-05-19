# tart-gnss-acquire

[![Crates.io](https://img.shields.io/crates/v/tart-gnss-acquire)](https://crates.io/crates/tart-gnss-acquire)

GNSS signal acquisition for the TART radio telescope — GPS L1 C/A, Galileo E1-C, and BeiDou B1C.

## Quick start

```bash
cargo install tart-gnss-acquire
tart-gnss-acquire --file observation.hdf --all
```

## Usage

```bash
# All three constellations
tart-gnss-acquire --file data/observation.hdf --all

# GPS only
tart-gnss-acquire --file data/observation.hdf --gps

# Single antenna
tart-gnss-acquire --file data/observation.hdf --galileo --ant 0

# Antenna cross-correlation
tart-gnss-acquire --file data/observation.hdf --i 0 --j 1
```

| Flag          | Description                                      |
|---------------|--------------------------------------------------|
| `--file PATH` | HDF5 observation file (required)                 |
| `--gps`       | GPS L1 C/A all-PRN search (1–38)                 |
| `--galileo`   | Galileo E1-C all-SV pilot search (1–50)          |
| `--beidou`    | BeiDou B1C all-SV pilot search (1–63)            |
| `--all`       | Run GPS + Galileo + BeiDou                       |
| `--ant IDX`   | Restrict acquisition to a single antenna         |
| `--i IDX`     | First antenna index (correlation mode)           |
| `--j IDX`     | Second antenna index (correlation mode)          |

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
