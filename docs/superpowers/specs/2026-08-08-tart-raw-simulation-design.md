# Design: TART Raw-Data Simulation (`--simulate`)

**Date:** 2026-08-08
**Status:** Approved

## Goal

Add a `--simulate` mode to `tart-gnss-acquire` that synthesizes a raw TART radio
observation (0/1 unipolar packed samples) from a list of catalogue sources
(`{name, az, el, r, jy}`), with correct geometric phase delays per antenna, and
writes it as an HDF5 observation file directly loadable by
`tart-gnss-acquire --file`.

## Background

`tart-gnss-acquire` is a Rust port of the TART Python toolbox. Its observation
format (0/1 unipolar samples at fs 16.368 MHz, IF 4.092 MHz, packed into HDF5
`config`/`timestamp`/`data` datasets) matches TART's digitized IF signal. The
Python `tart.simulation` package produces such observations; this feature ports
the essential signal model (`antennas.py` + `radio.Max2769B`) in a **minimal
coherent form**: each source is a sinusoid at a distinct frequency offset near
the IF, phase-corrected by the geometric delay, rather than full band-limited
noise.

The catalogue (`github.com/tart-telescope/catalogue`) `/catalog` endpoint
returns sources with `name`, `el` (deg), `az` (deg), `r` (m), `jy` (Jy).

## Architecture

- New **`--simulate` top-level CLI mode** in `main.rs`, alongside the existing
  acquisition/correlation/benchmark modes.
- New **`src/simulator.rs`** module containing all simulation logic
  (types, phase delays, signal synthesis, packing) so it is unit-testable.
- **Dependencies:** add `hdf5-writer` (pure-Rust HDF5 encoder, no C deps); bump
  `hdf5-reader` to 0.9.x for family/version alignment (API-compatible with the
  existing `observation.rs`).

### Signal model (minimal coherent, model B)

Given sources and antenna positions, for antenna `i` at ENU position `pos_i`:

1. **Direction unit vector** `ŝ = [sin(az)·cos(el), cos(az)·cos(el), sin(el)]`
   (East, North, Up). Azimuth in degrees, azimuth 0 = North, increasing toward
   East.
2. **Geometric delay** (TART convention): source position `S = r·ŝ`;
   `Δ = (|S − pos_i| − r)/c`, `c = 2.99793e8` m/s. For sources with `r > 1e4`
   m, `r` is treated as `1e4` (plane-wave approximation), matching TART.
3. **Frequency assignment:** source `k` (index in source list) gets
   `f_k = fc + band·(k/(n−1) − 0.5)` for `n > 1`, else `fc` — evenly spread
   distinct offsets across `[fc − band/2, fc + band/2]`.
4. **Amplitude:** `A_k = K·√jy_k`, where `K` is a global gain constant.
5. **Sample value** for sample `t = n/fs`:
   `sig_i[n] = Σ_k A_k·cos(2π f_k·(t + Δ_i,k))`,
   plus Gaussian noise `N(0, σ)` per antenna where
   `σ = √( (Σ_k A_k²/2) / 10^(snr/10) )` (SNR in dB = 10·log10(signal_power /
   noise_power), signal_power = Σ A_k²/2).
6. **1-bit quantization:** `bit = 1 if sig_i[n] >= 0 else 0` (TART NRZ: sign,
   zero → +1, then `(sig+1)/2`). Per-antenna bits are packed MSB-first into
   bytes (matching `numpy.packbits` and `observation::unpack_bits`).

### Antenna positions

- **Provided:** JSON list of `[east, north, up]` triples (meters). `up` (z) is
  forced to 0.
- **Generated** (when no positions given): `--N` antennas uniformly at random in
  the EN plane within a circle of radius `diameter/2`, z = 0, using a fixed
  seed (deterministic; seeded via a fixed `u64` constant unless `--seed` given).

### HDF5 output

Write an observation file the existing `Observation::from_hdf5` can read:
- dataset `config`: JSON string `{"num_antenna":N,"sampling_frequency":fs}`.
- dataset `timestamp`: ISO-8601 UTC string (now).
- dataset `data`: 2-D `u8` array, shape `[N, row_bytes]`, packed 1-bit samples
  (MSB-first) per antenna.

### CLI surface

```
tart-gnss-acquire --simulate \
    --sources <catalogue.json>      # [ {name, az, el, r, jy}, ... ]
    --positions <antennas.json>     # OR omit and use --N
    --N <int>                       # antennas to generate randomly (no --positions)
    --diameter <m>                  # max array diameter for random layout (default 3.0)
    --seed <u64>                    # RNG seed (fixed default when omitted)
    --output <obs.hdf>              # HDF5 file to write
    --samples <int>                 # samples per antenna (default 2^16 = 65536)
    --sample-rate <Hz>              # default 16.368e6
    --center-freq <Hz>              # IF, default 4.092e6
    --band <Hz>                     # distinct-offset spread, default 2.0e6
    --snr <dB>                      # signal-to-noise, default user-supplied (tunable)
```

`--simulate` prints a brief per-source summary (name, az, el, r, jy, assigned
frequency) to stderr and writes the observation file.

## Defaults

| Parameter | Default |
|---|---|
| sample rate | 16.368 MHz |
| center freq (IF) | 4.092 MHz |
| band | 2.0 MHz (±1 MHz) |
| samples | 65 536 |
| array diameter | 3.0 m |
| seed | fixed constant (overridable) |
| gain K | 1.0 |
| snr | tunable (no silent default) |

## Testing

Unit tests in `src/simulator.rs`:
- **Phase-delay correctness:** two antennas + far source; cross-correlate the
  two generated tone trains and assert the measured phase offset equals the
  predicted `2π f Δ` between antennas.
- **Amplitude ∝ √jy:** doubling `jy` doubles the recovered tone amplitude ratio
  (√2).
- **Frequency spreading:** `n` sources map to evenly spaced distinct offsets.
- **Near-field vs plane-wave:** delay formula for a close source
  (`r < 1e4`) uses exact `|S − pos|`, and plane-wave behaviour for `r > 1e4`.
- **Random position generation:** deterministic under a fixed seed and bounded
  by `diameter/2`.
- **HDF5 round-trip:** a simulated file written by `hdf5-writer` is read back
  by `Observation::from_hdf5` and matches the expected dimensions and
  config/sampling rate.

## Out of Scope

- Full TART-fidelity bandwidth/noise model (Butterworth filtering, LO mixing).
- Making plain-tone sources detectable by the current C/A-code acquisition
  (sources are coherent tones, not GNSS codes). Phase-delay correctness is
  validated by the dedicated test above.
- Real catalogue API fetching (input is a local JSON file).
