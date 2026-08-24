# Changes

## UNRELEASED

### Changed

- **`--test-antennas` reference rule relaxed** — a satellite is a reference
  satellite if its ACR C/N0 is ≥ `--test-min-cn0` on **at least one** antenna
  (previously: on every antenna).  A broken or degraded antenna no longer
  disqualifies satellites other antennas see strongly; every antenna is still
  scored (its low C/N0 values contribute to its median) and per-antenna
  `n_sats` now reports 0 for such an antenna.  Default threshold raised
  40 → 44 dB-Hz.

## v0.3.2

### Added

- **`--test-antennas` mode** — estimates the *relative* C/N0 quality of each
  radio/antenna from one acquisition snapshot.  It acquires with ACR C/N0
  estimation (implies `--cn0`, and `--all` unless a constellation flag is
  given), finds satellites commonly visible with good signal strength
  (ACR C/N0 ≥ `--test-min-cn0`, default 40 dB-Hz, on **every** selected
  antenna), then reports each antenna's median C/N0 over that reference set
  with an offset relative to the best antenna (`offset_db`) and an overall
  rank.  New flags `--test-antennas` / `--test-min-cn0`; new module
  `src/antenna_test.rs`.
- **Simulator per-antenna gains** — `antenna_signals_with_gains` and
  `synthesize_with_gains` scale each antenna's signal component
  independently, enabling end-to-end tests of antenna-quality ranking.

### Fixed

- **`stats::median` used unstable `is_multiple_of`** — replaced with
  `% 2 == 0` so the crate builds on toolchains older than Rust 1.87.
- **MSRV declared and code made compatible** — `rust-version = "1.81"` in
  `Cargo.toml`, matching the lowest dependency (hdf5-reader/writer 0.9.1).
  Removed `Option::is_none_or` (stabilised in 1.82) from the MAD filters;
  the other unstable-API use (`is_multiple_of`) is also gone.

### Fixed (carried over from UNRELEASED)

- **Galileo / BeiDou / L1C `rustfft` buffer-too-small panic** — the single-PRN
  acquisition functions resampled the local code replica to a *floored whole
  number of code periods*, making the code shorter than the FFT size whenever
  the antenna signal length was not an exact multiple of the code period. This
  triggered `Provided FFT buffer was too small` in `rustfft`. The code replica
  is now resampled to exactly the signal length (including a partial final
  period), matching the GPS / SBAS / QZSS path.

## v0.3.1

### Changed

- **`--ant` accepts comma-separated lists** — `--ant` now takes a
  comma-separated list of antenna indices (e.g. `--ant 0,2,5`) instead of a
  single index.  The old single-index form (`--ant 0`) still works.  The
  `ant_filter` parameter in all `acquire_all_*` functions changed from
  `Option<usize>` to `Option<Vec<usize>>`.

## v0.3.0

### Added

- **`--version` flag** — prints the program name and version (`tart-gnss-acquire
  v0.3.0`) then exits.
- **Parallel frequency-bin correlation** — each single-PRN acquisition now
  processes its ~41 Doppler frequency bins in parallel via `rayon::par_iter()`
  instead of sequentially.  Combined with the existing PRN-level parallelism,
  this gives a substantial speedup on multi-core systems, especially for the
  larger FFT sizes (BeiDou B1C and GPS L1C at 320K samples).

### Changed

- **FFT plan hoisting** — in `acquisition.rs`, `galileo.rs`, `beidou.rs`, and
  `l1c.rs`, forward and inverse FFT plans are now created once per call rather
  than being regenerated inside each frequency-bin iteration (~41 redundant
  plan constructions eliminated per PRN per antenna).  `sbas.rs` and `qzss.rs`
  already had this optimisation and are unchanged.
- **Extracted `src/correlate.rs`** — the identical ~80-line frequency-bin
  correlation loop (IQ wipe-off, forward FFT, multiply by code conjugate,
  inverse FFT, magnitude, argmax, second-peak search) previously duplicated
  across all six constellation modules is now a single shared
  `correlate_code()` function.  Each module's `acquire_*_single` function
  shrinks from ~115 lines to ~55 lines.

## v0.2.1

### Added

- **`--benchmark` flag** — runs one PRN per constellation and reports
  timing broken into startup (antenna-data extraction) and per-PRN
  search, with an aggregate acquisitions/second score.  Useful for
  profiling acquisition throughput without running a full all-PRN
  sweep.

### Changed

- **CI workflow** — `actions/checkout` updated from v4 to v6;
  `rust-lang/crates-io-auth-action` pinned from v1 to v1.0.4.

## v0.2.0

### Added

- **`--prn` flag** — restricts acquisition to specific PRN/SV numbers
  (comma-separated list).  Works across all constellation modes: GPS L1 C/A,
  Galileo E1-C, BeiDou B1C, SBAS, QZSS, and GPS L1C.  Each `acquire_all_*`
  function now takes an optional `prn_filter` parameter that feeds
  `into_par_iter` directly, avoiding work for unlisted PRNs.

### Fixed

- **Crate package size** — the published crate exceeded the 10 MB crates.io
  limit because `data/` (12 MB HDF5 observation file), `junk/` (5.4 MB of
  PDFs and reference scripts), and an untracked `test.json` were included
  in the package.  Added `exclude = ["data/*", "junk/*", "test.json"]` to
  `Cargo.toml`.  Package size dropped from ~18 MB to ~270 KB.

## v0.1.7

### Added

- **QZSS L1 C/A acquisition** — new `--qzss` flag for QZSS (Quasi-Zenith
  Satellite System) L1 C/A all-PRN search.  Covers PRNs 184–206 (23 PRNs)
  with G2 delay values from the L1 C/A PRN Code Assignments (Jan 2026).
  Reuses the existing GPS C/A gold-code generator (FFT circular
  cross-correlation).  The `--all` flag now includes QZSS alongside GPS,
  Galileo, BeiDou, SBAS, and L1C.
- **`--cn0` flag** — documented in README.  Enables per-antenna ACR C/N0
  estimation (added in v0.1.6 but not previously listed in the README
  flag table).

### Changed

- **SBAS PRN range extended** — SBAS now searches PRNs 120–158 (was
  120–138).  The 20 additional PRNs include operational satellites from
  MSAS, SDCM, KASS, BDSBAS, Pak-SBAS, ASECNA, A-SBAS, EGNOS, and UK SBAS.
  G2 delay values sourced from the Jan 2026 PRN Code Assignments document.
- **SBAS first-10-chips test** — added cross-validation of generated C/A
  codes for all 39 SBAS PRNs against the reference octal values from the
  PRN assignment document.
- **Performance optimisations** — several changes reduce per-PRN acquisition
  time and memory churn, giving ~2× speedup for `--all` mode:
  - Legendre sequences (L1C and BeiDou, length 10223/10243) now cached via
    `LazyLock`, computed once instead of per-PRN.
  - BOC(1,1) pilot codes for L1C and BeiDou (20460 samples × 63 PRNs ≈ 10 MB
    each) pre-computed and cached on first use.
  - Per-frequency-bin `iq` and `corr` buffers pre-allocated once and reused
    across the ~21-bin search loop, eliminating repeated large allocations.
  - `FftPlanner` moved to a thread-local static, reused across all PRN calls
    within a thread instead of being recreated per call.
  - Antenna data pre-converted from `f64` to `f32` and phase-ramp
    pre-computed once in each `acquire_all_*` function, shared across all PRN
    iterations instead of repeated per PRN.
  - `sort_by_key` replaced with `sort_by` to avoid per-element `String` clones
    during result ordering.
- **README** — updated with all CLI flags (`--l1c`, `--qzss`, `--cn0`,
  `--output`, `--debug`), expanded flag table, and current PRN ranges.

### Fixed

- **README SBAS range** — corrected from 120–138 to 120–158.

## v0.1.6

### Added

- **ACR C/N0 estimation** — per-antenna Carrier-to-Noise ratio (C/N0)
  estimates based on the Acquisition Correlation Ratio method from
  Ma, Li, Zhou & Lu (2024), GPS Solutions 28:143,
  doi:10.1007/s10291-024-01666-y.  Each per-PRN result now
  includes an optional `cn0_acr` field: a vector of C/N0 values in
  dB-Hz, one per antenna.  The estimate uses the squared ratio of the
  global correlation peak to the largest off-peak (>1 chip) correlation
  at the same frequency, mapped through a lookup table derived from the
  third-order hypothesis model for 1023-chip Gold codes.
- **`--cn0` flag** — C/N0 estimation is opt-in; pass `--cn0` to enable
  the ACR computation and include `cn0_acr` in the JSON output.

### Changed

- All acquisition functions (`acquire_full`, `acquire_galileo_single`,
  `acquire_beidou_single`, `acquire_sbas`, `acquire_l1c_single`) now
  also track and return the second-highest correlation peak at the
  best frequency bin, used for the ACR C/N0 computation.

### Fixed

- Second-peak search now uses code-phase distance (modulo code period)
  instead of raw index distance, correctly excluding periodic repeats
  of the main peak across multi-epoch correlation arrays.

## v0.1.5

### Added

- **GPS L1C signal acquisition** — new `--l1c` flag for GPS L1C pilot (L1Cp)
  acquisition. Uses 10230-chip Weil codes derived from a Legendre sequence
  (prime p=10223) with a 7-bit expansion sequence insertion, per IS-GPS-800H
  Section 3.2.2.1.1.  Supports all 63 PRNs (1–63) with Weil index and
  insertion point parameters from Table 3.2-2, cross-validated against the
  initial/final 24-chip octal reference vectors.  The `--all` flag now
  includes L1C alongside GPS, Galileo, BeiDou, and SBAS.

## v0.1.4

### Fixed

- **BeiDou & Galileo FFT buffer panic** — `acquire_beidou_single` and
  `acquire_galileo_single` could panic with `Provided FFT buffer was too
  small` when fewer signal samples were available than one full code period.
  The FFT was planned for `total_samples` (≥ 1 code period) but the `iq`
  buffer built from the signal only had `min(total_samples, x.len())`
  elements.  Now `num_samples` is clamped to the available data before
  planning the FFT and generating the local code replica.
