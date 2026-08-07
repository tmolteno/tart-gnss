# Design: Unit Tests for tart-gnss-acquire

**Date:** 2026-08-08
**Status:** Approved

## Goal

Add comprehensive unit tests to `tart-gnss-acquire` for the currently untested
signal-processing and configuration paths. Tests only — no optimization or
behavior changes, with one minimal, behavior-preserving refactor of `main.rs`
to make its logic testable.

## Background / Coverage Gap

The crate already has 67 unit tests concentrated on the code generators
(C/A gold codes, Weil/Legendre codes, SBAS/QZSS delay tables, hex parsing)
plus `stats`, `acr`, and `observation`. The following are effectively
untested:

- `correlate_code` in `src/correlate.rs` — the shared FFT cross-correlation
  core used by every constellation (zero tests).
- `acquire_full` / `acquire_*_single` — the per-PRN acquisition search math.
- `src/config.rs` — `Config::from_json` parsing.
- `src/main.rs` — CLI arg parsing and MAD filtering (embedded in `main`).
- `src/observation.rs` — some edge cases (`get_means`, `get_antenna`,
  validation).

## Approach

### 1. Shared synthetic-signal test helper

A shared test helper that builds a deterministic synthetic observation: a known
code modulated onto a carrier at a specified code-phase delay and Doppler
offset, embedded in weak noise using a seeded PRNG. This lets acquisition tests
verify the pipeline *recovers* known inputs rather than merely not panicking.

### 2. `correlate.rs` — core correlation tests

Test `correlate_code` directly:
- A code-delayed signal produces a peak at the expected code phase
  (`best_codephase`) and the expected frequency bin.
- Noise-only input runs without panic and returns bounded values.
- Corner cases: single frequency bin, code-period wrap (delay near the end of
  the period), all-equal bins.

### 3. Acquisition functions

- `acquire_full` (GPS C/A): end-to-end — inject code at known delay + Doppler,
  assert recovered `codephase_frac`, `frequency`, and `signal_strength` ≈
  injection.
- `acquire_galileo_single`, `acquire_beidou_single`, `acquire_l1c_single`,
  `acquire_sbas`, `acquire_qzss`: targeted tests — valid `codephase ∈ [0,1)`,
  frequency within the search band, nonzero peak on a synthetic signal, plus a
  couple of actual recovery checks where cheap.
- `acquire_all_*`: `prn_filter` restricts results, `ant_filter` selects
  antennas, and output fields (`sv` naming, median/MAD present for >1 antenna).

### 4. `config.rs`

`from_json`: valid parse, alias field names (`num_antenna`,
`sampling_frequency`), missing-field error, malformed JSON error, accessor
correctness.

### 5. `main.rs` refactor + tests (small, behavior-preserving)

Extract into testable units:
- **Arg parsing** → `parse_args(&[String]) -> Result<ParsedArgs, String>`
  returning a struct of all flags; `main` calls it and handles errors. Tests
  cover each flag, list parsing, `--all` implying sub-flags, unknown-arg error,
  and `--version`.
- **MAD filtering** → a function applying phase/freq MAD thresholds to a
  results list. Tests cover below-threshold filtering, boundary behavior, and
  `None` MAD (single antenna) passing.

### 6. `observation.rs` additions

Tests for `get_means`, `get_antenna` bipolar mapping + out-of-range panic,
sample-count mismatch panic in `correlate`, and `random()` dimensions.

## Testing Approach

- Deterministic, seeded data — no wall-clock dependence.
- Acquisition-path tests use small FFT sizes / short signals to keep the suite
  fast.
- Verify with `cargo test`; run `cargo clippy` if available.

## Out of Scope

- Any optimization changes.
- Changes to acquisition/search behavior.
- Real-HDF5 integration tests (would require the data files).
