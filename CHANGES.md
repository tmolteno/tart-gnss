# Changes

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
