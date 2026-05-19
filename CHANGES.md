# Changes

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
