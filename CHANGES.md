# Changes

## v0.1.4

### Fixed

- **BeiDou & Galileo FFT buffer panic** — `acquire_beidou_single` and
  `acquire_galileo_single` could panic with `Provided FFT buffer was too
  small` when fewer signal samples were available than one full code period.
  The FFT was planned for `total_samples` (≥ 1 code period) but the `iq`
  buffer built from the signal only had `min(total_samples, x.len())`
  elements.  Now `num_samples` is clamped to the available data before
  planning the FFT and generating the local code replica.
