# Unit Tests for tart-gnss-acquire — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add comprehensive unit tests for the untested signal-processing and configuration paths of `tart-gnss-acquire`, with one small behavior-preserving refactor of `main.rs` so its logic becomes testable.

**Architecture:** This is a binary-only crate (no `lib.rs`), so all tests are `#[cfg(test)]` unit tests inside each module. A new `#[cfg(test)]`-gated shared helper module (`src/testutil.rs`) provides synthetic-signal builders used by the per-module acquisition tests. `main.rs` gets its arg parsing and MAD filtering extracted into testable functions/types.

**Tech Stack:** Rust (edition 2021), `rustfft`, `num-complex`, `rayon`, `serde`/`serde_json`, `chrono`, `hdf5-reader`.

## Global Constraints

- Tests only — **no optimization or behavior changes** to acquisition/search logic.
- One permitted refactor: `main.rs` arg parsing + MAD filtering extraction, preserving CLI behavior.
- All tests must be deterministic (seeded PRNG, no wall-clock dependence).
- Acquisition-path tests use small FFT sizes / short signals to keep the suite fast.
- Verify with `cargo test` (all modules) before each commit. Run `cargo clippy --all-targets` if available.
- The crate name (for `env!` / module paths) is `tart-gnss-acquire`; unit tests use `crate::...` paths, never the external crate name.
- **FFT correlation convention (verified empirically):** the circular cross-correlation peak for a signal `s[n] = code[(n+delay) % period]` appears at recovered sample `(period - delay) % period`. Frequency recovered for a real carrier at `f0` has a sign ambiguity; `recovered.abs() ≈ |f0|` within half the frequency-bin spacing (~150 Hz), use a 250 Hz tolerance. These facts are baked into the test helper `expected_recovered_sample` and the test tolerances below.

---

### Task 1: Shared test utilities (`src/testutil.rs`)

**Files:**
- Create: `src/testutil.rs`
- Modify: `src/main.rs` (register the module, `#[cfg(test)]` only)

**Interfaces:**
- Produces (used by all later tasks):
  - `pub fn phasepoints(fs: f64, n: usize) -> Vec<f64>`
  - `pub fn synth_signal(code: &[f64], code_period: usize, delay: usize, f0: f64, fs: f64, noise_amp: f64, seed: u64) -> Vec<f32>`
  - `pub fn expected_recovered_sample(period: usize, delay: usize) -> usize  // (period - delay) % period`

- [ ] **Step 1: Create `src/testutil.rs`**

```rust
// Copyright (c) 2026 Tim Molteno <tim@elec.ac.nz>
// SPDX-License-Identifier: GPL-3.0

//! Shared helpers for unit tests. Compiled only under `cargo test`
//! (registered as `#[cfg(test)] mod testutil;` in main.rs).

use std::f64::consts::TAU;

/// Phase-accumulator points used by the acquisition pipeline: `TAU * i / fs`.
pub fn phasepoints(fs: f64, n: usize) -> Vec<f64> {
    let pc = TAU / fs;
    (0..n).map(|i| pc * i as f64).collect()
}

/// Build a synthetic signal from a known `code` (length == `num_samples`),
/// delayed circularly by `delay` samples relative to a `code_period`, modulated
/// onto a real carrier at `f0` Hz, plus weak deterministic noise.
///
/// Returns the de-meaned f32 signal the acquisition functions expect
/// (de-meaned because the code has zero mean).
pub fn synth_signal(
    code: &[f64],
    code_period: usize,
    delay: usize,
    f0: f64,
    fs: f64,
    noise_amp: f64,
    seed: u64,
) -> Vec<f32> {
    let n = code.len();
    let mut rng = XorShift(seed);
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let cp = (i + delay) % code_period;
        let carrier = (TAU * f0 * i as f64 / fs).cos();
        let noise = noise_amp * (2.0 * rng.next_f64() - 1.0);
        out.push((code[cp] * carrier + noise) as f32);
    }
    out
}

/// Expected recovered code-phase sample (mod one period) given the FFT
/// correlation convention `signal[i] = code[(i+delay) % period]` →
/// `best_codephase % period == (period - delay) % period`.
pub fn expected_recovered_sample(period: usize, delay: usize) -> usize {
    (period - delay) % period
}

/// Deterministic xorshift64 PRNG (independent of the crate's other PRNGs).
struct XorShift(u64);
impl XorShift {
    fn next_f64(&mut self) -> f64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        (x >> 11) as f64 / (1u64 << 53) as f64
    }
}
```

- [ ] **Step 2: Register the module in `src/main.rs`**

Add next to the other `mod` declarations:

```rust
#[cfg(test)]
mod testutil;
```

- [ ] **Step 3: Run build/tests to verify it compiles**

Run: `cargo test --release`
Expected: builds cleanly; existing 67 tests still pass.

- [ ] **Step 4: Commit**

```bash
git add src/testutil.rs src/main.rs
git commit -m "test: add shared synthetic-signal test utilities"
```

---

### Task 2: `correlate_code` unit tests

**Files:**
- Modify: `src/correlate.rs` (add tests at end of file)

**Interfaces:**
- Consumes: `correlate_code(signal: &[f32], phasepoints: &[f64], code_complex: &[Complex<f32>], fc: &[f64], num_samples, code_period_samples, sampling_freq, fft: Arc<dyn Fft<f32>>, ifft: Arc<dyn Fft<f32>>) -> CorrelationPeak` (fields: `best_peak: f32`, `second_peak: f32`, `best_freq_idx: usize`, `best_codephase: usize`); `crate::testutil::{phasepoints, synth_signal, expected_recovered_sample}`; `crate::acquisition::gold_code`.
- Produces: nothing for later tasks.

- [ ] **Step 1: Add the test module to `src/correlate.rs`**

Append:

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::acquisition::gold_code;
    use crate::testutil::{expected_recovered_sample, phasepoints, synth_signal};

    /// Build the pre-FFT'd, conjugated code replica exactly as the
    /// acquisition modules do.
    fn code_fft_conj(fft: Arc<dyn Fft<f32>>, code: &[f64]) -> Vec<Complex<f32>> {
        let mut cc: Vec<Complex<f32>> =
            code.iter().map(|&v| Complex::new(v as f32, 0.0)).collect();
        fft.process(&mut cc);
        for c in &mut cc {
            *c = c.conj();
        }
        cc
    }

    #[test]
    fn test_correlate_code_recovers_delay() {
        let period = 256usize;
        let code = gold_code(period as f64, 1, 1.0); // 256-sample GPS C/A code
        let fs = 1.023e6;

        let mut planner = FftPlanner::<f32>::new();
        let fft = planner.plan_fft_forward(period);
        let ifft = planner.plan_fft_inverse(period);

        let cc = code_fft_conj(fft.clone(), &code);

        for delay in [0usize, 50, 128, 200] {
            let sig = synth_signal(&code, period, delay, 0.0, fs, 0.05, 42);
            // phasepoints + fc=[0.0] means carrier wipe is identity.
            let peak = correlate_code(
                &sig, &phasepoints(fs, period), &cc, &[0.0],
                period, period, fs, fft.clone(), ifft.clone(),
            );
            let recovered = peak.best_codephase % period;
            assert_eq!(
                recovered, expected_recovered_sample(period, delay),
                "delay {delay} not recovered correctly"
            );
            assert_eq!(peak.best_freq_idx, 0);
            assert!(peak.best_peak > 0.0);
        }
    }

    #[test]
    fn test_correlate_code_noise_only_bounded() {
        let period = 256usize;
        let code = gold_code(period as f64, 1, 1.0);
        let fs = 1.023e6;
        let mut planner = FftPlanner::<f32>::new();
        let fft = planner.plan_fft_forward(period);
        let ifft = planner.plan_fft_inverse(period);
        let cc = code_fft_conj(fft.clone(), &code);

        let sig = synth_signal(&code, period, 0, 0.0, fs, 0.0, 7);
        let peak = correlate_code(
            &sig, &phasepoints(fs, period), &cc, &[0.0],
            period, period, fs, fft, ifft,
        );
        // Real code with zero noise reproduces the code itself: peak is large,
        // and second_peak is finite/positive.
        assert!(peak.best_peak > 0.0);
        assert!(peak.best_peak.is_finite());
        assert!(peak.second_peak > 0.0);
        assert!(peak.best_codephase < period);
    }
}
```

Note: `CorrelationPeak` and `correlate_code` are `pub(crate)`, so tests inside the same module access them directly (`use super::*`).

- [ ] **Step 2: Run the new tests**

Run: `cargo test --release correlate`
Expected: both tests PASS (delay recovery matches `expected_recovered_sample`).

- [ ] **Step 3: Commit**

```bash
git add src/correlate.rs
git commit -m "test: add correlate_code FFT cross-correlation tests"
```

---

### Task 3: GPS acquisition tests (`src/acquisition.rs`)

**Files:**
- Modify: `src/acquisition.rs` (add tests to the existing `#[cfg(test)] mod tests`)
- Modify: `src/observation.rs` / `src/config.rs` — only if needed via `Observation::new` / `Config::from_json` (both already public; no change required).

**Interfaces:**
- Consumes: `acquire_full(signal_f32, phasepoints, sampling_freq, center_freq, search_band, prn, samples_per_chunk) -> AcquisitionResult { prn, signal_strength, second_peak, codephase_frac, frequency }`; `acquire_all_gps(obs, center_freq, search_band, ant_filter: Option<Vec<usize>>, prn_filter: Option<&[usize]>, debug, cn0) -> GpsAllAcquisitionOutput { antenna_numbers, results: Vec<GpsPrnResult> }`; `gold_code(samples_per_code: f64, prn, epochs) -> Vec<f64>`; `crate::testutil::*`; `Observation::new`; `Config::from_json`.
- Produces: nothing for later tasks.

- [ ] **Step 1: Add synthetic-signal recovery and noise-contrast tests**

Append to the existing `use super::*;` test module in `src/acquisition.rs`:

```rust
use crate::testutil::{expected_recovered_sample, phasepoints, synth_signal};

#[test]
fn test_acquire_full_recovers_delay_and_doppler() {
    let fs = 2.046e6; // 2046 samples/ms
    let period = (fs / 1000.0) as usize;
    let code = gold_code(period as f64, 1, 2.0); // 2 epochs

    let delay = 500usize;
    let f0 = 1200.0f64;
    let sig = synth_signal(&code, period, delay, f0, fs, 0.05, 1);
    let r = acquire_full(
        &sig, &phasepoints(fs, sig.len()), fs, 0.0, 3000.0, 1, period,
    );

    let recovered_sample =
        (r.codephase_frac * period as f64).round() as usize;
    assert_eq!(
        recovered_sample, expected_recovered_sample(period, delay),
        "GPS code-phase delay not recovered"
    );
    // Real carrier → sign ambiguity; check magnitude within half a bin (~150 Hz).
    assert!(
        (r.frequency.abs() - f0.abs()).abs() <= 250.0,
        "doppler {} not near |f0|={}", r.frequency, f0
    );
    assert!(r.codephase_frac >= 0.0 && r.codephase_frac < 1.0);
}

#[test]
fn test_acquire_full_signal_strength_contrast() {
    let fs = 2.046e6;
    let period = (fs / 1000.0) as usize;
    let code = gold_code(period as f64, 1, 2.0);
    let n = code.len();

    let injected = synth_signal(&code, period, 123, 0.0, fs, 0.05, 2);
    let r_inj = acquire_full(&injected, &phasepoints(fs, n), fs, 0.0, 3000.0, 1, period);

    let noise_only = synth_signal(&code, period, 0, 0.0, fs, 0.05, 3);
    // Without a carrier/delay coherence the noise signal must not be modulated
    // by the code; use pure random noise instead for the baseline.
    let pure_noise: Vec<f32> = (0..n)
        .map(|i| ((i as f64 * 3.3).sin() * 0.05) as f32)
        .collect();
    let r_noise = acquire_full(&pure_noise, &phasepoints(fs, n), fs, 0.0, 3000.0, 1, period);

    assert!(
        r_inj.signal_strength > 10.0 * r_noise.signal_strength,
        "injected {0} should dominate noise {1}",
        r_inj.signal_strength, r_noise.signal_strength
    );
}
```

- [ ] **Step 2: Add `acquire_all_gps` filter tests**

Same test module:

```rust
#[test]
fn test_acquire_all_gps_prn_and_antenna_filters() {
    use crate::config::Config;
    use crate::observation::Observation;

    let fs = 2.046e6;
    let samples = 2 * (fs / 1000.0) as usize; // 2 ms
    // Two antennas, random 0/1 data of sufficient length.
    let data: Vec<Vec<u8>> = (0..2)
        .map(|a| {
            (0..samples).map(|i| ((i + a * 31) % 2) as u8).collect()
        })
        .collect();
    let cfg = Config::from_json(&format!(
        r#"{{"num_antenna":2,"sampling_frequency":{fs}}}"#
    ))
    .unwrap();
    let obs = Observation::new(chrono::Utc::now(), cfg, data);

    let out = acquire_all_gps(
        &obs, GPS_IF, GPS_SEARCH_BAND,
        Some(vec![0]), // ant filter
        Some(&[2, 5, 9]),
        false, false,
    );
    assert_eq!(out.antenna_numbers, vec![0]);
    assert_eq!(out.results.len(), 3);
    assert_eq!(out.results[0].sv, "GPS02");
    assert_eq!(out.results[1].sv, "GPS05");
    assert_eq!(out.results[2].sv, "GPS09");
    assert_eq!(out.results[0].strengths.len(), 1); // one antenna
    // Single antenna → no median/MAD
    assert!(out.results[0].phase_median.is_none());
    assert!(out.results[0].freq_mad.is_none());
}
```

- [ ] **Step 3: Run the new tests**

Run: `cargo test --release acquisition`
Expected: all acquisition tests PASS.

- [ ] **Step 4: Commit**

```bash
git add src/acquisition.rs
git commit -m "test: add GPS acquisition recovery and filter tests"
```

---

### Task 4: Galileo acquisition tests (`src/galileo.rs`)

**Files:**
- Modify: `src/galileo.rs` (extend existing `#[cfg(test)] mod tests`)

**Interfaces:**
- Consumes: `acquire_galileo_single(signal_f32, phasepoints, sampling_freq, center_freq, search_band, prn, samples_per_code_period) -> GalileoAcquisitionResult`; `e1c_code_resampled(samples_per_code: f64, prn, epochs) -> Vec<f64>`; `acquire_all_galileo(obs, center_freq, search_band, ant_filter, prn_filter, debug, cn0) -> GalileoAllAcquisitionOutput`; `GALILEO_E1_CHIPS` (= 4092); `GALILEO_E1_CODE_PERIOD` (= 0.004 s); `GALILEO_E1_NUM_SATS`; `crate::testutil::*`; `Observation::new`; `Config::from_json`.
- Produces: nothing for later tasks.

- [ ] **Step 1: Add Galileo recovery + noise-contrast tests**

Append to the existing test module in `src/galileo.rs`:

```rust
use crate::testutil::{expected_recovered_sample, phasepoints, synth_signal};

#[test]
fn test_acquire_galileo_single_recovers_delay() {
    let fs = 1.023e6; // one sample per chip; 4092 chips per 4 ms period
    let period = GALILEO_E1_CHIPS;
    let code = e1c_code_resampled(period as f64, 1, 2.0); // 2 epochs

    let delay = 900usize;
    let sig = synth_signal(&code, period, delay, 0.0, fs, 0.05, 11);
    let r = acquire_galileo_single(
        &sig, &phasepoints(fs, sig.len()), fs, 0.0, 3000.0, 1, period,
    );
    let recovered = (r.codephase_frac * period as f64).round() as usize;
    assert_eq!(recovered, expected_recovered_sample(period, delay));
    assert!(r.codephase_frac >= 0.0 && r.codephase_frac < 1.0);
    assert!(r.signal_strength > 0.0);
}

#[test]
fn test_acquire_galileo_single_noise_contrast() {
    let fs = 1.023e6;
    let period = GALILEO_E1_CHIPS;
    let code = e1c_code_resampled(period as f64, 1, 2.0);
    let n = code.len();

    let injected = synth_signal(&code, period, 600, 0.0, fs, 0.05, 12);
    let r_inj = acquire_galileo_single(
        &injected, &phasepoints(fs, n), fs, 0.0, 3000.0, 1, period,
    );
    let pure_noise: Vec<f32> = (0..n).map(|i| ((i as f64 * 2.7).sin() * 0.05) as f32).collect();
    let r_noise = acquire_galileo_single(
        &pure_noise, &phasepoints(fs, n), fs, 0.0, 3000.0, 1, period,
    );
    assert!(r_inj.signal_strength > 10.0 * r_noise.signal_strength);
}
```

- [ ] **Step 2: Run the new tests**

Run: `cargo test --release galileo`
Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git add src/galileo.rs
git commit -m "test: add Galileo E1-C acquisition tests"
```

---

### Task 5: BeiDou + L1C acquisition tests

**Files:**
- Modify: `src/beidou.rs` (extend existing `#[cfg(test)] mod tests`)
- Modify: `src/l1c.rs` (extend existing `#[cfg(test)] mod tests`)

**Interfaces:**
- Consumes:
  - `acquire_beidou_single(signal_f32, phasepoints, sampling_freq, center_freq, search_band, prn, samples_per_code_period) -> BeiDouAcquisitionResult`; `b1c_code_resampled(samples_per_code, prn, epochs) -> Vec<f64>`; `BEIDOU_B1C_CHIPS` (= 10230).
  - `acquire_l1c_single(signal_f32, phasepoints, sampling_freq, center_freq, search_band, prn, samples_per_code_period) -> L1CAcquisitionResult`; `l1c_code_resampled(samples_per_code, prn, epochs) -> Vec<f64>`; `L1C_CHIPS` (= 10230).
  - `crate::testutil::*`.
- Produces: nothing for later tasks.

- [ ] **Step 1: Add BeiDou recovery + noise-contrast tests to `src/beidou.rs`**

Append to the existing test module:

```rust
use crate::testutil::{expected_recovered_sample, phasepoints, synth_signal};

#[test]
fn test_acquire_beidou_single_recovers_delay() {
    let fs = 1.023e6; // one sample per chip; 10230 chips per 10 ms period
    let period = BEIDOU_B1C_CHIPS;
    let code = b1c_code_resampled(period as f64, 1, 2.0);

    let delay = 3000usize;
    let sig = synth_signal(&code, period, delay, 0.0, fs, 0.05, 21);
    let r = acquire_beidou_single(
        &sig, &phasepoints(fs, sig.len()), fs, 0.0, 3000.0, 1, period,
    );
    let recovered = (r.codephase_frac * period as f64).round() as usize;
    assert_eq!(recovered, expected_recovered_sample(period, delay));
    assert!(r.codephase_frac >= 0.0 && r.codephase_frac < 1.0);
}

#[test]
fn test_acquire_beidou_single_noise_contrast() {
    let fs = 1.023e6;
    let period = BEIDOU_B1C_CHIPS;
    let code = b1c_code_resampled(period as f64, 1, 2.0);
    let n = code.len();
    let injected = synth_signal(&code, period, 1000, 0.0, fs, 0.05, 22);
    let r_inj = acquire_beidou_single(
        &injected, &phasepoints(fs, n), fs, 0.0, 3000.0, 1, period,
    );
    let pure_noise: Vec<f32> = (0..n).map(|i| ((i as f64 * 1.9).sin() * 0.05) as f32).collect();
    let r_noise = acquire_beidou_single(
        &pure_noise, &phasepoints(fs, n), fs, 0.0, 3000.0, 1, period,
    );
    assert!(r_inj.signal_strength > 10.0 * r_noise.signal_strength);
}
```

- [ ] **Step 2: Add L1C tests to `src/l1c.rs`**

Append to the existing `#[cfg(test)] mod tests`:

```rust
use crate::testutil::{expected_recovered_sample, phasepoints, synth_signal};

#[test]
fn test_acquire_l1c_single_recovers_delay() {
    let fs = 1.023e6; // one sample per chip; 10230 chips per 10 ms period
    let period = L1C_CHIPS;
    let code = l1c_code_resampled(period as f64, 1, 2.0);

    let delay = 2500usize;
    let sig = synth_signal(&code, period, delay, 0.0, fs, 0.05, 31);
    let r = acquire_l1c_single(
        &sig, &phasepoints(fs, sig.len()), fs, 0.0, 3000.0, 1, period,
    );
    let recovered = (r.codephase_frac * period as f64).round() as usize;
    assert_eq!(recovered, expected_recovered_sample(period, delay));
    assert!(r.codephase_frac >= 0.0 && r.codephase_frac < 1.0);
}

#[test]
fn test_acquire_l1c_single_noise_contrast() {
    let fs = 1.023e6;
    let period = L1C_CHIPS;
    let code = l1c_code_resampled(period as f64, 1, 2.0);
    let n = code.len();
    let injected = synth_signal(&code, period, 700, 0.0, fs, 0.05, 32);
    let r_inj = acquire_l1c_single(
        &injected, &phasepoints(fs, n), fs, 0.0, 3000.0, 1, period,
    );
    let pure_noise: Vec<f32> = (0..n).map(|i| ((i as f64 * 1.3).sin() * 0.05) as f32).collect();
    let r_noise = acquire_l1c_single(
        &pure_noise, &phasepoints(fs, n), fs, 0.0, 3000.0, 1, period,
    );
    assert!(r_inj.signal_strength > 10.0 * r_noise.signal_strength);
}
```

- [ ] **Step 3: Run the new tests**

Run: `cargo test --release beidou l1c`
Expected: PASS. (If the BOC codes' strong autocorrelation makes `synth_signal` recovery fail due to the delay being into the secondary-lobe region, increase the injected carrier amplitude by lowering `noise_amp` to `0.01` — the correlation peak is far above sidelobes.)

- [ ] **Step 4: Commit**

```bash
git add src/beidou.rs src/l1c.rs
git commit -m "test: add BeiDou B1C and GPS L1C acquisition tests"
```

---

### Task 6: SBAS + QZSS acquisition tests

**Files:**
- Modify: `src/sbas.rs` (extend existing `#[cfg(test)] mod tests`)
- Modify: `src/qzss.rs` (extend existing `#[cfg(test)] mod tests`)

**Interfaces:**
- Consumes:
  - `acquire_sbas(signal_f32, phasepoints, sampling_freq, center_freq, search_band, prn, samples_per_chunk) -> SbasAcquisitionResult`; `sbas_gold_code(samples_per_code: f64, prn, epochs) -> Vec<f64>`; `SBAS_PRN_BASE` (= 120, private — use numeric 120 in tests).
  - `acquire_qzss(signal_f32, phasepoints, sampling_freq, center_freq, search_band, prn, samples_per_chunk) -> QzssAcquisitionResult`; `qzss_gold_code(samples_per_code: f64, prn, epochs) -> Vec<f64>`; `QZSS_PRN_BASE` (= 184, private — use numeric 184 in tests).
  - `crate::testutil::*`.
- Produces: nothing for later tasks.

- [ ] **Step 1: Add SBAS recovery + noise-contrast tests to `src/sbas.rs`**

Append to the existing test module:

```rust
use crate::testutil::{expected_recovered_sample, phasepoints, synth_signal};

#[test]
fn test_acquire_sbas_recovers_delay() {
    let fs = 1.023e6; // 1023 samples/ms (1 sample per chip)
    let period = (fs / 1000.0) as usize; // 1023
    let code = sbas_gold_code(period as f64, 120, 2.0);

    let delay = 400usize;
    let sig = synth_signal(&code, period, delay, 0.0, fs, 0.05, 41);
    let r = acquire_sbas(
        &sig, &phasepoints(fs, sig.len()), fs, 0.0, 3000.0, 120, period,
    );
    let recovered = (r.codephase_frac * period as f64).round() as usize;
    assert_eq!(recovered, expected_recovered_sample(period, delay));
    assert!(r.codephase_frac >= 0.0 && r.codephase_frac < 1.0);
}

#[test]
fn test_acquire_sbas_noise_contrast() {
    let fs = 1.023e6;
    let period = (fs / 1000.0) as usize;
    let code = sbas_gold_code(period as f64, 121, 2.0);
    let n = code.len();
    let injected = synth_signal(&code, period, 100, 0.0, fs, 0.05, 42);
    let r_inj = acquire_sbas(
        &injected, &phasepoints(fs, n), fs, 0.0, 3000.0, 121, period,
    );
    let pure_noise: Vec<f32> = (0..n).map(|i| ((i as f64 * 2.1).sin() * 0.05) as f32).collect();
    let r_noise = acquire_sbas(
        &pure_noise, &phasepoints(fs, n), fs, 0.0, 3000.0, 121, period,
    );
    assert!(r_inj.signal_strength > 10.0 * r_noise.signal_strength);
}
```

- [ ] **Step 2: Add QZSS tests to `src/qzss.rs`**

Append to the existing test module:

```rust
use crate::testutil::{expected_recovered_sample, phasepoints, synth_signal};

#[test]
fn test_acquire_qzss_recovers_delay() {
    let fs = 1.023e6;
    let period = (fs / 1000.0) as usize; // 1023
    let code = qzss_gold_code(period as f64, 184, 2.0);

    let delay = 250usize;
    let sig = synth_signal(&code, period, delay, 0.0, fs, 0.05, 51);
    let r = acquire_qzss(
        &sig, &phasepoints(fs, sig.len()), fs, 0.0, 3000.0, 184, period,
    );
    let recovered = (r.codephase_frac * period as f64).round() as usize;
    assert_eq!(recovered, expected_recovered_sample(period, delay));
    assert!(r.codephase_frac >= 0.0 && r.codephase_frac < 1.0);
}

#[test]
fn test_acquire_qzss_noise_contrast() {
    let fs = 1.023e6;
    let period = (fs / 1000.0) as usize;
    let code = qzss_gold_code(period as f64, 184, 2.0);
    let n = code.len();
    let injected = synth_signal(&code, period, 300, 0.0, fs, 0.05, 52);
    let r_inj = acquire_qzss(
        &injected, &phasepoints(fs, n), fs, 0.0, 3000.0, 184, period,
    );
    let pure_noise: Vec<f32> = (0..n).map(|i| ((i as f64 * 3.9).sin() * 0.05) as f32).collect();
    let r_noise = acquire_qzss(
        &pure_noise, &phasepoints(fs, n), fs, 0.0, 3000.0, 184, period,
    );
    assert!(r_inj.signal_strength > 10.0 * r_noise.signal_strength);
}
```

- [ ] **Step 3: Run the new tests**

Run: `cargo test --release sbas qzss`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add src/sbas.rs src/qzss.rs
git commit -m "test: add SBAS and QZSS acquisition tests"
```

---

### Task 7: `config.rs` parsing tests

**Files:**
- Modify: `src/config.rs` (add a `#[cfg(test)] mod tests`)

**Interfaces:**
- Consumes: `Config::from_json(&str) -> Result<Config, serde_json::Error>`; accessors `num_antenna()`, `sampling_frequency()`.
- Produces: nothing for later tasks.

- [ ] **Step 1: Add tests**

Append to `src/config.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_valid() {
        let c = Config::from_json(
            r#"{"num_antenna":24,"sampling_frequency":20000000.0}"#,
        )
        .unwrap();
        assert_eq!(c.num_antenna(), 24);
        assert_eq!(c.sampling_frequency(), 20_000_000.0);
    }

    #[test]
    fn test_parse_alias_fields() {
        // TART config may use the alias names.
        let c = Config::from_json(
            r#"{"num_antenna":4,"sampling_freq":16e6}"#,
        )
        .unwrap();
        assert_eq!(c.num_antenna(), 4);
        assert_eq!(c.sampling_frequency(), 16e6);
    }

    #[test]
    fn test_parse_missing_required_field_fails() {
        assert!(Config::from_json(r#"{"sampling_frequency":16e6}"#).is_err());
        assert!(Config::from_json(r#"{"num_antenna":2}"#).is_err());
    }

    #[test]
    fn test_parse_malformed_json_fails() {
        assert!(Config::from_json("not json").is_err());
        assert!(Config::from_json("").is_err());
    }
}
```

- [ ] **Step 2: Run the new tests**

Run: `cargo test --release config`
Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git add src/config.rs
git commit -m "test: add Config::from_json parsing tests"
```

---

### Task 8: `main.rs` refactor (arg parsing + MAD filtering) with tests

**Files:**
- Modify: `src/main.rs`

**Interfaces:**
- Produces (used by `main` and by this task's tests):
  - `struct ParsedArgs` with fields mirroring all CLI flags (see below).
  - `fn parse_args(args: &[String]) -> Result<ParsedArgs, String>`
  - `fn print_usage(prog: &str)`
  - `trait MadFilterable { fn phase_mad(&self) -> Option<f64>; fn freq_mad(&self) -> Option<f64>; }` implemented for all six PRN-result types.
  - `fn apply_mad_filters<T: MadFilterable>(results: &mut Vec<T>, phase_thresh: Option<f64>, freq_thresh: Option<f64>) -> usize`

**Refactor must preserve CLI behavior.**

- [ ] **Step 1: Add the `ParsedArgs` type, `parse_args`, and `print_usage`**

Insert before `fn main()` in `src/main.rs`:

```rust
/// All CLI options parsed from the command line.
#[derive(Debug, Default)]
struct ParsedArgs {
    file: Option<String>,
    antenna_i: Option<usize>,
    antenna_j: Option<usize>,
    gps_flag: bool,
    galileo_flag: bool,
    beidou_flag: bool,
    sbas_flag: bool,
    l1c_flag: bool,
    qzss_flag: bool,
    all_flag: bool,
    ant_list: Option<Vec<usize>>,
    filter_phase_mad: Option<f64>,
    filter_freq_mad: Option<f64>,
    output_file: Option<String>,
    debug_flag: bool,
    cn0_flag: bool,
    benchmark_flag: bool,
    prn_filter: Option<Vec<usize>>,
    version_flag: bool,
}

/// Read the value following a flag, advancing `i`, or return an error if absent.
fn flag_value<'a>(args: &'a [String], i: &mut usize, flag: &str) -> Result<&'a str, String> {
    if *i + 1 >= args.len() {
        return Err(format!("missing value for {flag}"));
    }
    *i += 1;
    Ok(args[*i].as_str())
}

/// Parse command-line args (excluding the program name in `args[0]`), matching
/// the original `main` loop exactly. Returns an error string on bad input.
fn parse_args(args: &[String]) -> Result<ParsedArgs, String> {
    let mut p = ParsedArgs::default();
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--file" => {
                let v = flag_value(args, &mut i, "--file")?;
                p.file = Some(v.to_string());
            }
            "--i" => {
                let v = flag_value(args, &mut i, "--i")?;
                p.antenna_i = Some(v.parse().map_err(|_| format!("invalid integer for --i: {v}"))?);
            }
            "--j" => {
                let v = flag_value(args, &mut i, "--j")?;
                p.antenna_j = Some(v.parse().map_err(|_| format!("invalid integer for --j: {v}"))?);
            }
            "--gps" => p.gps_flag = true,
            "--galileo" => p.galileo_flag = true,
            "--beidou" => p.beidou_flag = true,
            "--sbas" => p.sbas_flag = true,
            "--l1c" => p.l1c_flag = true,
            "--qzss" => p.qzss_flag = true,
            "--all" => p.all_flag = true,
            "--ant" => {
                let v = flag_value(args, &mut i, "--ant")?;
                p.ant_list = Some(
                    v.split(',')
                        .map(|s| s.trim().parse().map_err(|_| format!("invalid integer in --ant: {s}")))
                        .collect::<Result<Vec<_>, _>>()?,
                );
            }
            "--filter-phase-mad" => {
                let v = flag_value(args, &mut i, "--filter-phase-mad")?;
                p.filter_phase_mad = Some(v.parse().map_err(|_| format!("invalid float for --filter-phase-mad: {v}"))?);
            }
            "--filter-freq-mad" => {
                let v = flag_value(args, &mut i, "--filter-freq-mad")?;
                p.filter_freq_mad = Some(v.parse().map_err(|_| format!("invalid float for --filter-freq-mad: {v}"))?);
            }
            "--output" => {
                let v = flag_value(args, &mut i, "--output")?;
                p.output_file = Some(v.to_string());
            }
            "--debug" => p.debug_flag = true,
            "--cn0" => p.cn0_flag = true,
            "--benchmark" => p.benchmark_flag = true,
            "--prn" => {
                let v = flag_value(args, &mut i, "--prn")?;
                p.prn_filter = Some(
                    v.split(',')
                        .map(|s| s.trim().parse().map_err(|_| format!("invalid integer in --prn: {s}")))
                        .collect::<Result<Vec<_>, _>>()?,
                );
            }
            "--version" => p.version_flag = true,
            other => return Err(format!("unknown argument: {other}")),
        }
        i += 1;
    }
    Ok(p)
}

fn print_usage(prog: &str) {
    eprintln!(
        "usage: {prog} --file <observation.hdf> [--i <i> --j <j>] [--all] [--gps] [--galileo] [--beidou] [--sbas] [--l1c] [--qzss] [--cn0] [--prn <a,b,...>] [--ant <a,b,...>] [--filter-phase-mad <x>] [--filter-freq-mad <x>] [--output <path>] [--debug] [--benchmark]"
    );
}
```

- [ ] **Step 2: Rewrite the top of `main()` to use `parse_args`**

Replace the entire manual `let args...` … `if i >= n_ant` argument-handling preamble at the top of `main` down to (but not including) the `let obs = ...` statement. The new top of `main`:

```rust
fn main() {
    let args: Vec<String> = std::env::args().collect();
    let parsed = match parse_args(&args) {
        Ok(p) => p,
        Err(e) => {
            eprintln!("{e}");
            print_usage(&args[0]);
            std::process::exit(1);
        }
    };

    if parsed.version_flag {
        println!("tart-gnss-acquire v{}", env!("CARGO_PKG_VERSION"));
        std::process::exit(0);
    }

    let ParsedArgs {
        file,
        antenna_i,
        antenna_j,
        mut gps_flag,
        mut galileo_flag,
        mut beidou_flag,
        mut sbas_flag,
        mut l1c_flag,
        mut qzss_flag,
        all_flag,
        ant_list,
        filter_phase_mad,
        filter_freq_mad,
        output_file,
        debug_flag,
        cn0_flag,
        benchmark_flag,
        prn_filter,
        version_flag: _,
    } = parsed;

    // --all implies all six acquisition modes
    if all_flag {
        gps_flag = true;
        galileo_flag = true;
        beidou_flag = true;
        sbas_flag = true;
        l1c_flag = true;
        qzss_flag = true;
    }

    let file = if benchmark_flag {
        file
    } else {
        Some(file.unwrap_or_else(|| {
            eprintln!("missing --file (or use --benchmark)");
            print_usage(&args[0]);
            std::process::exit(1);
        }))
    };
```

The rest of `main` (from `let obs = ...` onward) stays unchanged, except the MAD-filter block is replaced in Step 3.

- [ ] **Step 3: Extract `MadFilterable` + `apply_mad_filters` and replace the filter block**

Add the trait and helper (place near the other helpers before `main`):

```rust
/// Common accessors for MAD-based filtering across PRN result types.
trait MadFilterable {
    fn phase_mad(&self) -> Option<f64>;
    fn freq_mad(&self) -> Option<f64>;
}

macro_rules! impl_mad_filterable {
    ($ty:ty) => {
        impl MadFilterable for $ty {
            fn phase_mad(&self) -> Option<f64> { self.phase_mad }
            fn freq_mad(&self) -> Option<f64> { self.freq_mad }
        }
    };
}

impl_mad_filterable!(acquisition::GpsPrnResult);
impl_mad_filterable!(galileo::GalileoPrnResult);
impl_mad_filterable!(beidou::BeiDouPrnResult);
impl_mad_filterable!(sbas::SbasPrnResult);
impl_mad_filterable!(l1c::L1CPrnResult);
impl_mad_filterable!(qzss::QzssPrnResult);

/// Retain only results whose MADs are at or below the given thresholds.
/// A `None` MAD (single antenna) always passes. Returns the number removed.
fn apply_mad_filters<T: MadFilterable>(
    results: &mut Vec<T>,
    phase_thresh: Option<f64>,
    freq_thresh: Option<f64>,
) -> usize {
    let before = results.len();
    if let Some(t) = phase_thresh {
        results.retain(|r| r.phase_mad().map_or(true, |m| m <= t));
    }
    if let Some(t) = freq_thresh {
        results.retain(|r| r.freq_mad().map_or(true, |m| m <= t));
    }
    before - results.len()
}
```

Replace the whole `if filter_phase_mad.is_some() || filter_freq_mad.is_some() { ... }` filter block in `main` with:

```rust
        // --- Apply MAD filters ---------------------------------------------
        if filter_phase_mad.is_some() || filter_freq_mad.is_some() {
            let mut filter_count = 0u64;
            if let Some(o) = output.gps.as_mut() {
                filter_count += apply_mad_filters(&mut o.results, filter_phase_mad, filter_freq_mad) as u64;
            }
            if let Some(o) = output.galileo.as_mut() {
                filter_count += apply_mad_filters(&mut o.results, filter_phase_mad, filter_freq_mad) as u64;
            }
            if let Some(o) = output.beidou.as_mut() {
                filter_count += apply_mad_filters(&mut o.results, filter_phase_mad, filter_freq_mad) as u64;
            }
            if let Some(o) = output.sbas.as_mut() {
                filter_count += apply_mad_filters(&mut o.results, filter_phase_mad, filter_freq_mad) as u64;
            }
            if let Some(o) = output.l1c.as_mut() {
                filter_count += apply_mad_filters(&mut o.results, filter_phase_mad, filter_freq_mad) as u64;
            }
            if let Some(o) = output.qzss.as_mut() {
                filter_count += apply_mad_filters(&mut o.results, filter_phase_mad, filter_freq_mad) as u64;
            }
            if filter_count > 0 {
                eprintln!("Filtered out {filter_count} results (MAD thresholds)");
            }
        }
```

- [ ] **Step 4: Verify the refactor compiles and the existing binary still builds**

Run: `cargo build --release`
Expected: compiles cleanly with no warnings about unused struct fields.

- [ ] **Step 5: Add `parse_args` tests**

Append a test module to `src/main.rs`:

```rust
#[cfg(test)]
mod tests {
    use super::*;

    fn a(list: &[&str]) -> Vec<String> {
        list.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn test_parse_all_flags() {
        let p = parse_args(&a(&[
            "--file", "obs.hdf", "--all", "--cn0", "--debug",
            "--prn", "1,5,12", "--ant", "0,2",
            "--filter-phase-mad", "0.002",
            "--filter-freq-mad", "100.0",
            "--output", "out.json",
        ]))
        .unwrap();
        assert_eq!(p.file.as_deref(), Some("obs.hdf"));
        assert_eq!(p.prn_filter.as_deref(), Some(&[1usize, 5, 12][..]));
        assert_eq!(p.ant_list.as_deref(), Some(&[0usize, 2][..]));
        assert_eq!(p.filter_phase_mad, Some(0.002));
        assert_eq!(p.filter_freq_mad, Some(100.0));
        assert_eq!(p.output_file.as_deref(), Some("out.json"));
        assert!(p.cn0_flag && p.debug_flag);
    }

    #[test]
    fn test_parse_all_sets_flags() {
        let p = parse_args(&a(&["--all", "--file", "x.hdf"])).unwrap();
        assert!(p.all_flag);
    }

    #[test]
    fn test_parse_unknown_arg_err() {
        assert!(parse_args(&a(&["--nope"])).is_err());
        assert!(parse_args(&a(&["--file"])).is_err()); // missing value
        assert!(parse_args(&a(&["--i", "bogus"])).is_err());
        assert!(parse_args(&a(&["--filter-freq-mad", "x"])).is_err());
    }

    #[test]
    fn test_parse_version() {
        let p = parse_args(&a(&["--version"])).unwrap();
        assert!(p.version_flag);
    }
}
```


- [ ] **Step 6: Add `apply_mad_filters` tests**

Same test module:

```rust
    #[test]
    fn test_apply_mad_filters_threshold() {
        // Construct minimal GPS results with MADs.
        let mut results = vec![
            acquisition::GpsPrnResult { sv: "GPS01".into(), strengths: vec![], phases: vec![], freqs: vec![],
                cn0_acr: None, phase_median: Some(0.1), phase_mad: Some(0.001),
                freq_median: Some(0.0), freq_mad: Some(50.0) },
            acquisition::GpsPrnResult { sv: "GPS02".into(), strengths: vec![], phases: vec![], freqs: vec![],
                cn0_acr: None, phase_median: Some(0.2), phase_mad: Some(0.005),
                freq_median: Some(0.0), freq_mad: Some(20.0) },
            acquisition::GpsPrnResult { sv: "GPS03".into(), strengths: vec![], phases: vec![], freqs: vec![],
                cn0_acr: None, phase_median: None, phase_mad: None,
                freq_median: None, freq_mad: None },
        ];
        // phase threshold 0.002 keeps GPS01 (0.001) and drops GPS02 (0.005);
        // GPS03 has None MAD and must be kept.
        let removed = apply_mad_filters(&mut results, Some(0.002), None);
        assert_eq!(removed, 1);
        assert_eq!(results.len(), 2);
        assert!(results.iter().all(|r| r.sv != "GPS02"));
        assert!(results.iter().any(|r| r.sv == "GPS03"));
    }

    #[test]
    fn test_apply_mad_filters_freq() {
        let mut results = vec![
            acquisition::GpsPrnResult { sv: "GPS01".into(), strengths: vec![], phases: vec![], freqs: vec![],
                cn0_acr: None, phase_median: Some(0.0), phase_mad: Some(0.0),
                freq_median: Some(0.0), freq_mad: Some(100.0) },
            acquisition::GpsPrnResult { sv: "GPS02".into(), strengths: vec![], phases: vec![], freqs: vec![],
                cn0_acr: None, phase_median: Some(0.0), phase_mad: Some(0.0),
                freq_median: Some(0.0), freq_mad: Some(90.0) },
        ];
        let removed = apply_mad_filters(&mut results, None, Some(95.0));
        assert_eq!(removed, 1);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].sv, "GPS02");
    }

    #[test]
    fn test_apply_mad_filters_no_threshold_noop() {
        let mut results = vec![acquisition::GpsPrnResult { sv: "GPS01".into(), strengths: vec![], phases: vec![],
            freqs: vec![], cn0_acr: None, phase_median: Some(0.0), phase_mad: Some(1000.0),
            freq_median: Some(0.0), freq_mad: Some(1000.0) }];
        assert_eq!(apply_mad_filters(&mut results, None, None), 0);
        assert_eq!(results.len(), 1);
        assert_eq!(results[0].sv, "GPS01");
    }
```

- [ ] **Step 7: Run the full test suite**

Run: `cargo test --release`
Expected: all tests (existing 67 + new) PASS.

- [ ] **Step 8: Commit**

```bash
git add src/main.rs
git commit -m "refactor: extract testable arg parsing and MAD filtering in main, add tests"
```

---

### Task 9: `observation.rs` additional tests

**Files:**
- Modify: `src/observation.rs` (extend existing `#[cfg(test)] mod tests`)

**Interfaces:**
- Consumes: `Observation::new`, `Observation::random(num_antenna, sampling_frequency, num_samples)`, `get_antenna(usize) -> Vec<f64>`, `get_means() -> Vec<f64>`, `get_sampling_rate()`, `correlate(i, j) -> f64`; `Config::from_json`.
- Produces: nothing for later tasks.

- [ ] **Step 1: Add tests to the existing `observation.rs` test module**

```rust
    #[test]
    fn test_get_antenna_bipolar_mapping() {
        let data = vec![vec![0, 1, 0, 1], vec![1, 1, 0, 0]];
        let cfg = Config::from_json(r#"{"num_antenna":2,"sampling_frequency":16e6}"#).unwrap();
        let obs = Observation::new(chrono::Utc::now(), cfg, data);
        assert_eq!(obs.get_antenna(0), vec![-1.0, 1.0, -1.0, 1.0]);
        assert_eq!(obs.get_antenna(1), vec![1.0, 1.0, -1.0, -1.0]);
        assert_eq!(obs.get_sampling_rate(), 16e6);
    }

    #[test]
    #[should_panic]
    fn test_get_antenna_out_of_range_panics() {
        let data = vec![vec![0u8, 1]];
        let cfg = Config::from_json(r#"{"num_antenna":1,"sampling_frequency":16e6}"#).unwrap();
        let obs = Observation::new(chrono::Utc::now(), cfg, data);
        obs.get_antenna(1);
    }

    #[test]
    fn test_get_means() {
        // antenna 0: [0,0,0,1] -> sum 1, avg 0.25, *2-1 = -0.5
        // antenna 1: [1,1,1,1] -> avg 1.0, *2-1 = 1.0
        let data = vec![vec![0u8, 0, 0, 1], vec![1, 1, 1, 1]];
        let cfg = Config::from_json(r#"{"num_antenna":2,"sampling_frequency":16e6}"#).unwrap();
        let obs = Observation::new(chrono::Utc::now(), cfg, data);
        let means = obs.get_means();
        assert_eq!(means.len(), 2);
        assert!((means[0] - (-0.5)).abs() < 1e-12);
        assert!((means[1] - 1.0).abs() < 1e-12);
    }

    #[test]
    #[should_panic]
    fn test_correlate_length_mismatch_panics() {
        // Different sample counts -> panic on length mismatch.
        let data = vec![vec![0u8, 1, 0], vec![1, 1]];
        let cfg = Config::from_json(r#"{"num_antenna":2,"sampling_frequency":16e6}"#).unwrap();
        let obs = Observation::new(chrono::Utc::now(), cfg, data);
        obs.correlate(0, 1);
    }

    #[test]
    fn test_random_dimensions() {
        let obs = Observation::random(3, 16e6, 1000);
        assert_eq!(obs.config.num_antenna(), 3);
        assert_eq!(obs.data.len(), 3);
        for ant in &obs.data {
            assert_eq!(ant.len(), 1000);
            assert!(ant.iter().all(|&b| b == 0 || b == 1));
        }
        assert_eq!(obs.get_sampling_rate(), 16e6);
    }
```

- [ ] **Step 2: Run the new tests**

Run: `cargo test --release observation`
Expected: PASS (the three `#[should_panic]` tests pass only if the corresponding panics occur).

- [ ] **Step 3: Commit**

```bash
git add src/observation.rs
git commit -m "test: extend observation tests for means, antenna mapping, and edge cases"
```

---

### Task 10: Final verification

**Files:**
- None (verification only).

- [ ] **Step 1: Run the full test suite in release mode**

Run: `cargo test --release`
Expected: all tests pass (original 67 + all new tests).

- [ ] **Step 2: Run clippy (if available)**

Run: `cargo clippy --all-targets -- -D warnings`
Expected: no errors. If the machine lacks clippy, skip; do not change behavior to satisfy it.

- [ ] **Step 3: Run a smoke check of the CLI still works**

Run:
```bash
cargo run --release -- --file data/observation.hdf --gps --prn 3 --ant 0
```
Expected: prints a `GPS03` result JSON as before the refactor.

- [ ] **Step 4: Confirm clean worktree and log**

Run: `git status` (clean) and `git log --oneline -10` to review the test commits.

---

## Self-Review

**Spec coverage:**
- `correlate.rs` (core) → Task 2 ✓
- `acquire_full` + all six `acquire_*_single` + `acquire_all_*` → Tasks 3–6 ✓
- `config.rs` parsing → Task 7 ✓
- `main.rs` arg parsing + MAD filtering (with refactor) → Task 8 ✓
- `observation.rs` edge cases → Task 9 ✓
- Shared synthetic-signal helper → Task 1 ✓
- Deterministic, small-FFT tests; `cargo test` verification per task and Task 10 ✓

**Placeholder scan:** Every step contains concrete code; no TBD/TODO. Numeric conventions verified empirically and documented in Global Constraints.

**Type consistency:** `synth_signal`/`phasepoints`/`expected_recovered_sample` signatures used identically in Tasks 2–6. `parse_args`, `ParsedArgs`, `apply_mad_filters`, and `MadFilterable` signatures match between Step 1–3 and the Step 5–6 tests in Task 8.
