// Copyright (c) 2026 Tim Molteno <tim@elec.ac.nz>
// SPDX-License-Identifier: GPL-3.0

//! BeiDou B1C primary codes (pilot component, B1Cp).
//!
//! 10230-chip Weil codes derived from a Legendre sequence of length 10243.
//! The pilot component uses BOC(1,1) modulation for acquisition.
//!
//! Reference implementation: Peter Monta's GNSS-DSP-tools
//!   https://github.com/pmonta/GNSS-DSP-tools/blob/master/gnsstools/beidou/b1cp.py

use std::sync::LazyLock;

/// Number of BeiDou B1C PRN codes available (SV 1–63).
pub const BEIDOU_B1C_NUM_SATS: usize = 63;

/// Chips per BeiDou B1C primary code period.
pub const BEIDOU_B1C_CHIPS: usize = 10_230;

/// Legendre sequence length (prime number).
const LEGENDRE_N: usize = 10_243;

/// Chip rate (Hz).
#[allow(dead_code)]
pub const BEIDOU_B1C_CHIP_RATE: f64 = 1.023e6;

/// Code period in seconds (10 ms).
pub const BEIDOU_B1C_CODE_PERIOD: f64 = 0.010;

/// BOC(1,1) subcarrier frequency (Hz).
#[allow(dead_code)]
pub const BOC_SUBCARRIER_FREQ: f64 = 1.023e6;

// ---------------------------------------------------------------------------
// B1Cp code parameters (Weil index w, phase offset p)
// ---------------------------------------------------------------------------

/// (w, p) parameters for the B1C pilot component, PRNs 1–63.
///
/// `w` — Weil shift index for the Legendre-sequence-based Weil code.
/// `p` — phase offset: the code is extracted from the Weil sequence starting
///       at index `p-1` and wrapping modulo N=10243.
///
/// From GNSS-DSP-tools `gnsstools/beidou/b1cp.py`.
#[rustfmt::skip]
const B1CP_PARAMS: [(usize, usize); BEIDOU_B1C_NUM_SATS] = [
    ( 796,  7575), ( 156,  2369), (4198,  5688), (3941,   539),
    (1374,  2270), (1338,  7306), (1833,  6457), (2521,  6254),
    (3175,  5644), ( 168,  7119), (2715,  1402), (4408,  5557),
    (3160,  5764), (2796,  1073), ( 459,  7001), (3594,  5910),
    (4813, 10060), ( 586,  2710), (1428,  1546), (2371,  6887),
    (2285,  1883), (3377,  5613), (4965,  5062), (3779,  1038),
    (4547, 10170), (1646,  6484), (1430,  1718), ( 607,  2535),
    (2118,  1158), (4709,   526), (1149,  7331), (3283,  5844),
    (2473,  6423), (1006,  6968), (3670,  1280), (1817,  1838),
    ( 771,  1989), (2173,  6468), ( 740,  2091), (1433,  1581),
    (2458,  1453), (3459,  6252), (2155,  7122), (1205,  7711),
    ( 413,  7216), ( 874,  2113), (2463,  1095), (1106,  1628),
    (1590,  1713), (3873,  6102), (4026,  6123), (4272,  6070),
    (3556,  1115), ( 128,  8047), (1200,  6795), ( 130,  2575),
    (4494,    53), (1871,  1729), (3073,  6388), (4386,   682),
    (4098,  5565), (1923,  7160), (1176,  2277),
];

// ---------------------------------------------------------------------------
// Legendre sequence generation
// ---------------------------------------------------------------------------

/// Generate the Legendre sequence of length `LEGENDRE_N` (10243).
///
/// For prime N:
/// - L(0) = 0
/// - L(i) = 1 if i is a quadratic residue modulo N, else 0
fn legendre_sequence() -> [u8; LEGENDRE_N] {
    let mut seq = [0u8; LEGENDRE_N];
    let mut residues = [false; LEGENDRE_N];
    for x in 1..LEGENDRE_N {
        let idx = (x * x) % LEGENDRE_N;
        residues[idx] = true;
    }
    // seq[0] stays 0
    for i in 1..LEGENDRE_N {
        seq[i] = if residues[i] { 1 } else { 0 };
    }
    seq
}

/// Cached Legendre sequence — computed once, shared by all PRNs.
static LEGENDRE: LazyLock<[u8; LEGENDRE_N]> = LazyLock::new(legendre_sequence);

/// Cached BOC(1,1)-modulated pilot codes for all 63 PRNs.
/// Pre-computed once at first use; ~10 MB total.
static B1C_BOC_CODES: LazyLock<Vec<Vec<f64>>> = LazyLock::new(|| {
    (1..=BEIDOU_B1C_NUM_SATS)
        .map(generate_b1c_pilot_code)
        .collect()
});

// ---------------------------------------------------------------------------
// Weil code generation
// ---------------------------------------------------------------------------

/// Generate the raw B1Cp chips (0/1 values) for a given PRN (1-based).
///
/// Returns a `[u8; BEIDOU_B1C_CHIPS]` array of 0/1 chips matching the
/// GNSS-DSP-tools `b1cp_code(prn)` output.
///
/// Algorithm:
///   1. Compute Legendre sequence L of length N=10243.
///   2. For each k in 0..N-1, Weil code W[k] = L[k] XOR L[(k+w) mod N].
///   3. Extract code_length=10230 chips starting at offset (p-1):
///      c[n] = W[(n + p - 1) mod N],  for n = 0..code_length-1.
fn generate_b1cp_chips(prn: usize) -> [u8; BEIDOU_B1C_CHIPS] {
    let (w, p) = B1CP_PARAMS[prn - 1];
    let legendre = &*LEGENDRE;
    let mut chips = [0u8; BEIDOU_B1C_CHIPS];

    for (n, chip) in chips.iter_mut().enumerate() {
        let k = (n + p - 1) % LEGENDRE_N;
        let l_k = legendre[k];
        let l_shifted = legendre[(k + w) % LEGENDRE_N];
        *chip = l_k ^ l_shifted;
    }

    chips
}

/// Generate the ±1 B1Cp code chips for a given PRN (1-based).
///
/// Maps 0→+1, 1→-1, matching the `1.0 - 2.0*x` convention in GNSS-DSP-tools.
pub fn generate_b1c_code(prn: usize) -> [f64; BEIDOU_B1C_CHIPS] {
    let chips = generate_b1cp_chips(prn);
    let mut code = [0.0f64; BEIDOU_B1C_CHIPS];
    for (i, &chip) in chips.iter().enumerate() {
        code[i] = if chip == 0 { 1.0 } else { -1.0 };
    }
    code
}

/// Generate the BOC(1,1)-modulated B1C pilot code for a given PRN (1-based).
///
/// Each chip is split into two half-chips: [+chip, -chip], giving 20460
/// samples per code period.  The BOC(1,1) envelope is `[+1, -1]` per chip,
/// which is equivalent to GNSS-DSP-tools `boc11 = [1.0, -1.0]`.
pub fn generate_b1c_pilot_code(prn: usize) -> Vec<f64> {
    assert!(
        (1..=BEIDOU_B1C_NUM_SATS).contains(&prn),
        "BeiDou PRN must be in 1..={BEIDOU_B1C_NUM_SATS}, got {prn}"
    );
    let weil = generate_b1c_code(prn);

    let mut boc = Vec::with_capacity(BEIDOU_B1C_CHIPS * 2);
    for &chip in &weil {
        boc.push(chip);
        boc.push(-chip);
    }
    boc
}

/// Resample the BOC(1,1)-modulated B1C pilot code to exactly `num_samples`
/// samples.
///
/// The code has `samples_per_code` samples per 10-ms code period; it is
/// repeated (with a possible partial final period) to fill exactly
/// `num_samples` samples, so its length always matches the FFT size.
pub fn b1c_code_resampled(samples_per_code: f64, prn: usize, num_samples: usize) -> Vec<f64> {
    let boc = &B1C_BOC_CODES[prn - 1];
    let boc_len = boc.len() as f64; // 20460
    let samples_per_chip = samples_per_code / boc_len;

    (0..num_samples)
        .map(|n| {
            let idx = ((n as f64 / samples_per_chip).floor() as usize) % boc.len();
            boc[idx]
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Testbench — cross-validated against GNSS-DSP-tools reference
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- Legendre sequence -------------------------------------------------

    #[test]
    fn test_legendre_length() {
        let seq = legendre_sequence();
        assert_eq!(seq.len(), LEGENDRE_N);
    }

    #[test]
    fn test_legendre_zero() {
        let seq = legendre_sequence();
        assert_eq!(seq[0], 0);
    }

    #[test]
    fn test_legendre_residue_count() {
        let seq = legendre_sequence();
        let ones: usize = seq.iter().map(|&b| b as usize).sum();
        assert_eq!(ones, (LEGENDRE_N - 1) / 2);
    }

    // --- B1Cp chips (0/1) cross-checked against Python reference ----------

    /// Convert 24 chips (0/1, MSB-first, groups of 3) to an octal string.
    /// Matches the Python `chips2octal` function in GNSS-DSP-tools.
    fn chips24_to_octal(chips: &[u8]) -> String {
        assert_eq!(chips.len(), 24);
        let mut s = String::with_capacity(8);
        for i in 0..8 {
            let d = 4 * chips[3 * i] + 2 * chips[3 * i + 1] + chips[3 * i + 2];
            s.push(char::from_digit(d as u32, 8).unwrap());
        }
        s
    }

    /// Reference octal strings for the first 24 chips of each PRN 1..=63,
    /// generated by GNSS-DSP-tools `acquire-beidou-b1cp.py` code generator.
    #[rustfmt::skip]
    const REF_FIRST24: [&str; BEIDOU_B1C_NUM_SATS] = [
        "71676756", "60334021", "24562714", "61011650", "67337730", "23762642",
        "25365366", "57226722", "72643175", "00236125", "12071371", "61136116",
        "36261215", "13607013", "31010541", "73163062", "30250537", "56226421",
        "26205736", "02450570", "66511327", "06323465", "10633350", "10544206",
        "43714115", "55641056", "26572456", "75123401", "70041254", "53034467",
        "50733517", "73077145", "55454316", "37137206", "45724432", "55560467",
        "13467065", "24245150", "22265044", "10003471", "36537736", "57706617",
        "76411007", "61643153", "50125760", "66657234", "01350500", "43621551",
        "42435620", "74327566", "44553226", "52231514", "46576047", "46312270",
        "04717127", "50407031", "10044104", "36610123", "73470741", "24072445",
        "07765425", "32242545", "03210227",
    ];

    /// Reference octal strings for the last 24 chips of each PRN 1..=63.
    #[rustfmt::skip]
    const REF_LAST24: [&str; BEIDOU_B1C_NUM_SATS] = [
        "13053205", "46604773", "60007065", "23616424", "66243127", "33630334",
        "43456307", "76521063", "52465264", "76142064", "60232627", "05607727",
        "77737367", "16031533", "55416670", "33076260", "73355574", "42437243",
        "66470710", "54366756", "23666556", "74622250", "16402734", "54230354",
        "37167223", "56136734", "62211315", "40615033", "63213062", "03066540",
        "30062510", "34360276", "45431517", "47647044", "33773217", "77620561",
        "17327352", "62223375", "67665257", "27515010", "37705710", "76736116",
        "77202566", "25334277", "70220333", "22376763", "31043217", "20166102",
        "16423062", "31245527", "37160613", "03414402", "04003162", "54703562",
        "25225202", "31643432", "27063234", "40756155", "24774305", "51507057",
        "12225744", "62104320", "56250500",
    ];

    #[test]
    fn test_b1cp_first24_octal() {
        for prn in 1..=BEIDOU_B1C_NUM_SATS {
            let chips = generate_b1cp_chips(prn);
            let octal = chips24_to_octal(&chips[0..24]);
            assert_eq!(
                octal,
                REF_FIRST24[prn - 1],
                "PRN {prn}: first 24 chips mismatch"
            );
        }
    }

    #[test]
    fn test_b1cp_last24_octal() {
        for prn in 1..=BEIDOU_B1C_NUM_SATS {
            let chips = generate_b1cp_chips(prn);
            let octal = chips24_to_octal(&chips[BEIDOU_B1C_CHIPS - 24..]);
            assert_eq!(
                octal,
                REF_LAST24[prn - 1],
                "PRN {prn}: last 24 chips mismatch"
            );
        }
    }

    // --- B1Cp bipolar (±1) code -------------------------------------------

    #[test]
    fn test_generate_b1c_code_length() {
        let code = generate_b1c_code(1);
        assert_eq!(code.len(), BEIDOU_B1C_CHIPS);
    }

    #[test]
    fn test_generate_b1c_code_bipolar() {
        let code = generate_b1c_code(5);
        for &v in &code {
            assert!(v == 1.0 || v == -1.0, "unexpected value {v}");
        }
    }

    #[test]
    fn test_generate_b1c_pilot_code_length() {
        let code = generate_b1c_pilot_code(1);
        assert_eq!(code.len(), BEIDOU_B1C_CHIPS * 2);
    }

    #[test]
    fn test_b1c_code_resampled_length() {
        let samples_per_code = (BEIDOU_B1C_CHIPS * 2) as f64;
        let code = b1c_code_resampled(samples_per_code, 1, BEIDOU_B1C_CHIPS * 2);
        assert_eq!(code.len(), BEIDOU_B1C_CHIPS * 2);
    }

    #[test]
    fn test_all_prns_generate() {
        for prn in 1..=BEIDOU_B1C_NUM_SATS {
            let code = generate_b1c_pilot_code(prn);
            assert_eq!(
                code.len(),
                BEIDOU_B1C_CHIPS * 2,
                "PRN {prn} wrong pilot code length"
            );
            for &v in &code {
                assert!(v == 1.0 || v == -1.0, "PRN {prn}: unexpected value {v}");
            }
        }
    }

    #[test]
    fn test_all_params_in_range() {
        for (prn, &(w, p)) in B1CP_PARAMS.iter().enumerate() {
            assert!(
                w > 0 && w < LEGENDRE_N,
                "PRN {} pilot w={} out of range 1..{}",
                prn + 1,
                w,
                LEGENDRE_N - 1
            );
            assert!(
                p > 0 && p <= LEGENDRE_N,
                "PRN {} pilot p={} out of range 1..={}",
                prn + 1,
                p,
                LEGENDRE_N
            );
        }
    }

    /// Cross-check: verify 0→+1, 1→-1 mapping matches Python
    /// `1.0 - 2.0*x` convention.
    #[test]
    fn test_chip_mapping_matches_python() {
        let chips = generate_b1cp_chips(1);
        let code = generate_b1c_code(1);
        for i in 0..BEIDOU_B1C_CHIPS {
            let expected = if chips[i] == 0 { 1.0 } else { -1.0 };
            assert_eq!(
                code[i], expected,
                "chip {i}: mapping mismatch"
            );
        }
    }
}
