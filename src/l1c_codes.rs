// Copyright (c) 2026 Tim Molteno <tim@elec.ac.nz>
// SPDX-License-Identifier: GPL-3.0

//! GPS L1C primary codes (pilot component, L1Cp).
//!
//! 10230-chip Weil codes derived from a Legendre sequence of length 10223
//! with a 7-bit expansion sequence insertion.
//!
//! Reference: IS-GPS-800H, Section 3.2.2.1.1 and Table 3.2-2.

/// Number of GPS L1C PRN codes available (SV 1–63).
pub const L1C_NUM_SATS: usize = 63;

/// Chips per GPS L1C primary code period (10230 chips after insertion).
pub const L1C_CHIPS: usize = 10_230;

/// Legendre sequence length (prime number p = 10223).
const LEGENDRE_N: usize = 10_223;

/// Chip rate (Hz).
#[allow(dead_code)]
pub const L1C_CHIP_RATE: f64 = 1.023e6;

/// Code period in seconds (10 ms).
pub const L1C_CODE_PERIOD: f64 = 0.010;

/// BOC(1,1) subcarrier frequency (Hz).
#[allow(dead_code)]
const BOC_SUBCARRIER_FREQ: f64 = 1.023e6;

// ---------------------------------------------------------------------------
// L1Cp code parameters (Weil index w, insertion index p)
// ---------------------------------------------------------------------------

/// (w, p) parameters for the L1C pilot component, PRNs 1–63.
///
/// `w` — Weil shift index for the Legendre-sequence-based Weil code.
/// `p` — Insertion index (1-based). The expansion sequence is inserted
///       *before* the p-th value of the Weil code.
///
/// From IS-GPS-800H Table 3.2-2, L1CP column.
#[rustfmt::skip]
const L1CP_PARAMS: [(usize, usize); L1C_NUM_SATS] = [
    (5111,   412), (5109,   161), (5108,     1), (5106,   303),
    (5103,   207), (5101,  4971), (5100,  4496), (5098,     5),
    (5095,  4557), (5094,   485), (5093,   253), (5091,  4676),
    (5090,     1), (5081,    66), (5080,  4485), (5069,   282),
    (5068,   193), (5054,  5211), (5044,   729), (5027,  4848),
    (5026,   982), (5014,  5955), (5004,  9805), (4980,   670),
    (4915,   464), (4909,    29), (4893,   429), (4885,   394),
    (4832,   616), (4824,  9457), (4591,  4429), (3706,  4771),
    (5092,   365), (4986,  9705), (4965,  9489), (4920,  4193),
    (4917,  9947), (4858,   824), (4847,   864), (4790,   347),
    (4770,   677), (4318,  6544), (4126,  6312), (3961,  9804),
    (3790,   278), (4911,  9461), (4881,   444), (4827,  4839),
    (4795,  4144), (4789,  9875), (4725,   197), (4675,  1156),
    (4539,  4674), (4535, 10035), (4458,  4504), (4197,     5),
    (4096,  9937), (3484,   430), (3481,     5), (3393,   355),
    (3175,   909), (2360,  1622), (1852,  6284),
];

/// Expansion sequence inserted into the Weil code.
const EXPANSION_SEQUENCE: [u8; 7] = [0, 1, 1, 0, 1, 0, 0];

// ---------------------------------------------------------------------------
// Legendre sequence generation
// ---------------------------------------------------------------------------

/// Generate the Legendre sequence of length `LEGENDRE_N` (10223).
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

// ---------------------------------------------------------------------------
// Weil code + expansion insertion
// ---------------------------------------------------------------------------

/// Generate the raw L1Cp chips (0/1 values) for a given PRN (1-based).
///
/// Returns a `[u8; L1C_CHIPS]` array of 0/1 chips.
///
/// Algorithm (IS-GPS-800H Section 3.2.2.1.1):
///   1. Compute Legendre sequence L of length N=10223.
///   2. Weil code: W[k] = L[k] XOR L[(k + w) mod N] for k = 0..N-1.
///   3. Insert expansion sequence [0,1,1,0,1,0,0] before index p-1:
///      - First p-1 chips: W[0..p-2]
///      - 7-bit expansion sequence
///      - Remaining: W[p-1..N-1]
///
///      Result: 10230 chips.
fn generate_l1cp_chips(prn: usize) -> [u8; L1C_CHIPS] {
    let (w, p) = L1CP_PARAMS[prn - 1];
    let legendre = legendre_sequence();
    let mut chips = [0u8; L1C_CHIPS];

    // Build the Weil code of length N
    let mut weil = [0u8; LEGENDRE_N];
    for k in 0..LEGENDRE_N {
        let l_k = legendre[k];
        let l_shifted = legendre[(k + w) % LEGENDRE_N];
        weil[k] = l_k ^ l_shifted;
    }

    // Insert expansion sequence before index p-1
    // p is 1-based, so we insert before the (p-1)-th element of weil
    let split_point = p - 1;

    // Copy W[0..split_point-1]
    let mut out_idx = 0;
    for &chip in weil.iter().take(split_point) {
        chips[out_idx] = chip;
        out_idx += 1;
    }

    // Insert expansion sequence
    for &bit in EXPANSION_SEQUENCE.iter() {
        chips[out_idx] = bit;
        out_idx += 1;
    }

    // Copy remaining W[split_point..N-1]
    for &chip in weil.iter().skip(split_point) {
        chips[out_idx] = chip;
        out_idx += 1;
    }

    debug_assert_eq!(out_idx, L1C_CHIPS);

    chips
}

/// Generate the ±1 L1Cp code chips for a given PRN (1-based).
///
/// Maps 0→+1, 1→-1.
pub fn generate_l1c_code(prn: usize) -> [f64; L1C_CHIPS] {
    let chips = generate_l1cp_chips(prn);
    let mut code = [0.0f64; L1C_CHIPS];
    for (i, &chip) in chips.iter().enumerate() {
        code[i] = if chip == 0 { 1.0 } else { -1.0 };
    }
    code
}

/// Generate the BOC(1,1)-modulated L1C pilot code for a given PRN (1-based).
///
/// Each chip is split into two half-chips: [+chip, -chip], giving 20460
/// samples per code period. Uses BOC(1,1) modulation for acquisition.
pub fn generate_l1c_pilot_code(prn: usize) -> Vec<f64> {
    assert!(
        (1..=L1C_NUM_SATS).contains(&prn),
        "GPS L1C PRN must be in 1..={}, got {prn}",
        L1C_NUM_SATS
    );
    let weil = generate_l1c_code(prn);

    let mut boc = Vec::with_capacity(L1C_CHIPS * 2);
    for &chip in &weil {
        boc.push(chip);
        boc.push(-chip);
    }
    boc
}

/// Resample the BOC(1,1)-modulated L1C pilot code to `samples_per_code`
/// samples per 10~ms code period, repeating for `epochs` full periods.
pub fn l1c_code_resampled(samples_per_code: f64, prn: usize, epochs: f64) -> Vec<f64> {
    let boc = generate_l1c_pilot_code(prn);
    let boc_len = boc.len() as f64; // 20460
    let samples_per_chip = samples_per_code / boc_len;
    let num_samples = (samples_per_code * epochs).floor() as usize;

    (0..num_samples)
        .map(|n| {
            let idx = ((n as f64 / samples_per_chip).floor() as usize) % boc.len();
            boc[idx]
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Testbench — cross-validated against IS-GPS-800H Table 3.2-2 reference data
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
        // exactly (N-1)/2 ones for prime N
        assert_eq!(ones, (LEGENDRE_N - 1) / 2);
    }

    // --- L1Cp chips (0/1) cross-checked against IS-GPS-800H ---------------

    /// Convert 24 chips (0/1, MSB-first, groups of 3) to an octal string.
    fn chips24_to_octal(chips: &[u8]) -> String {
        assert_eq!(chips.len(), 24);
        let mut s = String::with_capacity(8);
        for i in 0..8 {
            let d = 4 * chips[3 * i] + 2 * chips[3 * i + 1] + chips[3 * i + 2];
            s.push(char::from_digit(d as u32, 8).unwrap());
        }
        s
    }

    /// Reference octal strings for the first 24 chips of each PRN 1..=63
    /// (L1CP column) from IS-GPS-800H Table 3.2-2 "Initial 24 Chips (Octal)".
    #[rustfmt::skip]
    const REF_FIRST24: [&str; L1C_NUM_SATS] = [
        "05752067", "70146401", "32066222", "72125121", "42323273",
        "01650642", "21303446", "35504263", "66434311", "52631623",
        "04733076", "50352603", "32026612", "07476042", "22210746",
        "30706376", "75764610", "73202225", "47227426", "16064126",
        "66415734", "27600270", "66101627", "17717055", "47500232",
        "52057615", "76153566", "22444670", "62330044", "13674337",
        "60635146", "73527653", "63772350", "33564215", "52236055",
        "64506521", "73561133", "12647121", "16640265", "11161337",
        "22055260", "11546064", "24765004", "14042504", "53512265",
        "15317006", "16151224", "67454561", "47542743", "65057230",
        "77415771", "75364651", "75664330", "44600202", "23211425",
        "51504740", "47712554", "67325233", "61517015", "43217554",
        "52520062", "77073716", "56350460",
    ];

    /// Reference octal strings for the last 24 chips of each PRN 1..=63
    /// (L1CP column) from IS-GPS-800H Table 3.2-2 "Final 24 Chips (Octal)".
    #[rustfmt::skip]
    const REF_LAST24: [&str; L1C_NUM_SATS] = [
        "20173742", "35437154", "00161056", "71435437", "15035661",
        "32606570", "03475644", "11316575", "23047575", "07355246",
        "15210113", "72643606", "63457333", "46623624", "35467322",
        "70116567", "62731643", "14040613", "07750525", "37171211",
        "01302134", "37672235", "32201230", "37437553", "23310544",
        "07152415", "02571041", "52270664", "61317104", "43137330",
        "20336467", "40745656", "50272475", "75604301", "52550266",
        "15334214", "53445703", "71136024", "01607455", "73467421",
        "54372454", "11526534", "16522173", "74053703", "52211303",
        "72655147", "01212152", "10410122", "22473073", "63145220",
        "65734110", "25167435", "17524136", "47064764", "14016156",
        "11723025", "76760325", "04724615", "72504743", "51215201",
        "00630473", "71217605", "50200707",
    ];

    #[test]
    fn test_l1cp_first24_octal() {
        for prn in 1..=L1C_NUM_SATS {
            let chips = generate_l1cp_chips(prn);
            let octal = chips24_to_octal(&chips[0..24]);
            assert_eq!(
                octal,
                REF_FIRST24[prn - 1],
                "PRN {prn}: first 24 chips mismatch"
            );
        }
    }

    #[test]
    fn test_l1cp_last24_octal() {
        for prn in 1..=L1C_NUM_SATS {
            let chips = generate_l1cp_chips(prn);
            let octal = chips24_to_octal(&chips[L1C_CHIPS - 24..]);
            assert_eq!(
                octal,
                REF_LAST24[prn - 1],
                "PRN {prn}: last 24 chips mismatch"
            );
        }
    }

    // --- L1Cp bipolar (±1) code -------------------------------------------

    #[test]
    fn test_generate_l1c_code_length() {
        let code = generate_l1c_code(1);
        assert_eq!(code.len(), L1C_CHIPS);
    }

    #[test]
    fn test_generate_l1c_code_bipolar() {
        let code = generate_l1c_code(5);
        for &v in &code {
            assert!(v == 1.0 || v == -1.0, "unexpected value {v}");
        }
    }

    #[test]
    fn test_generate_l1c_pilot_code_length() {
        let code = generate_l1c_pilot_code(1);
        assert_eq!(code.len(), L1C_CHIPS * 2);
    }

    #[test]
    fn test_l1c_code_resampled_length() {
        let samples_per_code = (L1C_CHIPS * 2) as f64;
        let code = l1c_code_resampled(samples_per_code, 1, 1.0);
        assert_eq!(code.len(), L1C_CHIPS * 2);
    }

    #[test]
    fn test_all_prns_generate() {
        for prn in 1..=L1C_NUM_SATS {
            let code = generate_l1c_pilot_code(prn);
            assert_eq!(
                code.len(),
                L1C_CHIPS * 2,
                "PRN {prn} wrong pilot code length"
            );
            for &v in &code {
                assert!(v == 1.0 || v == -1.0, "PRN {prn}: unexpected value {v}");
            }
        }
    }

    #[test]
    fn test_all_params_in_range() {
        for (prn, &(w, p)) in L1CP_PARAMS.iter().enumerate() {
            assert!(
                w > 0 && w < LEGENDRE_N,
                "PRN {} pilot w={} out of range 1..{}",
                prn + 1,
                w,
                LEGENDRE_N - 1
            );
            assert!(
                p >= 1 && p <= LEGENDRE_N,
                "PRN {} pilot p={} out of range 1..={}",
                prn + 1,
                p,
                LEGENDRE_N
            );
        }
    }

    /// Cross-check: verify 0→+1, 1→-1 mapping.
    #[test]
    fn test_chip_mapping() {
        let chips = generate_l1cp_chips(1);
        let code = generate_l1c_code(1);
        for i in 0..L1C_CHIPS {
            let expected = if chips[i] == 0 { 1.0 } else { -1.0 };
            assert_eq!(code[i], expected, "chip {i}: mapping mismatch");
        }
    }
}
