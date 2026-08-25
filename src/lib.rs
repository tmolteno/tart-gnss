// Copyright (c) 2026 Tim Molteno <tim@elec.ac.nz>
// SPDX-License-Identifier: GPL-3.0

//! GNSS signal acquisition for the TART radio telescope — GPS L1 C/A,
//! GPS L1C, Galileo E1-C, BeiDou B1C, SBAS L1 C/A and QZSS L1 C/A.
//!
//! Library crate for the `tart-gnss-acquire` binary.  See the crate
//! README for command-line usage.

pub mod acr;
pub mod acquisition;
pub mod antenna_test;
pub mod beidou;
pub mod config;
pub mod correlate;
pub mod galileo;
pub mod l1c;
pub mod observation;
pub mod qzss;
pub mod sbas;
pub mod simulator;
pub mod stats;

#[cfg(test)]
mod testutil;

use serde::Serialize;

/// Combined acquisition output across all selected constellations.
#[derive(Serialize)]
pub struct CombinedOutput {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gps: Option<acquisition::GpsAllAcquisitionOutput>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub galileo: Option<galileo::GalileoAllAcquisitionOutput>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub beidou: Option<beidou::BeiDouAllAcquisitionOutput>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sbas: Option<sbas::SbasAllAcquisitionOutput>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub l1c: Option<l1c::L1CAllAcquisitionOutput>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub qzss: Option<qzss::QzssAllAcquisitionOutput>,
}
