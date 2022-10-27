// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::eyre::ensure;
use color_eyre::Result;

use crate::errorhandling::Validation;

/// Enum for the different kinds of Meshes available
#[derive(Debug, Copy, Clone)]
pub enum MeshMode {
    Cartesian,
}

/// Carries information about how the mesh should shaped
#[derive(Debug)]
pub struct MeshConfig {
    /// The type of Mesh that should constructed
    pub mode: MeshMode,

    /// The number of cells in the computational area
    pub n_comp: usize,

    /// The number of ghost cells per edge
    pub n_gc: usize,

    /// xi coordinate at the inner/western boundary of the computational area
    pub xi_in: f64,

    /// xi coordinate at the outer/eastern boundary of the computational area
    pub xi_out: f64,

    /// Ratio in terms of radial thickness between the massive disk and the computational area
    pub ratio_disk: f64,
}

impl Validation for MeshConfig {
    fn validate(&self) -> Result<()> {
        ensure!(
            self.xi_in < self.xi_out,
            "This must hold: xi_in < xi_out! Got xi_in = {} ; xi_out = {}",
            self.xi_in,
            self.xi_out
        );
        ensure!(
            self.ratio_disk > 0.0 && self.ratio_disk <= 1.0,
            "This must hold: 0.0 < ratio_disk <= 1.0! Got ratio_disk = {}",
            self.ratio_disk
        );
        return Ok(());
    }
}
