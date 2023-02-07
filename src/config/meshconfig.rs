// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports [MeshConfig] for configuring the mesh of the simulation.

use color_eyre::{eyre::ensure, Result};
use serde::Serialize;
use crate::errorhandling::Validation;

/// Enum for the different kinds of Meshes available; defaults to cartesian mesh
#[derive(Debug, Serialize, Copy, Clone, Default, PartialEq, Eq)]
pub enum MeshMode {
    /// Cartesian mesh (default)
    #[default]
    Cartesian,
}

unsafe impl Send for MeshMode {}
unsafe impl Sync for MeshMode {}

/// Carries information about how the mesh should shaped
#[derive(Debug, Serialize, Clone, Default)]
pub struct MeshConfig {
    /// The type of Mesh that should constructed
    pub mode: MeshMode,

    /// xi coordinate at the inner/western boundary of the computational area
    pub xi_in: f64,

    /// xi coordinate at the outer/eastern boundary of the computational area
    pub xi_out: f64,

    /// Ratio in terms of radial thickness between the massive disk and the computational area
    pub ratio_disk: f64,
}

unsafe impl Send for MeshConfig {}
unsafe impl Sync for MeshConfig {}

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
            "This must hold: 0.0 < ratio_disk <= 1.0! Got {}",
            self.ratio_disk
        );
        Ok(())
    }
}
