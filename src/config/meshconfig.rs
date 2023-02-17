// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports [MeshConfig] for configuring the [Mesh](crate::mesh::Mesh) of the simulation.

use crate::errorhandling::Validation;
use color_eyre::{eyre::ensure, Result};
use serde::Serialize;

/// Enum for the different kinds of Meshes available
#[derive(Debug, Serialize, Copy, Clone, PartialEq, Eq)]
pub enum MeshMode {
    /// Cartesian mesh
    Cartesian,
}

unsafe impl Send for MeshMode {}
unsafe impl Sync for MeshMode {}

/// Carries information about how the mesh should shaped
#[derive(Debug, Serialize, Clone)]
pub struct MeshConfig {
    /// The type of Mesh that should constructed
    pub mode: MeshMode,

    /// xi coordinate at the inner/western boundary of the computational area
    pub xi_in: f64,

    /// xi coordinate at the outer/eastern boundary of the computational area
    pub xi_out: f64,
}

impl MeshConfig {
    /// Sets up a MeshConfig that can be used for most Riemann tests, which uses a Cartesian mesh
    /// with coordinate boundaries set to 1.0 and 2.0.
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    ///
    /// let meshconfig = MeshConfig::default_riemann_test();
    /// assert_eq!(meshconfig.mode, MeshMode::Cartesian);
    /// assert_eq!(meshconfig.xi_in, 1.0);
    /// assert_eq!(meshconfig.xi_out, 2.0);
    /// ```
    pub fn default_riemann_test() -> Self {
        Self {
            mode: MeshMode::Cartesian,
            xi_in: 1.0,
            xi_out: 2.0,
        }
    }
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
        Ok(())
    }
}
