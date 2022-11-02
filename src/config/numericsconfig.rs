// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports [NumericsConfig] for configuring the numerical schemes

use color_eyre::Result;

use crate::errorhandling::Validation;

/// Enum for the different kinds of numerical flux schemes available
#[derive(Debug, Copy, Clone)]
pub enum NumFluxMode {
    /// HLL scheme
    Hll,
}

/// Carries information about how the mesh should shaped
#[derive(Debug)]
pub struct NumericsConfig {
    /// The type of numerical flux scheme to use
    pub numflux_mode: NumFluxMode,
}

impl Validation for NumericsConfig {
    fn validate(&self) -> Result<()> {
        return Ok(());
    }
}
