// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports [PhysicsConfig] for configuring the physics equations and conditions.

use color_eyre::eyre::ensure;
use color_eyre::Result;
use serde::Serialize;

use crate::errorhandling::Validation;

/// Carries information about how the mesh should shaped
#[derive(Debug, Serialize)]
pub struct PhysicsConfig {
    /// Ratio of specific heats
    pub adiabatic_index: f64,
}

impl Validation for PhysicsConfig {
    fn validate(&self) -> Result<()> {
        ensure!(
            self.adiabatic_index > 1.0 && self.adiabatic_index < 2.0,
            "This must hold: 1 < adiabatic_index < 2 ! Got {}",
            self.adiabatic_index
        );
        return Ok(());
    }
}
