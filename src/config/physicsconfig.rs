// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports [PhysicsConfig] for configuring the physics equations and conditions.

use color_eyre::{eyre::ensure, Result};
use serde::Serialize;
use crate::{errorhandling::Validation, UnitsMode};

/// Carries information about how the mesh should shaped
#[derive(Debug, Serialize, Clone, Default)]
pub struct PhysicsConfig {
    /// Ratio of specific heats
    pub adiabatic_index: f64,

    /// The units system
    pub units_mode: UnitsMode,
}

unsafe impl Send for PhysicsConfig {}
unsafe impl Sync for PhysicsConfig {}

impl Validation for PhysicsConfig {
    fn validate(&self) -> Result<()> {
        ensure!(
            self.adiabatic_index > 1.0 && self.adiabatic_index < 2.0,
            "This must hold: 1 < adiabatic_index < 2 ! Got {}",
            self.adiabatic_index
        );
        Ok(())
    }
}
