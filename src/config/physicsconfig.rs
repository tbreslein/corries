// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports [PhysicsConfig] for configuring [Physics](crate::state::Physics) objects.

use crate::{errorhandling::Validation, UnitsMode};
use color_eyre::{eyre::ensure, Result};
use serde::Serialize;

/// Carries information about how to construct [Physics](crate::state::Physics) objects.
///
/// Defaults to setting:
///
/// * `adiabatic_index` to 5.0 / 3.0
/// * `units_mode` to [UnitsMode::default()]
#[derive(Debug, Serialize, Clone)]
pub struct PhysicsConfig {
    /// Ratio of specific heats
    pub adiabatic_index: f64,

    /// The units system
    pub units_mode: UnitsMode,
}

impl Default for PhysicsConfig {
    fn default() -> Self {
        Self {
            adiabatic_index: 5.0 / 3.0,
            units_mode: UnitsMode::default(),
        }
    }
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
