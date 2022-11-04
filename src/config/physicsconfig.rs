// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports [PhysicsConfig] for configuring the physics equations and conditions.

use color_eyre::eyre::ensure;
use color_eyre::Result;

use crate::errorhandling::Validation;
use crate::units::UnitsMode;

/// Enum for the different kinds of Physics available
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum PhysicsMode {
    /// 1D isothermal Euler equations
    Euler1DIsot,

    /// 2D isothermal Euler equations
    Euler2DIsot,

    /// 1D adiabatic Euler equations
    Euler1DAdiabatic,
    // Euler2DAdiabatic,
    // NavStoIsot,
    // NavStoAdiabatic,
}

/// Carries information about how the mesh should shaped
#[derive(Debug)]
pub struct PhysicsConfig {
    /// The type of Mesh that should constructed
    pub mode: PhysicsMode,

    /// The number of cells in the computational area
    pub units_mode: UnitsMode,

    /// Ratio of specific heats
    pub adiabatic_index: f64,

    /// Global initial speed of sound
    pub c_sound_0: f64,
}

impl Validation for PhysicsConfig {
    fn validate(&self) -> Result<()> {
        ensure!(
            self.adiabatic_index > 1.0 && self.adiabatic_index < 2.0,
            "This must hold: 1 < adiabatic_index < 2 ! Got {}",
            self.adiabatic_index
        );
        ensure!(
            self.c_sound_0 > 0.0,
            "This must hold: c_sound_0 > 0! Got {}",
            self.c_sound_0
        );
        return Ok(());
    }
}
