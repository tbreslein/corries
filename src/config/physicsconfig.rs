// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports [PhysicsConfig] for configuring the physics equations and conditions.

use color_eyre::eyre::ensure;
use color_eyre::Result;

use crate::errorhandling::Validation;
use crate::units::UnitsMode;

/// Enum for the different kinds of Physics available
#[derive(Debug, Copy, Clone)]
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
        // ensure!(
        //     self.mu_si > 0.0,
        //     "This must hold: mu_si > 0! Got mu_si = {}",
        //     self.mu_si
        // );
        // ensure!(
        //     self.critical_reynolds > 0.0,
        //     "This must hold: critical_reynolds > 0! Got critical_reynolds = {}",
        //     self.critical_reynolds
        // );
        ensure!(
            self.adiabatic_index > 0.0 && self.adiabatic_index < 1.0,
            "This must hold: 0 < adiabatic_index < 1! Got adiabatic_index = {}",
            self.adiabatic_index
        );
        ensure!(
            self.c_sound_0 > 0.0,
            "This must hold: c_sound_0 > 0! Got c_sound_0 = {}",
            self.c_sound_0
        );
        // ensure!(
        //     self.temperature_0 > 0.0,
        //     "This must hold: temperature_0 > 0! Got temperature_0 = {}",
        //     self.temperature_0
        // );
        // ensure!(
        //     self.m_central_0 > 0.0,
        //     "This must hold: m_central_0 > 0! Got m_central_0 = {}",
        //     self.m_central_0
        // );
        // ensure!(
        //     self.m_disk_0 > 0.0,
        //     "This must hold: m_disk_0 > 0! Got m_disk_0 = {}",
        //     self.m_disk_0
        // );
        // ensure!(
        //     self.accretion_efficiency > 0.0 && self.accretion_efficiency <= 1.0,
        //     "This must hold: 0 < accretion_efficiency <= 1! Got accretion_efficiency = {}",
        //     self.accretion_efficiency
        // );
        return Ok(());
    }
}
