// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::eyre::ensure;
use color_eyre::Result;

use crate::errorhandling::Validation;

/// Enum for the different kinds of Physics available
#[derive(Debug, Copy, Clone)]
pub enum PhysicsMode {
    Euler1DIsot,
}

/// Enum for the different kinds of unit systems available
#[derive(Debug, Copy, Clone)]
pub enum UnitsMode {
    SI,
}

/// Carries information about how the mesh should shaped
#[derive(Debug)]
pub struct PhysicsConfig {
    /// The type of Mesh that should constructed
    pub mode: PhysicsMode,

    /// The number of cells in the computational area
    pub units_mode: UnitsMode,

    /// Mean molecular mass in SI units
    pub mu_si: f64,

    /// Critical Reynolds number
    pub critical_reynolds: f64,

    /// Ratio of specific heats
    pub adiabatic_index: f64,

    /// Global initial speed of sound
    pub c_sound_0: f64,
    
    /// Global initial temperature
    pub temperature_0: f64,

    /// Whether to init pressure and temperature using c_sound_0 as a starting point; set to false
    /// to start with temperature_0 instead.
    pub init_with_c_sound: bool,

    /// Starting mass for the central object
    pub m_central_0: f64,

    /// Starting mass for the disk
    pub m_disk_0: f64,

    /// Accretion efficiency
    pub accretion_efficiency: f64,

    /// Whether to use the eddington limit to limit the growth rate of the central object
    pub use_eddington_limit: bool,
}

impl Validation for PhysicsConfig {
    fn validate(&self) -> Result<()> {
        ensure!(
            self.mu_si > 0.0,
            "This must hold: mu_si > 0! Got mu_si = {}",
            self.mu_si
        );
        ensure!(
            self.critical_reynolds > 0.0,
            "This must hold: critical_reynolds > 0! Got critical_reynolds = {}",
            self.critical_reynolds
        );
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
        ensure!(
            self.temperature_0 > 0.0,
            "This must hold: temperature_0 > 0! Got temperature_0 = {}",
            self.temperature_0
        );
        ensure!(
            self.m_central_0 > 0.0,
            "This must hold: m_central_0 > 0! Got m_central_0 = {}",
            self.m_central_0
        );
        ensure!(
            self.m_disk_0 > 0.0,
            "This must hold: m_disk_0 > 0! Got m_disk_0 = {}",
            self.m_disk_0
        );
        ensure!(
            self.accretion_efficiency > 0.0 && self.accretion_efficiency <= 1.0,
            "This must hold: 0 < accretion_efficiency <= 1! Got accretion_efficiency = {}",
            self.accretion_efficiency
        );
        return Ok(());
    }
}

