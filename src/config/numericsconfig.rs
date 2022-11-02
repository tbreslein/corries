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

/// Enum for the different kinds of time integration schemes available
#[derive(Debug, Copy, Clone)]
pub enum TimeIntegrationMode {
    /// Runge-Kutta-Fehlberg
    RKF,
}

/// Carries information about how the mesh should shaped
#[derive(Debug)]
pub struct NumericsConfig {
    /// The type of numerical flux scheme to use
    pub numflux_mode: NumFluxMode,

    /// The type of time integration scheme to use
    pub time_integration_mode: TimeIntegrationMode,

    /// Maximum number of time step iterations throughout the simulation
    pub iter_max: usize,

    /// Time coordinate at the start of the simulation
    pub t0: f64,

    /// Time coordinate at which the simulation should stop
    pub t_end: f64,

    /// Minimal value for dt at which the simulation errors out if the time step dips below it
    pub dt_min: f64,

    /// Upper limit for the dt
    pub dt_max: f64,

    /// CFL parameter
    pub dt_cfl_param: f64,
}

impl Validation for NumericsConfig {
    fn validate(&self) -> Result<()> {
        return Ok(());
    }
}
