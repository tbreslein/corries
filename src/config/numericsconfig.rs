// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports [NumericsConfig] for configuring the numerical schemes

use color_eyre::{
    eyre::{ensure, Context},
    Result,
};
use serde::Serialize;

mod numfluxconfig;
pub use numfluxconfig::*;

use crate::errorhandling::Validation;

/// Enum for the different kinds of Runge-Kutta-Fehlberg schemes available
#[derive(Debug, Serialize, Copy, Clone)]
pub enum RKFMode {
    /// Runge-Kutta 1
    RK1,

    /// Runge-Kutta 2
    RK2,

    /// Runge-Kutta 3
    RK3,

    /// Runge-Kutta 4
    RK4,

    /// Second-order Heun
    Heun2,

    /// Runge-Kutta-Fehlberg 1(2)
    RKF12,

    /// Runge-Kutta-Fehlberg 4(5)
    RKF45,

    /// Strong stability preserving Runge-Kutta 3
    SSPRK3,

    /// Strong stability preserving Runge-Kutta 5
    SSPRK5,
}

/// Enumerates the different types of configuration for time integration schemes
///
/// TODO: Export this into its own module
#[derive(Debug, Serialize)]
pub enum TimeIntegrationConfig {
    /// Runge-Kutta-Fehlberg
    Rkf(RkfConfig),
}

impl Validation for TimeIntegrationConfig {
    fn validate(&self) -> Result<()> {
        match self {
            TimeIntegrationConfig::Rkf(c) => c.validate().context("Validating RkfConfig"),
        }
    }
}

/// Configures the Runge-Kutta-Fehlberg scheme
#[derive(Debug, Serialize)]
pub struct RkfConfig {
    /// Type of RKF scheme to use
    pub rkf_mode: RKFMode,

    /// Whether to use automated time step control
    pub asc: bool,

    /// Relative tolerance for automated time step control
    pub asc_relative_tolerance: f64,

    /// Absolte tolerance for automated time step control
    pub asc_absolute_tolerance: f64,

    /// Timestep "friction" for automated time step control
    pub asc_timestep_friction: f64,
}

/// Carries information about how the mesh should shaped
#[derive(Debug, Serialize)]
pub struct NumericsConfig {
    /// Configures the numerical flux schemes
    pub numflux_config: NumFluxConfig,

    /// The type of time integration scheme to use
    pub time_integration_config: TimeIntegrationConfig,

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
        self.numflux_config.validate()?;
        ensure!(
            self.iter_max > 0,
            "This must hold: iter_max > 0 ! Got {}",
            self.iter_max
        );
        ensure!(self.t0 >= 0.0, "This must hold: t0 > 0.0 ! Got {}", self.t0);
        ensure!(
            self.t_end > self.t0,
            "This must hold: t_end > t0 ! Got t_end = {} ; t0 = {}",
            self.t_end,
            self.t0
        );
        ensure!(self.dt_min > 0.0, "This must hold: dt_min > 0.0 ! Got {}", self.dt_min);
        ensure!(self.dt_max > 0.0, "This must hold: dt_max > 0.0 ! Got {}", self.dt_max);
        ensure!(
            self.dt_cfl_param > 0.0,
            "This must hold: dt_cfl_param > 0.0 ! Got {}",
            self.dt_cfl_param
        );
        self.time_integration_config
            .validate()
            .context("Validating config.numericsconfig.time_integration_config")?;
        Ok(())
    }
}

impl Validation for RkfConfig {
    fn validate(&self) -> Result<()> {
        ensure!(
            self.asc_absolute_tolerance > 0.0,
            "This must hold: asc_absolute_tolerance > 0.0 ! Got {}",
            self.asc_absolute_tolerance
        );
        ensure!(
            self.asc_relative_tolerance > 0.0,
            "This must hold: asc_relative_tolerance > 0.0 ! Got {}",
            self.asc_relative_tolerance
        );
        ensure!(
            self.asc_timestep_friction > 0.0,
            "This must hold: asc_timestep_friction > 0.0 ! Got {}",
            self.asc_timestep_friction
        );
        Ok(())
    }
}
