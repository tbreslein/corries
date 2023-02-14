// Copyright (c) 2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports [TimeIntegrationConfig] for configuring [TimeSolver](crate::time::TimeSolver) objects.

use crate::errorhandling::Validation;
use color_eyre::{
    eyre::{ensure, Context},
    Result,
};
use serde::Serialize;

/// Enumerates the different types of configuration [TimeSolver](crate::time::TimeSolver) objects.
#[derive(Debug, Serialize, Copy, Clone)]
pub enum TimeIntegrationConfig {
    /// Configuration for the Runge-Kutta-Fehlberg solver, i.e.
    /// [RungeKuttaFehlberg](crate::time::rkf::RungeKuttaFehlberg).
    ///
    /// The payload for this variant is a [RkfConfig] object.
    Rkf(RkfConfig),
}

unsafe impl Send for TimeIntegrationConfig {}
unsafe impl Sync for TimeIntegrationConfig {}

impl TimeIntegrationConfig {
    /// Constructs a [TimeIntegrationConfig] for the default Runge-Kutta-Fehlberg setup.
    ///
    /// Check the asserts in the example for defaults.
    ///
    /// ```
    /// use corries::prelude::*;
    ///
    /// // define the config instance
    /// let timeintegration_config = TimeIntegrationConfig::default_rkf();
    ///
    /// if let TimeIntegrationConfig::Rkf(c) = timeintegration_config {
    ///     assert_eq!(c.rkf_mode, RKFMode::SSPRK5);
    ///     assert_eq!(c.asc, false);
    ///     assert_eq!(c.asc_relative_tolerance, 0.001);
    ///     assert_eq!(c.asc_absolute_tolerance, 0.001);
    ///     assert_eq!(c.asc_timestep_friction, 0.08);
    /// }
    /// ```
    pub fn default_rkf() -> Self {
        Self::Rkf(RkfConfig::default())
    }
}

impl Validation for TimeIntegrationConfig {
    fn validate(&self) -> Result<()> {
        match self {
            TimeIntegrationConfig::Rkf(c) => c.validate().context("Validating RkfConfig"),
        }
    }
}

/// Enum for the different kinds of Runge-Kutta-Fehlberg schemes available
///
/// Defaults to [SSPRK5](RKFMode::SSPRK5)
#[derive(Debug, Serialize, Copy, Clone, Default, PartialEq, Eq)]
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
    #[default]
    SSPRK5,
}

unsafe impl Send for RKFMode {}
unsafe impl Sync for RKFMode {}

/// Configures the [RungeKuttaFehlberg](crate::time::rkf::RungeKuttaFehlberg) objects.
#[derive(Debug, Serialize, Copy, Clone)]
pub struct RkfConfig {
    /// Type of Runge-Kutta-Fehlberg scheme to use
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

unsafe impl Send for RkfConfig {}
unsafe impl Sync for RkfConfig {}

impl Default for RkfConfig {
    fn default() -> Self {
        Self {
            rkf_mode: RKFMode::default(),
            asc: false,
            asc_relative_tolerance: 0.001,
            asc_absolute_tolerance: 0.001,
            asc_timestep_friction: 0.08,
        }
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
