// Copyright (c) 2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports [TimeIntegrationConfig] for configuring [TimeSolver](crate::time::TimeSolver) objects.

use crate::{check_positive_double, errorhandling::Validation};
use color_eyre::{eyre::ensure, Result};
use serde::Serialize;

/// Enumerates the different types of configuration [TimeSolver](crate::time::TimeSolver) objects.
#[derive(Debug, Serialize, Copy, Clone)]
pub enum TimeIntegrationConfig {
    /// Configuration for the Runge-Kutta-Fehlberg solver, i.e.
    /// [RungeKuttaFehlberg](crate::time::rkf::RungeKuttaFehlberg).
    ///
    /// This variant disables automated step control
    Rkf {
        /// Type of RungeKuttaFehlberg scheme to use
        rkf_mode: RKFMode,
    },

    /// Configuration for the Runge-Kutta-Fehlberg solver, i.e.
    /// [RungeKuttaFehlberg](crate::time::rkf::RungeKuttaFehlberg).
    ///
    /// This variant includes automated step control
    RkfASC {
        /// Type of RungeKuttaFehlberg scheme to use
        rkf_mode: RKFMode,

        /// Relative tolerance for automated time step control
        relative_tolerance: f64,

        /// Absolute tolerance for automated time step control
        absolute_tolerance: f64,

        /// Timestep "friction" for automated time step control
        timestep_friction: f64,
    },
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
    /// if let TimeIntegrationConfig::Rkf { rkf_mode } = timeintegration_config {
    ///     assert_eq!(rkf_mode, RKFMode::SSPRK5);
    /// }
    /// ```
    pub fn default_rkf() -> Self {
        Self::Rkf {
            rkf_mode: RKFMode::default(),
        }
    }
}

impl Validation for TimeIntegrationConfig {
    fn validate(&self) -> Result<()> {
        match self {
            TimeIntegrationConfig::Rkf { rkf_mode: _ } => Ok(()),
            TimeIntegrationConfig::RkfASC {
                rkf_mode: _,
                relative_tolerance,
                absolute_tolerance,
                timestep_friction,
            } => {
                check_positive_double!(*absolute_tolerance, *relative_tolerance, *timestep_friction);
                Ok(())
            },
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
