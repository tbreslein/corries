// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports [NumericsConfig] for configuring [NumFlux] and [Time](crate::time::Time) objects

use std::any::TypeId;

use crate::{errorhandling::Validation, Hll, Kt, NumFlux};
use color_eyre::{
    eyre::{ensure, Context},
    Result,
};
pub use numfluxconfig::*;
use serde::Serialize;
pub use timeintegrationconfig::*;

mod numfluxconfig;
mod timeintegrationconfig;

/// Carries information about how the mesh should shaped
#[derive(Debug, Serialize, Clone)]
pub struct NumericsConfig {
    /// Configures [NumFlux] objects
    pub numflux_config: NumFluxConfig,

    /// Configures [TimeSolver](crate::time::TimeSolver)
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

impl NumericsConfig {
    /// Sets up a [NumericsConfig] for most Riemann tests.
    ///
    /// It needs to be passed the type of [NumFlux] for the simulation as a template parameter so
    /// that the `numflux_config` field can be set accordingly.
    /// It also sets up the default Runge-Kutta-Fehlberg setup.
    ///
    /// You can check the other defaults that are being set by checking the asserts in the example.
    ///
    /// # Arguments
    ///
    /// * `t_end` - The time coordinate at which to terminate the simulation
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    ///
    /// // set up constants
    /// set_Physics_and_E!(Euler1DAdiabatic);
    /// const S: usize = 100;
    /// type N = Kt<E,S>;
    ///
    /// let t_end = 1.0;
    ///
    /// // define the config instance
    /// let numerics_config = NumericsConfig::default_riemann_test::<N, E, S>(t_end);
    ///
    /// assert_eq!(numerics_config.numflux_config, NumFluxConfig::Kt { limiter_mode: LimiterMode::VanLeer });
    /// assert_eq!(numerics_config.iter_max, usize::MAX - 2);
    /// assert_eq!(numerics_config.t0, 0.0);
    /// assert_eq!(numerics_config.t_end, t_end);
    /// assert_eq!(numerics_config.dt_min, 1.0e-12);
    /// assert_eq!(numerics_config.dt_max, f64::MAX);
    /// assert_eq!(numerics_config.dt_cfl_param, 0.4);
    /// ```
    pub fn default_riemann_test<N: NumFlux<E, S> + 'static, const E: usize, const S: usize>(t_end: f64) -> Self {
        Self {
            numflux_config: if TypeId::of::<N>() == TypeId::of::<Hll<E, S>>() {
                NumFluxConfig::Hll
            } else if TypeId::of::<N>() == TypeId::of::<Kt<E, S>>() {
                NumFluxConfig::Kt {
                    limiter_mode: LimiterMode::VanLeer,
                }
            } else {
                panic!("Tried constructing NumericsConfig::default_riemann_test, but the case for N was not covered!")
            },
            time_integration_config: TimeIntegrationConfig::default_rkf(),
            iter_max: usize::MAX - 2,
            t0: 0.0,
            t_end,
            dt_min: 1.0e-12,
            dt_max: f64::MAX,
            dt_cfl_param: 0.4,
        }
    }
}

unsafe impl Send for NumericsConfig {}
unsafe impl Sync for NumericsConfig {}

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
