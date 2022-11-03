// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [TimeIntegration] struct, and the [DtKind] enum.

use color_eyre::{eyre::bail, Result};

use crate::{
    config::{numericsconfig::TimeIntegrationConfig, CorriesConfig},
    errorhandling::Validation,
    mesh::Mesh,
    physics::Physics,
    rhs::Rhs,
};

use self::{rkf::RungeKuttaFehlberg, timestep::TimeStep};

mod rkf;
mod timestep;

/// Enumerates the different kinds of effects that can limit the time step width.
pub enum DtKind {
    /// The initial value of such a variable at the beginning of the simulation
    Init,

    /// CFL limited
    Cfl,
}

/// Trait for objects that solve the time integration step and produce new solutions
trait TimeSolver<const S: usize, const EQ: usize> {
    /// Calculates the next solution for the physical state
    ///
    /// # Arguments
    ///
    /// * `time` - Information about the time coordinate
    /// * `u` - The current physical state
    /// * `rhs` - Solves the right-hand side
    /// * `mesh` - Information about spatial properties
    fn next_solution(
        &mut self,
        time: &mut TimeStep,
        u: &mut Physics<S, EQ>,
        rhs: &mut Rhs<S, EQ>,
        mesh: &Mesh<S>,
    ) -> Result<()>;
}

/// Struct for everything regarding time step integration.
///
/// Carries a [TimeStep] for keeping track of the time coordinate and data regarding it, as well as
/// a `Box<dyn TimeSolver>` for calculating new solutions.
pub struct TimeIntegration<const S: usize, const EQ: usize> {
    /// Keeps track of the time coordinate and related data
    pub time: TimeStep,

    /// Calculates solutions for the [Physics] state
    solver: Box<dyn TimeSolver<S, EQ>>,
}

impl<const S: usize, const EQ: usize> TimeIntegration<S, EQ> {
    /// Constructs a new [TimeIntegration] struct.
    ///
    /// # Arguments
    ///
    /// * `config` - a [CorriesConfig] configuration object
    pub fn new(config: &CorriesConfig) -> Self {
        return Self {
            time: TimeStep::new(&config.numericsconfig, config.output_counter_max),
            solver: Box::new(match &config.numericsconfig.time_integration_config {
                TimeIntegrationConfig::Rkf(rkfconfig) => RungeKuttaFehlberg::new(rkfconfig, &config.physicsconfig),
            }),
        };
    }

    /// Calculates the next state for the [Physics] object `u`.
    ///
    /// # Arguments
    ///
    /// * `u` - the [Physics] state being modified to transition between current and next state
    /// * `rhs` - solves the right-hand side
    /// * `mesh` - Information about spatial properties
    pub fn next_solution(&mut self, u: &mut Physics<S, EQ>, rhs: &mut Rhs<S, EQ>, mesh: &Mesh<S>) -> Result<()> {
        self.solver.next_solution(&mut self.time, u, rhs, mesh)?;
        if self.time.iter >= self.time.iter_max {
            bail!(
                "time.iter reached time.iter_max! time.iter = {}, time.iter_max = {}",
                self.time.iter,
                self.time.iter_max
            );
        }
        return Ok(());
    }
}

impl<const S: usize, const EQ: usize> Validation for TimeIntegration<S, EQ> {
    fn validate(&self) -> Result<()> {
        todo!();
    }
}
