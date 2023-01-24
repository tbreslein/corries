// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [TimeSolver] trait, the [Time] struct, and the [DtKind] enum.

use std::fmt::Display;
use std::marker::PhantomData;

use color_eyre::{
    eyre::{bail, Context},
    Result,
};

use crate::{config::CorriesConfig, mesh::Mesh, rhs::Rhs, state::Physics, NumFlux};

pub mod rkf;
pub mod timestep;
pub use self::rkf::RungeKuttaFehlberg;
use self::timestep::TimeStep;

/// Enumerates the different kinds of effects that can limit the time step width.
pub enum DtKind {
    /// The initial value of such a variable at the beginning of the simulation
    Init,

    /// CFL limited
    Cfl,

    /// Used when dumping state because of an error
    ErrorDump,
}

impl Display for DtKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        return match self {
            DtKind::Init => write!(f, "init"),
            DtKind::Cfl => write!(f, "cfl"),
            DtKind::ErrorDump => write!(f, "err"),
        };
    }
}

/// Trait for objects that solve the time integration step and produce new solutions
pub trait TimeSolver<P: Physics<E, S>, const E: usize, const S: usize> {
    /// Constructs a new [TimeSolver] object
    ///
    /// # Arguments
    ///
    /// * `rkfconfig` - Configuration specifically for [RungeKuttaFehlberg] objects
    /// * `physicsconfig` - Configuration for [Physics] objects, needed because `utilde`
    fn new(config: &CorriesConfig, u: &P) -> Result<Self>
    where
        Self: Sized;

    /// Calculates the next solution for the physical state
    ///
    /// # Arguments
    ///
    /// * `time` - Information about the time coordinate
    /// * `u` - The current physical state
    /// * `rhs` - Solves the right-hand side
    /// * `mesh` - Information about spatial properties
    fn next_solution<N: NumFlux<E, S>>(
        &mut self,
        time: &mut TimeStep,
        u: &mut P,
        rhs: &mut Rhs<N, E, S>,
        mesh: &Mesh<S>,
    ) -> Result<()>;
}

/// Struct for everything regarding time step integration.
///
/// Carries a [TimeStep] for keeping track of the time coordinate and data regarding it, as well as
/// a `Box<dyn TimeSolver>` for calculating new solutions.
pub struct Time<P: Physics<E, S>, T: TimeSolver<P, E, S>, const E: usize, const S: usize> {
    /// Keeps track of the time coordinate and related data
    pub timestep: TimeStep,

    /// Calculates solutions for the [Physics] state
    solver: T,

    embedded_type: PhantomData<P>,
}

impl<P: Physics<E, S>, T: TimeSolver<P, E, S>, const E: usize, const S: usize> Time<P, T, E, S> {
    /// Constructs a new [Time] struct.
    ///
    /// # Arguments
    ///
    /// * `config` - a [CorriesConfig] configuration object
    pub fn new(config: &CorriesConfig, u: &P) -> Result<Self> {
        let solver = T::new(config, u)?;
        return Ok(Self {
            timestep: TimeStep::new(&config.numerics_config, config.output_counter_max),
            solver,
            embedded_type: PhantomData,
        });
    }

    /// Calculates the next state for the [Physics] object `u`.
    ///
    /// # Arguments
    ///
    /// * `u` - the [Physics] state being modified to transition between current and next state
    /// * `rhs` - solves the right-hand side
    /// * `mesh` - Information about spatial properties
    pub fn next_solution<N: NumFlux<E, S>>(&mut self, u: &mut P, rhs: &mut Rhs<N, E, S>, mesh: &Mesh<S>) -> Result<()> {
        self.solver
            .next_solution::<N>(&mut self.timestep, u, rhs, mesh)
            .context("Calling TimeIntegration::solver.next_solution in TimeIntegration::next_solution")?;
        if self.timestep.iter >= self.timestep.iter_max {
            bail!(
                "time.iter reached time.iter_max! time.iter = {}, time.iter_max = {}",
                self.timestep.iter,
                self.timestep.iter_max
            );
        }
        return Ok(());
    }
}
