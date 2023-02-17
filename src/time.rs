// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [TimeSolver] trait and the [DtKind] enum.

use crate::{config::CorriesConfig, mesh::Mesh, rhs::Rhs, state::Physics, NumFlux, State};
use color_eyre::Result;
use std::fmt::Display;

pub mod rkf;
pub mod timestep;
pub use self::rkf::RungeKuttaFehlberg;
use self::timestep::TimeStep;

/// Enumerates the different kinds of effects that can limit the time step width.
#[derive(Debug, Copy, Clone, Default, PartialEq, Eq)]
pub enum DtKind {
    /// The initial value of such a variable at the beginning of the simulation.
    #[default]
    Init,

    /// Denotes that the time step was limited by the CFL criterium.
    Cfl,

    /// Used when dumping state because of an error.
    ErrorDump,
}

unsafe impl Send for DtKind {}
unsafe impl Sync for DtKind {}

impl Display for DtKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DtKind::Init => write!(f, "init"),
            DtKind::Cfl => write!(f, "cfl"),
            DtKind::ErrorDump => write!(f, "err"),
        }
    }
}

/// Trait for objects that solve the time integration step and produce new solutions
pub trait TimeSolver<P: Physics<E, S>, const E: usize, const S: usize> {
    /// Constructs a new [TimeSolver] object
    ///
    /// # Arguments
    ///
    /// * `config` - Configuration for the [corries](crate) simulation
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    ///
    /// const S: usize = 100;
    /// set_Physics_and_E!(Euler1DIsot);
    /// type N = Hll<E,S>;
    /// type T = RungeKuttaFehlberg<P, E, S>;
    ///
    /// // use the default config for Riemann tests
    /// let t_end = 0.5;
    /// let folder_name = "results";
    /// let file_name = "noh";
    ///
    /// // define the config instance
    /// let config = CorriesConfig::default_riemann_test::<N, E, S>(t_end, folder_name, file_name);
    ///
    /// let time = T::new(&config).unwrap();
    /// ```
    fn new(config: &CorriesConfig) -> Result<Self>
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
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use corries::prelude::*;
    ///
    /// // Set up a simulation using isothermal Euler physics on a 100 cell mesh, using the Hll and
    /// // Runge-Kutta-Fehlberg schemes for solving the equations.
    /// const S: usize = 100;
    /// set_Physics_and_E!(Euler1DIsot);
    /// type N = Hll<E, S>;
    /// type T = RungeKuttaFehlberg<P, E, S>;
    ///
    /// // use the default config for Riemann tests
    /// let t_end = 0.5;
    /// let folder_name = "results";
    /// let file_name = "noh";
    ///
    /// let config = CorriesConfig::default_riemann_test::<N, E, S>(t_end, folder_name, file_name);
    /// let (mut u, mut solver, mesh, _) = config.init_corries::<P, N, T, E, S>(|_,_,_| Ok(())).unwrap();
    ///
    /// // solver carries the an instance of [TimeSolver]
    /// solver.next_solution(&mut u, &mesh).unwrap();
    /// ```
    fn next_solution<N: NumFlux<E, S>>(
        &mut self,
        time: &mut TimeStep,
        u: &mut State<P, E, S>,
        rhs: &mut Rhs<N, E, S>,
        mesh: &Mesh<S>,
    ) -> Result<()>;
}
