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

use crate::{config::CorriesConfig, mesh::Mesh, rhs::Rhs, state::Physics, NumFlux, State};

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
    /// // define the config instance
    /// let config = CorriesConfig::default_riemann_test::<N, E, S>(t_end, folder_name, file_name);
    ///
    /// let (mut u, mut rhs, mut time, mesh, _) = init_corries::<P, N, T, E, S>(&config).unwrap();
    /// let mut time_solver = T::new(&config).unwrap();
    /// time_solver.next_solution(&mut time.timestep, &mut u, &mut rhs, &mesh).unwrap();
    /// ```
    fn next_solution<N: NumFlux<E, S>>(
        &mut self,
        time: &mut TimeStep,
        u: &mut State<P, E, S>,
        rhs: &mut Rhs<N, E, S>,
        mesh: &Mesh<S>,
    ) -> Result<()>;
}

/// Struct for everything regarding time step integration.
///
/// Carries a [TimeStep] for keeping track of the time coordinate and data regarding it, as well as
/// an embedded type implementing [TimeSolver] that provides methods for the exact time integration
/// scheme.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Time<P: Physics<E, S>, T: TimeSolver<P, E, S>, const E: usize, const S: usize> {
    /// Keeps track of the time coordinate and related data
    pub timestep: TimeStep,

    /// Calculates solutions for the [State]
    solver: T,

    embedded_type: PhantomData<P>,
}

unsafe impl<P: Physics<E, S>, T: TimeSolver<P, E, S>, const E: usize, const S: usize> Send for Time<P, T, E, S> {}
unsafe impl<P: Physics<E, S>, T: TimeSolver<P, E, S>, const E: usize, const S: usize> Sync for Time<P, T, E, S> {}

impl<P: Physics<E, S>, T: TimeSolver<P, E, S>, const E: usize, const S: usize> Time<P, T, E, S> {
    /// Constructs a new [Time] struct.
    ///
    /// # Arguments
    ///
    /// * `config` - a [CorriesConfig] configuration object
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
    /// // define the config instance
    /// let config = CorriesConfig::default_riemann_test::<N, E, S>(t_end, folder_name, file_name);
    ///
    /// let time = Time::<P,T,E,S>::new(&config).unwrap();
    /// ```
    pub fn new(config: &CorriesConfig) -> Result<Self> {
        Ok(Self {
            timestep: TimeStep::new(&config.numerics_config, config.output_counter_max),
            solver: T::new(config)?,
            embedded_type: PhantomData,
        })
    }

    /// Calculates the next state for the [State] object `u`.
    ///
    /// # Arguments
    ///
    /// * `u` - the [State] being modified to transition between current and next state
    /// * `rhs` - solves the right-hand side
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
    /// // define the config instance
    /// let config = CorriesConfig::default_riemann_test::<N, E, S>(t_end, folder_name, file_name);
    /// let (mut u, mut rhs, mut time, mesh, _) = init_corries::<P, N, T, E, S>(&config).unwrap();
    ///
    /// time.next_solution(&mut u, &mut rhs, &mesh);
    /// ```
    pub fn next_solution<N: NumFlux<E, S>>(
        &mut self,
        u: &mut State<P, E, S>,
        rhs: &mut Rhs<N, E, S>,
        mesh: &Mesh<S>,
    ) -> Result<()> {
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
        Ok(())
    }
}
