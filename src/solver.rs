// Copyright (c) 2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [Solver] struct, responsible for generating new solutions for [State] objects

use crate::{timestep::TimeStep, CorriesConfig, Mesh, NumFlux, Physics, Rhs, State, TimeSolver};
use color_eyre::{
    eyre::{bail, Context},
    Result,
};
use std::marker::PhantomData;

/// Struct that is responsible for generating new solutions for [State] objects.
///
/// Acts as a bundler struct for [TimeStep], [Rhs], and [TimeSolver].
pub struct Solver<P: Physics<E, S>, N: NumFlux<E, S>, T: TimeSolver<P, E, S>, const E: usize, const S: usize> {
    /// Keeps track of the time coordinate and related data
    pub timestep: TimeStep,

    /// Solves time integration
    time_solver: T,

    /// Solves the right-hand side
    pub rhs: Rhs<N, E, S>,

    phantom_type: PhantomData<P>,
}

unsafe impl<P, N, T, const E: usize, const S: usize> Send for Solver<P, N, T, E, S>
where
    P: Physics<E, S>,
    N: NumFlux<E, S>,
    T: TimeSolver<P, E, S>,
{
}
unsafe impl<P, N, T, const E: usize, const S: usize> Sync for Solver<P, N, T, E, S>
where
    P: Physics<E, S>,
    N: NumFlux<E, S>,
    T: TimeSolver<P, E, S>,
{
}

impl<P: Physics<E, S>, N: NumFlux<E, S>, T: TimeSolver<P, E, S>, const E: usize, const S: usize> Solver<P, N, T, E, S> {
    /// Constructs a new [Solver] struct.
    ///
    /// # Arguments
    ///
    /// * `config` - configuration for this simulation
    /// * `mesh` - Information about spatial properties
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    ///
    /// // Set up a simulation using isothermal Euler physics on a 100 cell mesh, using the Hll and
    /// // Runge-Kutta-Fehlberg schemes for solving the equations.
    /// const S: usize = 100;
    /// set_Physics_and_E!(Euler1DIsot);
    /// type N = Hll<E, S>;
    /// type T = RungeKuttaFehlberg<P,E,S>;
    ///
    /// // use the default config for Riemann tests
    /// let t_end = 0.5;
    /// let folder_name = "results";
    /// let file_name = "noh";
    ///
    /// // define the config instance
    /// let config = CorriesConfig::default_riemann_test::<N, E, S>(t_end, folder_name, file_name);
    /// let mesh = Mesh::<S>::new(&config.mesh_config).unwrap();
    ///
    /// let solver = Solver::<P,N,T,E,S>::new(&config, &mesh).unwrap();
    /// ```
    pub fn new(config: &CorriesConfig, mesh: &Mesh<S>) -> Result<Self> {
        Ok(Self {
            timestep: TimeStep::new(&config.numerics_config, config.output_counter_max),
            time_solver: T::new(config)?,
            rhs: Rhs::<N, E, S>::new(config, mesh)?,
            phantom_type: PhantomData,
        })
    }

    /// Calculates the next state for the [State] object `u`.
    ///
    /// # Arguments
    ///
    /// * `u` - the [State] being modified to transition between current and next state
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
    /// type T = RungeKuttaFehlberg<P,E,S>;
    ///
    /// // use the default config for Riemann tests
    /// let t_end = 0.5;
    /// let folder_name = "results";
    /// let file_name = "noh";
    ///
    /// // define the config instance
    /// let config = CorriesConfig::default_riemann_test::<N, E, S>(t_end, folder_name, file_name);
    /// let (mut u, mut solver, mesh, _) = config.init_corries::<P, N, T, E, S>(|_,_,_| Ok(())).unwrap();
    ///
    /// solver.next_solution(&mut u, &mesh).unwrap();
    /// ```
    pub fn next_solution(&mut self, u: &mut State<P, E, S>, mesh: &Mesh<S>) -> Result<()> {
        self.time_solver
            .next_solution::<N>(&mut self.timestep, u, &mut self.rhs, mesh)
            .context("Calling TimeIntegration::solver.next_solution in Solver::next_solution")?;
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
