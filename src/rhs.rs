// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [Rhs] struct that carries objects and methods for solving the right-hand side of a
//! set of equations.

pub use self::numflux::{hll::Hll,kt::Kt, NumFlux};
use crate::{
    boundaryconditions::{init_boundary_condition, BoundaryCondition},
    errorhandling::Validation,
    prelude::*,
};
use color_eyre::{
    eyre::{ensure, Context},
    Result,
};
use ndarray::Array2;

pub mod numflux;

/// Carries objects and methods for solving the right-hand side of a set of equations.
pub struct Rhs<N: NumFlux<E, S>, const E: usize, const S: usize> {
    /// Full summed up rhs
    pub full_rhs: Array2<f64>,

    /// Calculates the numerical flux
    numflux: N,

    /// Boundary condition operator for the west boundary
    pub boundary_west: Box<dyn BoundaryCondition<E, S>>,

    /// Boundary condition operator for the east boundary
    pub boundary_east: Box<dyn BoundaryCondition<E, S>>,
}

unsafe impl<N: NumFlux<E, S>, const E: usize, const S: usize> Send for Rhs<N, E, S> {}
unsafe impl<N: NumFlux<E, S>, const E: usize, const S: usize> Sync for Rhs<N, E, S> {}

impl<N: NumFlux<E, S>, const E: usize, const S: usize> Rhs<N, E, S> {
    /// Constructs a new [Rhs] object.
    ///
    /// # Arguments
    ///
    /// * `config` - Configuration for the whole simulation
    /// * `u` - The main [State] object for this simulation
    /// * `numflux` - The [NumFlux] object about to be stored in this [Rhs] object
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
    ///
    /// // use the default config for Riemann tests
    /// let t_end = 0.5;
    /// let folder_name = "results";
    /// let file_name = "noh";
    ///
    /// // define the config instance
    /// let config = CorriesConfig::default_riemann_test::<N, E, S>(t_end, folder_name, file_name);
    /// let mesh: Mesh<S> = Mesh::<S>::new(&config.mesh_config).unwrap();
    ///
    /// let rhs: Rhs<N, E, S> = Rhs::<N,E,S>::new(&config, &mesh).unwrap();
    /// ```
    pub fn new(config: &CorriesConfig, mesh: &Mesh<S>) -> Result<Self> {
        Ok(Rhs {
            full_rhs: Array2::zeros((E, S)),
            numflux: N::new(&config.numerics_config.numflux_config, mesh)?,
            boundary_west: Box::new(init_boundary_condition::<E, S>(Direction::West, config)),
            boundary_east: Box::new(init_boundary_condition::<E, S>(Direction::East, config)),
        })
    }

    /// Solves the right-hand side and updates the `full_rhs` field
    ///
    /// # Arguments
    ///
    /// * `u` - The current [Physics] state
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
    /// let (mut u, mut rhs, _, mesh, _) = init_corries::<P, N, T, E, S>(&config).unwrap();
    ///
    /// rhs.update(&mut u, &mesh).unwrap();
    /// ```
    pub fn update<P: Physics<E, S>>(&mut self, u: &mut State<P, E, S>, mesh: &Mesh<S>) -> Result<()> {
        // this assumes that u.cons is up-to-date
        u.update_vars_from_cons(&mut self.boundary_west, &mut self.boundary_east, mesh);
        self.numflux
            .calc_dflux_dxi(&mut self.full_rhs, u, mesh)
            .context("Calling Rhs::numflux::calc_dflux_dxi in Rhs::update_dflux_dxi")?;
        Ok(())
    }
}

impl<N: NumFlux<E, S>, const E: usize, const S: usize> Validation for Rhs<N, E, S> {
    fn validate(&self) -> Result<()> {
        ensure!(
            self.full_rhs.fold(true, |acc, x| acc && x.is_finite()),
            "Rhs::full_rhs must be finite! Got: {}",
            self.full_rhs
        );
        Ok(())
    }
}
