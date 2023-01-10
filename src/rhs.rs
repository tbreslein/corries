// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [Rhs] struct that carries objects and methods for solving the right-hand side of a
//! set of equations.

use color_eyre::{
    eyre::{ensure, Context},
    Result,
};
use ndarray::Array2;

use crate::errorhandling::Validation;
use crate::{
    boundaryconditions::{init_boundary_condition, BoundaryCondition, Direction},
    prelude::*,
};

pub mod numflux;
pub use self::numflux::{hll::Hll, init_numflux, NumFlux};

/// Carries objects and methods for solving the right-hand side of a set of equations.
pub struct Rhs<P: Physics, N: NumFlux, const S: usize> {
    /// Full summed up rhs
    pub full_rhs: Array2<f64>,

    /// Calculates the numerical flux
    numflux: N,

    /// Boundary condition operator for the west boundary
    pub boundary_west: Box<dyn BoundaryCondition<P, S>>,

    /// Boundary condition operator for the east boundary
    pub boundary_east: Box<dyn BoundaryCondition<P, S>>,
}

impl<P: Physics + 'static, N: NumFlux, const S: usize> Rhs<P, N, S> {
    /// Constructs a new [Rhs] object.
    ///
    /// # Arguments
    ///
    /// * `config` - Configuration for the whole simulation
    /// * `u` - The main [Physics] object for this simulation
    /// * `numflux` - The [dyn NumFlux] object about to be stored in this [Rhs] object
    pub fn new<const E: usize>(config: &CorriesConfig) -> Self {
        return Rhs {
            full_rhs: Array2::zeros((E, S)),
            numflux: init_numflux::<N, E, S>(&config.numerics_config),
            boundary_west: Box::new(init_boundary_condition::<P, S>(Direction::West, config)),
            boundary_east: Box::new(init_boundary_condition::<P, S>(Direction::East, config)),
        };
    }

    /// Solves the right-hand side and updates the `full_rhs` field
    ///
    /// # Arguments
    ///
    /// * `u` - The current [Physics] state
    /// * `mesh` - Information about spatial properties
    pub fn update<const E: usize>(&mut self, u: &mut P, mesh: &Mesh<S>) -> Result<()> {
        // this assumes that u.cons is up-to-date
        update_everything_from_cons(u, &mut self.boundary_west, &mut self.boundary_east, mesh);
        self.numflux
            .calc_dflux_dxi::<P, E, S>(&mut self.full_rhs, u, mesh)
            .context("Calling Rhs::numflux::calc_dflux_dxi in Rhs::update_dflux_dxi")?;
        return Ok(());
    }
}

impl<P: Physics, N: NumFlux, const S: usize> Validation for Rhs<P, N, S> {
    fn validate(&self) -> Result<()> {
        ensure!(
            self.full_rhs.fold(true, |acc, x| acc && x.is_finite()),
            "Rhs::full_rhs must be finite! Got: {}",
            self.full_rhs
        );
        return Ok(());
    }
}
