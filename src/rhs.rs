// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [Rhs] struct that carries objects and methods for solving the right-hand side of a
//! set of equations.

use color_eyre::{
    eyre::{ensure, Context},
    Result,
};
use ndarray::Array2;

use crate::{config::CorriesConfig, errorhandling::Validation, mesh::Mesh, physics::Physics, boundaryconditions::{BoundaryCondition, init_boundary_condition, Direction}, update_everything_from_cons};

// use self::numflux::{init_numflux, NumFlux};

// pub mod numflux;

/// Carries objects and methods for solving the right-hand side of a set of equations.
pub struct Rhs<P: Physics, const S: usize> {
    /// Full summed up rhs
    pub full_rhs: Array2<f64>,

    /// Calculates the numerical flux
    // numflux: Box<dyn NumFlux<S, EQ>>,

    /// Stores the numerical flux derivative along xi
    dflux_dxi: Array2<f64>,

    /// Boundary condition operator for the west boundary
    pub boundary_west: Box<dyn BoundaryCondition<P, S>>,

    /// Boundary condition operator for the east boundary
    pub boundary_east: Box<dyn BoundaryCondition<P, S>>,
}

impl<P: Physics + 'static, const S: usize> Rhs<P, S> {
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
            // numflux: Box::new(init_numflux(&config.numericsconfig)),
            dflux_dxi: Array2::zeros((E, S)),
            boundary_west: Box::new(init_boundary_condition::<P, S>(Direction::West, &config)),
            boundary_east: Box::new(init_boundary_condition::<P, S>(Direction::East, &config)),
        };
    }

    /// Solves the right-hand side and updates the `full_rhs` field
    ///
    /// # Arguments
    ///
    /// * `u` - The current [Physics] state
    /// * `mesh` - Information about spatial properties
    pub fn update(&mut self, u: &mut P, mesh: &Mesh<S>) -> Result<()> {
        self.update_dflux_dxi(u, mesh)
            .context("Calling Rhs::update_dflux_dxi in Rhs::update")?;
        self.full_rhs.assign(&self.dflux_dxi);
        return Ok(());
    }

    /// Updates the numerical flux derivative `dflux_dxi`
    ///
    /// # Arguments
    ///
    /// * `u` - The current [Physics] state
    /// * `mesh` - Information about spatial properties
    fn update_dflux_dxi(&mut self, u: &mut P, mesh: &Mesh<S>) -> Result<()> {
        // this assumes that u.cons is up-to-date
        update_everything_from_cons(u, &mut self.boundary_west, &mut self.boundary_east, mesh);
        // self.numflux
        //     .calc_dflux_dxi(&mut self.dflux_dxi, u, mesh)
        //     .context("Calling Rhs::numflux::calc_dflux_dxi in Rhs::update_dflux_dxi")?;
        return Ok(());
    }
}

impl<P: Physics, const S: usize> Validation for Rhs<P, S> {
    fn validate(&self) -> Result<()> {
        ensure!(
            self.dflux_dxi.fold(true, |acc, x| acc && x.is_finite()),
            "Rhs::dflux_dxi must be finite! Got: {}",
            self.dflux_dxi
        );
        ensure!(
            self.full_rhs.fold(true, |acc, x| acc && x.is_finite()),
            "Rhs::full_rhs must be finite! Got: {}",
            self.full_rhs
        );
        return Ok(());
    }
}
