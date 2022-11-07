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

use crate::{
    boundaryconditions::BoundaryConditionContainer, config::CorriesConfig, errorhandling::Validation, mesh::Mesh,
    physics::Physics,
};

use self::numflux::{init_numflux, NumFlux};

pub mod numflux;

/// Carries objects and methods for solving the right-hand side of a set of equations.
pub struct Rhs<const S: usize, const EQ: usize> {
    /// Full summed up rhs
    pub full_rhs: Array2<f64>,

    /// Calculates the numerical flux
    numflux: Box<dyn NumFlux<S, EQ>>,

    /// Stores the numerical flux derivative along xi
    dflux_dxi: Array2<f64>,

    /// Stores the boundary conditions for this simulation
    pub boundary_conditions: BoundaryConditionContainer<S, EQ>,
}

impl<const S: usize, const EQ: usize> Rhs<S, EQ> {
    /// Constructs a new [Rhs] object.
    ///
    /// # Arguments
    ///
    /// * `config` - Configuration for the whole simulation
    /// * `u` - The main [Physics] object for this simulation
    /// * `numflux` - The [dyn NumFlux] object about to be stored in this [Rhs] object
    pub fn new(config: &CorriesConfig, u: &Physics<S, EQ>) -> Self {
        return Rhs {
            full_rhs: Array2::zeros((EQ, S)),
            numflux: Box::new(init_numflux(&config.numericsconfig)),
            dflux_dxi: Array2::zeros((EQ, S)),
            boundary_conditions: BoundaryConditionContainer::new(config, u),
        };
    }

    /// Solves the right-hand side and updates the `full_rhs` field
    ///
    /// # Arguments
    ///
    /// * `u` - The current [Physics] state
    /// * `mesh` - Information about spatial properties
    pub fn update(&mut self, u: &mut Physics<S, EQ>, mesh: &Mesh<S>) -> Result<()> {
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
    fn update_dflux_dxi(&mut self, u: &mut Physics<S, EQ>, mesh: &Mesh<S>) -> Result<()> {
        // this assumes that u.cons is up-to-date; self.update_physics updates u.prim and the
        // derived variables anyways.
        u.update_everything_from_cons(&mut self.boundary_conditions, mesh)
            .context("Calling Physics::update_everything_from_prim in Rhs::update_dflux_dxi")?;
        self.numflux
            .calc_dflux_dxi(&mut self.dflux_dxi, u, mesh)
            .context("Calling Rhs::numflux::calc_dflux_dxi in Rhs::update_dflux_dxi")?;
        return Ok(());
    }
}

impl<const S: usize, const EQ: usize> Validation for Rhs<S, EQ> {
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
