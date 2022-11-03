// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [Rhs] struct that carries objects and methods for solving the right-hand side of a
//! set of equations.

use ndarray::Array2;

use crate::{boundaryconditions::BoundaryConditionContainer, config::CorriesConfig, mesh::Mesh, physics::Physics};

use self::numflux::NumFlux;

pub mod numflux;

/// Carries objects and methods for solving the right-hand side of a set of equations.
pub struct Rhs<'a, const S: usize, const EQ: usize> {
    /// Full summed up rhs
    pub full_rhs: Array2<f64>,

    /// Calculates the numerical flux
    numflux: &'a mut dyn numflux::NumFlux<S, EQ>,

    /// Stores the numerical flux derivative along xi
    dflux_dxi: Array2<f64>,

    /// Stores the boundary conditions for this simulation
    pub boundary_conditions: BoundaryConditionContainer<S, EQ>,
}

impl<'a, const S: usize, const EQ: usize> Rhs<'a, S, EQ> {
    /// Constructs a new [Rhs] object.
    pub fn new(config: &CorriesConfig, u: &Physics<S, EQ>, numflux: &'a mut dyn NumFlux<S, EQ>) -> Self {
        return Rhs {
            full_rhs: Array2::zeros((EQ, S)),
            numflux,
            dflux_dxi: Array2::zeros((EQ, S)),
            boundary_conditions: BoundaryConditionContainer::new(config, u),
        };
    }

    pub fn update(&mut self, u: &mut Physics<S, EQ>, mesh: &Mesh<S>) {
        self.update_dflux_dxi(u, mesh);
        self.full_rhs.assign(&self.dflux_dxi);
    }

    fn update_dflux_dxi(&mut self, u: &mut Physics<S, EQ>, mesh: &Mesh<S>) {
        // this assumes that u.cons is up-to-date; self.update_physics updates u.prim and the
        // derived variables anyways.
        self.update_physics(u, mesh);
        self.numflux.calc_dflux_dxi(&mut self.dflux_dxi, u, mesh);
    }

    pub fn update_physics(&mut self, u: &mut Physics<S, EQ>, mesh: &Mesh<S>) {
        u.update_prim();
        u.update_derived_variables();
        self.boundary_conditions.apply(u, mesh);
        u.update_cons();
    }
}
