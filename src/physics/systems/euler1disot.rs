// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

use color_eyre::{eyre::ensure, Result};
use ndarray::{ArrayView1, ArrayView2, ArrayViewMut2, Zip};

use crate::{
    config::physicsconfig::PhysicsConfig, data::Data, errorhandling::Validation, physics::Physics,
    variables::Variables, Collectable,
};

const E: usize = 2;
const JRHO: usize = 0;
const JXI: usize = 1;
const J_EIGENMIN: usize = 0;
const J_EIGENMAX: usize = 1;

/// Test struct for using a trait for Physics
#[derive(Debug)]
pub struct Euler1DIsot<const S: usize> {
    /// Variables at the centre of the mesh's cells
    pub cent: Variables<E, S>,

    /// Variables at the centre of the mesh's cells
    pub west: Variables<E, S>,

    /// Variables at the centre of the mesh's cells
    pub east: Variables<E, S>,
}

impl<const S: usize> Physics for Euler1DIsot<S> {
    const E: usize = E;
    const S: usize = S;
    type Vars = Variables<E, S>;

    fn new(physics_config: &PhysicsConfig) -> Self {
        return Self {
            cent: Variables::new(physics_config),
            west: Variables::new(physics_config),
            east: Variables::new(physics_config),
        };
    }

    fn cent<'a>(&self) -> &'a Self::Vars {
        &self.cent
    }

    fn west<'a>(&self) -> &'a Self::Vars {
        &self.west
    }

    fn east<'a>(&self) -> &'a Self::Vars {
        &self.east
    }

    fn is_adiabatic(&self) -> bool {
        return false;
    }

    fn update_prim(&mut self) {
        cons_to_prim(
            &mut self.cent().prim.view_mut(),
            JRHO,
            JXI,
            &self.cent().cons.view(),
            JRHO,
            JXI,
        );
    }

    fn update_cons(&mut self) {
        prim_to_cons(
            &mut self.cent().cons.view_mut(),
            JRHO,
            JXI,
            &self.cent().prim.view(),
            JRHO,
            JXI,
        );
    }

    fn update_derived_values(&mut self) {
        update_eigen_vals(
            self.cent().eigen_vals.view_mut(),
            J_EIGENMIN,
            J_EIGENMAX,
            self.cent().prim.row(JXI),
            self.cent().c_sound.view(),
        );
    }

    fn update_flux_cent(&mut self) {
        update_flux(
            self.cent().flux.view_mut(),
            self.cent().prim.view(),
            self.cent().cons.view(),
            JRHO,
            JXI,
            self.cent().c_sound.view(),
        );
    }
}

/// Updates eigen values
#[inline(always)]
pub fn update_eigen_vals(
    mut eigen_vals: ArrayViewMut2<f64>,
    j_eigen_min: usize,
    j_eigen_max: usize,
    xi_vel: ArrayView1<f64>,
    c_sound: ArrayView1<f64>,
) {
    Zip::from(eigen_vals.row_mut(j_eigen_min))
        .and(xi_vel)
        .and(c_sound)
        .for_each(|e, &v, &cs| *e = v - cs);

    for j in j_eigen_min + 1..j_eigen_max {
        eigen_vals.row_mut(j).assign(&xi_vel);
    }

    Zip::from(eigen_vals.row_mut(j_eigen_max))
        .and(xi_vel)
        .and(c_sound)
        .for_each(|e, &v, &cs| *e = v + cs);
}

/// Updates physical flux
#[inline(always)]
pub fn update_flux(
    mut flux: ArrayViewMut2<f64>,
    prim: ArrayView2<f64>,
    cons: ArrayView2<f64>,
    j_rho: usize,
    j_xi: usize,
    c_sound: ArrayView1<f64>,
) {
    Zip::from(flux.row_mut(j_rho))
        .and(cons.row(j_xi))
        .for_each(|f, &xi_mom| *f = xi_mom);

    Zip::from(flux.row_mut(j_xi))
        .and(cons.row(j_rho))
        .and(prim.row(j_xi))
        .and(c_sound)
        .for_each(|f, &rho, &xi_vel, &cs| *f = rho * (xi_vel * xi_vel + cs * cs));
}

/// Converts conservative to primitive variables
#[inline(always)]
pub fn cons_to_prim(
    prim: &mut ArrayViewMut2<f64>,
    j_rho_prim: usize,
    j_xi_vel: usize,
    cons: &ArrayView2<f64>,
    j_rho_cons: usize,
    j_xi_mom: usize,
) {
    Zip::from(prim.row_mut(j_rho_prim))
        .and(cons.row(j_rho_cons))
        .for_each(|rho_prim, &rho_cons| *rho_prim = rho_cons);
    Zip::from(prim.row_mut(j_xi_vel))
        .and(cons.row(j_xi_mom))
        .and(cons.row(j_rho_cons))
        .for_each(|xi_vel, &xi_mom, &rho_cons| *xi_vel = xi_mom / rho_cons);
}

/// Converts primitive to conservative variables
#[inline(always)]
pub fn prim_to_cons(
    cons: &mut ArrayViewMut2<f64>,
    j_rho_cons: usize,
    j_xi_mom: usize,
    prim: &ArrayView2<f64>,
    j_rho_prim: usize,
    j_xi_vel: usize,
) {
    Zip::from(cons.row_mut(j_rho_cons))
        .and(prim.row(j_rho_prim))
        .for_each(|rho_cons, &rho_prim| *rho_cons = rho_prim);
    Zip::from(cons.row_mut(j_xi_mom))
        .and(prim.row(j_xi_vel))
        .and(prim.row(j_rho_prim))
        .for_each(|xi_mom, &xi_vel, &rho_cons| *xi_mom = xi_vel * rho_cons);
}

impl<const S: usize> Collectable for Euler1DIsot<S> {
    fn collect_data(&self, name: &mut Data, mesh_offset: usize) -> Result<()> {
        return self.cent().collect_data(name, mesh_offset);
    }
}

impl<const S: usize> Validation for Euler1DIsot<S> {
    fn validate(&self) -> Result<()> {
        self.cent().validate()?;
        validate(self, JRHO)?;
        return Ok(());
    }
}

#[inline(always)]
pub fn validate<P: Physics>(u: &P, j_rho: usize) -> Result<()> {
    ensure!(
        u.cent().prim_row(j_rho).fold(true, |acc, x| acc && x > &0.0),
        "Mass density must be positive! Got: {}",
        u.cent().prim_row(j_rho)
    );
    return Ok(());
}
