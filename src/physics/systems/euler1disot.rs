// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

use color_eyre::{eyre::ensure, Result};
use ndarray::{Array1, Array2, ArrayView1, ArrayView2, ArrayViewMut2, Zip};

use crate::{
    config::physicsconfig::PhysicsConfig, data::Data, errorhandling::Validation, physics::Physics, Collectable,
};

const E: usize = 2;
const JRHO: usize = 0;
const JXI: usize = 1;
const J_EIGENMIN: usize = 0;
const J_EIGENMAX: usize = 1;

/// Test struct for using a trait for Physics
#[derive(Debug)]
pub struct Euler1DIsot<const S: usize> {
    /// Primitive variables
    pub prim: Array2<f64>,

    /// Conservative variables
    pub cons: Array2<f64>,

    /// Speed of sound
    pub c_sound: Array1<f64>,

    /// Eigen values
    pub eigen_vals: Array2<f64>,

    /// Physical flux
    pub flux: Array2<f64>,
}

impl<const S: usize> Physics for Euler1DIsot<S> {
    fn new(_: &PhysicsConfig) -> Self {
        return Self {
            prim: Array2::zeros((E, S)),
            cons: Array2::zeros((E, S)),
            c_sound: Array1::zeros(S),
            eigen_vals: Array2::zeros((E, S)),
            flux: Array2::zeros((E, S)),
        };
    }

    fn is_adiabatic(&self) -> bool {
        return false;
    }

    #[inline(always)]
    fn prim_entry(&self, j: usize, i: usize) -> f64 {
        return self.prim[[j, i]];
    }

    #[inline(always)]
    fn prim_row(&self, j: usize) -> ArrayView1<f64> {
        return self.prim.row(j);
    }

    #[inline(always)]
    fn prim(&self) -> ArrayView2<f64> {
        return self.prim.view();
    }

    #[inline(always)]
    fn cons_entry(&self, j: usize, i: usize) -> f64 {
        return self.cons[[j, i]];
    }

    #[inline(always)]
    fn cons_row(&self, j: usize) -> ArrayView1<f64> {
        return self.cons.row(j);
    }

    #[inline(always)]
    fn cons(&self) -> ArrayView2<f64> {
        return self.cons.view();
    }

    #[inline(always)]
    fn eigen_vals(&self) -> ArrayView2<f64> {
        return self.eigen_vals.view();
    }

    #[inline(always)]
    fn eigen_min(&self) -> ArrayView1<f64> {
        return self.eigen_vals.row(J_EIGENMIN);
    }

    #[inline(always)]
    fn eigen_max(&self) -> ArrayView1<f64> {
        return self.eigen_vals.row(J_EIGENMAX);
    }

    #[inline(always)]
    fn c_sound(&self) -> ArrayView1<f64> {
        return self.c_sound.view();
    }

    #[inline(always)]
    fn flux_entry(&self, j: usize, i: usize) -> f64 {
        return self.flux[[j, i]];
    }

    #[inline(always)]
    fn flux_row(&self, j: usize) -> ArrayView1<f64> {
        return self.flux.row(j);
    }

    #[inline(always)]
    fn flux(&self) -> ArrayView2<f64> {
        return self.flux.view();
    }

    fn update_prim(&mut self) {
        cons_to_prim(&mut self.prim.view_mut(), JRHO, JXI, &self.cons.view(), JRHO, JXI);
    }

    fn update_cons(&mut self) {
        prim_to_cons(&mut self.cons.view_mut(), JRHO, JXI, &self.prim.view(), JRHO, JXI);
    }

    fn update_derived_values(&mut self) {
        update_eigen_vals(
            self.eigen_vals.view_mut(),
            J_EIGENMIN,
            J_EIGENMAX,
            self.prim.row(JXI),
            self.c_sound.view(),
        );
        update_flux(
            self.flux.view_mut(),
            self.prim.view(),
            self.cons.view(),
            JRHO,
            JXI,
            self.c_sound.view(),
        );
    }

    #[inline(always)]
    fn assign_prim_entry(&mut self, j: usize, i: usize, rhs: f64) {
        self.prim[[j, i]] = rhs;
    }

    #[inline(always)]
    fn assign_cons_entry(&mut self, j: usize, i: usize, rhs: f64) {
        self.cons[[j, i]] = rhs;
    }

    #[inline(always)]
    fn assign_prim(&mut self, rhs: &ArrayView2<f64>) {
        self.prim.assign(rhs);
    }

    #[inline(always)]
    fn assign_cons(&mut self, rhs: &ArrayView2<f64>) {
        self.cons.assign(rhs);
    }

    #[inline(always)]
    fn assign_c_sound(&mut self, rhs: &ArrayView1<f64>) {
        self.c_sound.assign(rhs);
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
        return super::super::collect_data(self, name, mesh_offset);
    }
}

impl<const S: usize> Validation for Euler1DIsot<S> {
    fn validate(&self) -> Result<()> {
        super::super::validate(self)?;
        validate(self, JRHO)?;
        return Ok(());
    }
}

#[inline(always)]
pub fn validate<P: Physics>(u: &P, j_rho: usize) -> Result<()> {
    ensure!(
        u.prim_row(j_rho).fold(true, |acc, x| acc && x > &0.0),
        "Mass density must be positive! Got: {}",
        u.prim_row(j_rho)
    );
    return Ok(());
}
