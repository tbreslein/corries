// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

use ndarray::{Array1, Array2, ArrayView1, ArrayView2, ArrayViewMut2, Zip};

use crate::{config::physicsconfig::PhysicsConfig, physics::Physics};

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
            prim: Array2::zeros((2, S)),
            cons: Array2::zeros((2, S)),
            c_sound: Array1::zeros(S),
            eigen_vals: Array2::zeros((2, S)),
            flux: Array2::zeros((2, S)),
        };
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
        return self.eigen_vals.row(0);
    }

    #[inline(always)]
    fn eigen_max(&self) -> ArrayView1<f64> {
        return self.eigen_vals.row(1);
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
        cons_to_prim(self.prim.view_mut(), self.cons.view());
    }

    fn update_cons(&mut self) {
        prim_to_cons(self.cons.view_mut(), self.prim.view());
    }

    fn update_derived_values(&mut self) {
        update_eigen_vals(self.eigen_vals.view_mut(), self.prim.row(1), self.c_sound.view());
        update_flux(
            self.flux.view_mut(),
            self.prim.view(),
            self.cons.view(),
            self.c_sound.view(),
        );
    }

    #[inline(always)]
    fn assign_prim_entry(&mut self, rhs: f64, j: usize, i: usize) {
        self.prim[[j, i]] = rhs;
    }

    #[inline(always)]
    fn assign_cons_entry(&mut self, rhs: f64, j: usize, i: usize) {
        self.cons[[j, i]] = rhs;
    }

    #[inline(always)]
    fn assign_prim(&mut self, rhs: &Array2<f64>) {
        self.prim.assign(rhs);
    }

    #[inline(always)]
    fn assign_cons(&mut self, rhs: &Array2<f64>) {
        self.cons.assign(rhs);
    }
}

/// Updates eigen values
#[inline(always)]
pub fn update_eigen_vals(mut eigen_vals: ArrayViewMut2<f64>, xi_vel: ArrayView1<f64>, c_sound: ArrayView1<f64>) {
    Zip::from(eigen_vals.row_mut(0))
        .and(xi_vel)
        .and(c_sound)
        .for_each(|e, &v, &cs| *e = v - cs);

    Zip::from(eigen_vals.row_mut(1))
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
    c_sound: ArrayView1<f64>,
) {
    Zip::from(flux.row_mut(0))
        .and(cons.row(1))
        .for_each(|f, &xi_mom| *f = xi_mom);

    Zip::from(flux.row_mut(1))
        .and(cons.row(0))
        .and(prim.row(1))
        .and(c_sound)
        .for_each(|f, &rho, &xi_vel, &cs| *f = rho * (xi_vel * xi_vel + cs * cs));
}

/// Converts conservative to primitive variables
#[inline(always)]
pub fn cons_to_prim(mut prim: ArrayViewMut2<f64>, cons: ArrayView2<f64>) {
    Zip::from(prim.row_mut(0))
        .and(cons.row(0))
        .for_each(|rho_prim, &rho_cons| *rho_prim = rho_cons);

    Zip::from(prim.row_mut(1))
        .and(cons.row(0))
        .and(cons.row(1))
        .for_each(|xi_vel, &xi_mom, &rho_cons| *xi_vel = xi_mom / rho_cons);
}

/// Converts primitive to conservative variables
#[inline(always)]
pub fn prim_to_cons(mut cons: ArrayViewMut2<f64>, prim: ArrayView2<f64>) {
    Zip::from(cons.row_mut(0))
        .and(prim.row(0))
        .for_each(|rho_cons, &rho_prim| *rho_cons = rho_prim);

    Zip::from(cons.row_mut(1))
        .and(prim.row(0))
        .and(prim.row(1))
        .for_each(|xi_mom, &xi_vel, &rho_cons| *xi_mom = xi_vel / rho_cons);
}
