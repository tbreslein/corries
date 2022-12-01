// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

use ndarray::{Array1, Array2, ArrayView1, ArrayView2, ArrayViewMut1, ArrayViewMut2, Zip};

use crate::{config::physicsconfig::PhysicsConfig, physics::Physics};

/// Test struct for using a trait for Physics
#[derive(Debug)]
pub struct Euler1DAdiabatic<const S: usize> {
    /// Primitive variables
    pub prim: Array2<f64>,

    /// Conservative variables
    pub cons: Array2<f64>,

    /// Speed of sound
    c_sound: Array1<f64>,

    /// Eigen values
    eigen_vals: Array2<f64>,

    /// Physical flux
    flux: Array2<f64>,

    /// Adiabatic index
    gamma: f64,
}

impl<const S: usize> Physics for Euler1DAdiabatic<S> {
    fn new(physics_config: &PhysicsConfig) -> Self {
        return Self {
            prim: Array2::zeros((3, S)),
            cons: Array2::zeros((3, S)),
            c_sound: Array1::zeros(S),
            eigen_vals: Array2::zeros((3, S)),
            flux: Array2::zeros((3, S)),
            gamma: physics_config.adiabatic_index,
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
        cons_to_prim(&mut self.prim.view_mut(), &self.cons.view(), self.gamma);
    }

    fn update_cons(&mut self) {
        prim_to_cons(&mut self.cons.view_mut(), &self.prim.view(), self.gamma);
    }

    fn update_derived_values(&mut self) {
        update_c_sound(self.c_sound.view_mut(), self.gamma, self.prim.row(2), self.prim.row(0));
        update_eigen_vals(self.eigen_vals.view_mut(), self.prim.row(1), self.c_sound.view());
        update_flux(self.flux.view_mut(), self.prim.view(), self.cons.view());
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
    fn assign_prim(&mut self, rhs: &Array2<f64>) {
        self.prim.assign(rhs);
    }

    #[inline(always)]
    fn assign_cons(&mut self, rhs: &Array2<f64>) {
        self.cons.assign(rhs);
    }
}

/// Update c_sound
#[inline(always)]
pub fn update_c_sound(c_sound: ArrayViewMut1<f64>, gamma: f64, pressure: ArrayView1<f64>, rho: ArrayView1<f64>) {
    Zip::from(c_sound)
        .and(pressure)
        .and(rho)
        .for_each(|cs, &p, &rho| *cs = (gamma * p / rho).sqrt());
}

/// Updates eigen values
#[inline(always)]
pub fn update_eigen_vals(mut eigen_vals: ArrayViewMut2<f64>, xi_vel: ArrayView1<f64>, c_sound: ArrayView1<f64>) {
    Zip::from(eigen_vals.row_mut(0))
        .and(xi_vel)
        .and(c_sound)
        .for_each(|e, &v, &cs| *e = v - cs);

    eigen_vals.row_mut(1).assign(&xi_vel);

    Zip::from(eigen_vals.row_mut(2))
        .and(xi_vel)
        .and(c_sound)
        .for_each(|e, &v, &cs| *e = v + cs);
}

/// Updates physical flux
#[inline(always)]
pub fn update_flux(mut flux: ArrayViewMut2<f64>, prim: ArrayView2<f64>, cons: ArrayView2<f64>) {
    Zip::from(flux.row_mut(0))
        .and(cons.row(1))
        .for_each(|f, &xi_mom| *f = xi_mom);

    Zip::from(flux.row_mut(1))
        .and(prim.row(1))
        .and(cons.row(1))
        .and(prim.row(2))
        .for_each(|f, &xi_vel, &xi_mom, &pressure| *f = xi_vel * xi_mom + pressure);

    Zip::from(flux.row_mut(2))
        .and(cons.row(2))
        .and(prim.row(2))
        .and(prim.row(1))
        .for_each(|f, &energy, &pressure, &xi_vel| *f = (energy + pressure) * xi_vel);
}

/// Converts conservative to primitive variables
#[inline(always)]
pub fn cons_to_prim(prim: &mut ArrayViewMut2<f64>, cons: &ArrayView2<f64>, gamma: f64) {
    super::euler1disot::cons_to_prim(prim, cons);
    Zip::from(prim.row_mut(2))
        .and(cons.row(2))
        .and(cons.row(0))
        .and(cons.row(1))
        .for_each(|pressure, &energy, &rho_cons, &xi_mom| {
            *pressure = (gamma - 1.0) * (energy - 0.5 / rho_cons * xi_mom * xi_mom)
        });
}

/// Converts primitive to conservative variables
#[inline(always)]
pub fn prim_to_cons(cons: &mut ArrayViewMut2<f64>, prim: &ArrayView2<f64>, gamma: f64) {
    super::euler1disot::prim_to_cons(cons, prim);
    Zip::from(cons.row_mut(2))
        .and(prim.row(2))
        .and(prim.row(0))
        .and(prim.row(1))
        .for_each(|energy, &pressure, &rho_prim, &xi_vel| {
            *energy = pressure * (gamma - 1.0).recip() + 0.5 * rho_prim * xi_vel * xi_vel
        });
}
