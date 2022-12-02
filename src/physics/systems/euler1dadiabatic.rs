// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

use ndarray::{Array1, Array2, ArrayView1, ArrayView2, ArrayViewMut1, ArrayViewMut2, Zip};

use crate::{config::physicsconfig::PhysicsConfig, physics::Physics};

const E: usize = 3;
const JRHO: usize = 0;
const JXI: usize = 1;
const JP: usize = 2;
const J_EIGENMIN: usize = 0;
const J_EIGENMAX: usize = 2;

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
            prim: Array2::zeros((E, S)),
            cons: Array2::zeros((E, S)),
            c_sound: Array1::zeros(S),
            eigen_vals: Array2::zeros((E, S)),
            flux: Array2::zeros((E, S)),
            gamma: physics_config.adiabatic_index,
        };
    }

    fn is_adiabatic(&self) -> bool {
        return true;
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
        cons_to_prim(
            &mut self.prim.view_mut(),
            JRHO,
            JXI,
            JP,
            &self.cons.view(),
            JRHO,
            JXI,
            JP,
            self.gamma,
        );
    }

    fn update_cons(&mut self) {
        prim_to_cons(
            &mut self.cons.view_mut(),
            JRHO,
            JXI,
            JP,
            &self.prim.view(),
            JRHO,
            JXI,
            JP,
            self.gamma,
        );
    }

    fn update_derived_values(&mut self) {
        update_c_sound(
            self.c_sound.view_mut(),
            self.gamma,
            self.prim.row(JP),
            self.prim.row(JRHO),
        );
        super::euler1disot::update_eigen_vals(
            self.eigen_vals.view_mut(),
            J_EIGENMIN,
            J_EIGENMAX,
            self.prim.row(JXI),
            self.c_sound.view(),
        );
        update_flux(self.flux.view_mut(), self.prim.view(), self.cons.view(), JRHO, JXI, JP);
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

/// Update c_sound
#[inline(always)]
pub fn update_c_sound(c_sound: ArrayViewMut1<f64>, gamma: f64, pressure: ArrayView1<f64>, rho: ArrayView1<f64>) {
    Zip::from(c_sound)
        .and(pressure)
        .and(rho)
        .for_each(|cs, &p, &rho| *cs = (gamma * p / rho).sqrt());
}

/// Updates physical flux
#[inline(always)]
pub fn update_flux(
    mut flux: ArrayViewMut2<f64>,
    prim: ArrayView2<f64>,
    cons: ArrayView2<f64>,
    j_rho: usize,
    j_xi: usize,
    j_p: usize,
) {
    Zip::from(flux.row_mut(j_rho))
        .and(cons.row(j_xi))
        .for_each(|f, &xi_mom| *f = xi_mom);

    Zip::from(flux.row_mut(j_xi))
        .and(prim.row(j_xi))
        .and(cons.row(j_xi))
        .and(prim.row(j_p))
        .for_each(|f, &xi_vel, &xi_mom, &pressure| *f = xi_vel * xi_mom + pressure);

    Zip::from(flux.row_mut(j_p))
        .and(cons.row(j_p))
        .and(prim.row(j_p))
        .and(prim.row(j_xi))
        .for_each(|f, &energy, &pressure, &xi_vel| *f = (energy + pressure) * xi_vel);
}

/// Converts conservative to primitive variables
#[inline(always)]
pub fn cons_to_prim(
    prim: &mut ArrayViewMut2<f64>,
    j_rho_prim: usize,
    j_xi_vel: usize,
    j_pressure: usize,
    cons: &ArrayView2<f64>,
    j_rho_cons: usize,
    j_xi_mom: usize,
    j_energy: usize,
    gamma: f64,
) {
    super::euler1disot::cons_to_prim(prim, j_rho_prim, j_xi_vel, cons, j_rho_cons, j_xi_mom);
    Zip::from(prim.row_mut(j_pressure))
        .and(cons.row(j_energy))
        .and(cons.row(j_rho_cons))
        .and(cons.row(j_xi_mom))
        .for_each(|pressure, &energy, &rho_cons, &xi_mom| {
            *pressure = (gamma - 1.0) * (energy - 0.5 / rho_cons * xi_mom * xi_mom)
        });
}

/// Converts primitive to conservative variables
#[inline(always)]
pub fn prim_to_cons(
    cons: &mut ArrayViewMut2<f64>,
    j_rho_cons: usize,
    j_xi_mom: usize,
    j_energy: usize,
    prim: &ArrayView2<f64>,
    j_rho_prim: usize,
    j_xi_vel: usize,
    j_pressure: usize,
    gamma: f64,
) {
    super::euler1disot::prim_to_cons(cons, j_rho_cons, j_xi_mom, prim, j_rho_prim, j_xi_vel);
    Zip::from(cons.row_mut(j_energy))
        .and(prim.row(j_pressure))
        .and(prim.row(j_rho_prim))
        .and(prim.row(j_xi_vel))
        .for_each(|energy, &pressure, &rho_prim, &xi_vel| {
            *energy = pressure * (gamma - 1.0).recip() + 0.5 * rho_prim * xi_vel * xi_vel
        });
}
