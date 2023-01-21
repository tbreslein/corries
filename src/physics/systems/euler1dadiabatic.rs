// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

use color_eyre::{eyre::ensure, Result};
use ndarray::{ArrayView1, ArrayView2, ArrayViewMut1, ArrayViewMut2, Zip};

use crate::{
    config::physicsconfig::PhysicsConfig, errorhandling::Validation, physics::Physics, variables::Variables,
    Collectable, Data,
};

const E: usize = 3;
const JRHO: usize = 0;
const JXI: usize = 1;
const JP: usize = 2;
const J_EIGENMIN: usize = 0;
const J_EIGENMAX: usize = 2;

/// Test struct for using a trait for Physics
#[derive(Debug)]
pub struct Euler1DAdiabatic<const S: usize> {
    /// Variables at the centre of the mesh's cells
    pub cent: Variables<E, S>,

    /// Variables at the centre of the mesh's cells
    pub west: Variables<E, S>,

    /// Variables at the centre of the mesh's cells
    pub east: Variables<E, S>,
}

impl<const S: usize> Physics for Euler1DAdiabatic<S> {
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
        return true;
    }

    fn update_prim(&mut self) {
        cons_to_prim(
            &mut self.cent().prim.view_mut(),
            JRHO,
            JXI,
            JP,
            &self.cent().cons.view(),
            JRHO,
            JXI,
            JP,
            self.cent().gamma,
        );
    }

    fn update_cons(&mut self) {
        prim_to_cons(
            &mut self.cent().cons.view_mut(),
            JRHO,
            JXI,
            JP,
            &self.cent().prim.view(),
            JRHO,
            JXI,
            JP,
            self.cent().gamma,
        );
    }

    fn update_derived_values(&mut self) {
        update_c_sound(
            self.cent().c_sound.view_mut(),
            self.cent().gamma,
            self.cent().prim.row(JP),
            self.cent().prim.row(JRHO),
        );
        super::euler1disot::update_eigen_vals(
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
            JP,
        );
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

impl<const S: usize> Collectable for Euler1DAdiabatic<S> {
    fn collect_data(&self, name: &mut Data, mesh_offset: usize) -> Result<()> {
        return self.cent().collect_data(name, mesh_offset);
    }
}

impl<const S: usize> Validation for Euler1DAdiabatic<S> {
    fn validate(&self) -> Result<()> {
        self.cent().validate()?;
        super::euler1disot::validate(self, JRHO)?;
        validate(self, JP)?;
        return Ok(());
    }
}

#[inline(always)]
pub fn validate<P: Physics>(u: &P, j_pressure: usize) -> Result<()> {
    ensure!(
        u.cent().prim_row(j_pressure).fold(true, |acc, x| acc && x > &0.0),
        "Pressure must be positive! Got: {}",
        u.cent().prim_row(j_pressure)
    );
    return Ok(());
}
