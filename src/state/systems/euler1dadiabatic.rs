// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

use color_eyre::{eyre::ensure, Result};
use ndarray::{ArrayView2, ArrayViewMut2, Zip};

use crate::{state::Physics, variables::Variables};

const E: usize = 3;

/// Test struct for using a trait for Physics
#[derive(Debug)]
pub struct Euler1DAdiabatic<const S: usize>;

impl<const S: usize> Physics<E, S> for Euler1DAdiabatic<S> {
    const IS_ADIABATIC: bool = true;
    const JRHO: usize = 0;
    const JXI: usize = 1;
    const JETA: usize = std::usize::MAX;
    const JPRESSURE: usize = 2;

    fn new() -> Self {
        return Self;
    }

    #[inline(always)]
    fn update_prim(vars: &mut Variables<E, S>) {
        cons_to_prim(
            &mut vars.prim.view_mut(),
            Self::JRHO,
            Self::JXI,
            Self::JPRESSURE,
            &vars.cons.view(),
            Self::JRHO,
            Self::JXI,
            Self::JPRESSURE,
            vars.gamma,
        );
    }

    #[inline(always)]
    fn update_cons(vars: &mut Variables<E, S>) {
        prim_to_cons(
            &mut vars.cons.view_mut(),
            Self::JRHO,
            Self::JXI,
            Self::JPRESSURE,
            &vars.prim.view(),
            Self::JRHO,
            Self::JXI,
            Self::JPRESSURE,
            vars.gamma,
        );
    }

    #[inline(always)]
    fn update_flux(vars: &mut Variables<E, S>) {
        update_flux(
            vars.flux.view_mut(),
            vars.prim.view(),
            vars.cons.view(),
            Self::JRHO,
            Self::JXI,
            Self::JPRESSURE,
        );
    }

    #[inline(always)]
    fn validate(vars: &Variables<E, S>) -> Result<()> {
        super::euler1disot::validate(vars, Self::JRHO)?;
        return validate(vars, Self::JPRESSURE);
    }
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

#[inline(always)]
pub fn validate<const E: usize, const S: usize>(vars: &Variables<E, S>, j_pressure: usize) -> Result<()> {
    ensure!(
        vars.prim.row(j_pressure).fold(true, |acc, x| acc && x > &0.0),
        "Pressure must be positive! Got: {}",
        vars.prim.row(j_pressure)
    );
    return Ok(());
}
