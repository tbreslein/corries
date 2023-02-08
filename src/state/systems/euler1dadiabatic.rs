// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [Euler1DAdiabatic] struct, which is an implementer for [Physics] for 1-dimensional
//! adiabatic Euler equations

use color_eyre::{eyre::ensure, Result};
use ndarray::{ArrayView2, ArrayViewMut2, Zip};

use crate::{state::Physics, variables::Variables};

const E: usize = 3;

/// Test struct for using a trait for Physics
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Euler1DAdiabatic<const S: usize>;

unsafe impl<const S: usize> Send for Euler1DAdiabatic<S> {}
unsafe impl<const S: usize> Sync for Euler1DAdiabatic<S> {}

impl<const S: usize> Physics<E, S> for Euler1DAdiabatic<S> {
    const IS_ADIABATIC: bool = true;
    const JRHO: usize = 0;
    const JXI: usize = 1;
    const JETA: usize = std::usize::MAX;
    const JPRESSURE: usize = 2;

    fn new() -> Self {
        Self
    }

    #[inline(always)]
    fn update_prim(vars: &mut Variables<E, S>) {
        // PERF: This was benchmarked between the following options:
        // - raw index loop over columns calculating the tuple of the primitive values per column
        // (~20% speed up over the zip)
        // - ndarray::Zip!, where each row is calculated in a separate zip

        for i in 0..S {
            (
                vars.prim[[Self::JRHO, i]],
                vars.prim[[Self::JXI, i]],
                vars.prim[[Self::JPRESSURE, i]],
            ) = cons_to_prim(
                vars.cons[[Self::JRHO, i]],
                vars.cons[[Self::JXI, i]],
                vars.cons[[Self::JPRESSURE, i]],
                vars.gamma,
            );
        }
    }

    #[inline(always)]
    fn update_cons(vars: &mut Variables<E, S>) {
        // PERF: This was benchmarked between the following options:
        // - raw index loop over columns calculating the tuple of the primitive values per column
        // (~20% speed up over the zip)
        // - ndarray::Zip!, where each row is calculated in a separate zip

        for i in 0..S {
            (
                vars.cons[[Self::JRHO, i]],
                vars.cons[[Self::JXI, i]],
                vars.cons[[Self::JPRESSURE, i]],
            ) = prim_to_cons(
                vars.prim[[Self::JRHO, i]],
                vars.prim[[Self::JXI, i]],
                vars.prim[[Self::JPRESSURE, i]],
                vars.gamma,
            );
        }
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
        validate(vars, Self::JPRESSURE)
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
pub fn cons_to_prim(rho_cons: f64, xi_mom: f64, energy: f64, gamma: f64) -> (f64, f64, f64) {
    let (rho_prim, xi_vel) = super::euler1disot::cons_to_prim(rho_cons, xi_mom);
    (
        rho_prim,
        xi_vel,
        (gamma - 1.0) * (energy - 0.5 / rho_cons * xi_mom * xi_mom),
    )
}

/// Converts primitive to conservative variables
#[inline(always)]
pub fn prim_to_cons(rho_prim: f64, xi_vel: f64, pressure: f64, gamma: f64) -> (f64, f64, f64) {
    let (rho_cons, xi_mom) = super::euler1disot::prim_to_cons(rho_prim, xi_vel);
    (
        rho_cons,
        xi_mom,
        pressure * (gamma - 1.0).recip() + 0.5 * rho_prim * xi_vel * xi_vel,
    )
}

/// Checks vars for inconsistency, like negative pressure
#[inline(always)]
pub fn validate<const E: usize, const S: usize>(vars: &Variables<E, S>, j_pressure: usize) -> Result<()> {
    ensure!(
        vars.prim.row(j_pressure).fold(true, |acc, x| acc && x > &0.0),
        "Pressure must be positive! Got: {}",
        vars.prim.row(j_pressure)
    );
    Ok(())
}
