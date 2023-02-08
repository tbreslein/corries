// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [Euler1DIsot] struct, which is an implementer for [Physics] for 1-dimensional
//! isothermal Euler equations

use color_eyre::{eyre::ensure, Result};
use ndarray::{ArrayView1, ArrayView2, ArrayViewMut2, Zip};

use crate::{state::Physics, variables::Variables};

const E: usize = 2;

/// Test struct for using a trait for Physics
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Euler1DIsot<const S: usize>;

unsafe impl<const S: usize> Send for Euler1DIsot<S> {}
unsafe impl<const S: usize> Sync for Euler1DIsot<S> {}

impl<const S: usize> Physics<E, S> for Euler1DIsot<S> {
    const IS_ADIABATIC: bool = false;
    const JRHO: usize = 0;
    const JXI: usize = 1;
    const JETA: usize = std::usize::MAX;
    const JPRESSURE: usize = std::usize::MAX;

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
            (vars.prim[[Self::JRHO, i]], vars.prim[[Self::JXI, i]]) =
                cons_to_prim(vars.cons[[Self::JRHO, i]], vars.cons[[Self::JXI, i]]);
        }
    }

    #[inline(always)]
    fn update_cons(vars: &mut Variables<E, S>) {
        // PERF: This was benchmarked between the following options:
        // - raw index loop over columns calculating the tuple of the primitive values per column
        // (~20% speed up over the zip)
        // - ndarray::Zip!, where each row is calculated in a separate zip
        for i in 0..S {
            (vars.cons[[Self::JRHO, i]], vars.cons[[Self::JXI, i]]) =
                prim_to_cons(vars.prim[[Self::JRHO, i]], vars.prim[[Self::JXI, i]]);
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
            vars.c_sound.view(),
        );
    }

    #[inline(always)]
    fn validate(vars: &Variables<E, S>) -> Result<()> {
        validate(vars, Self::JRHO)
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
pub fn cons_to_prim(rho_cons: f64, xi_mom: f64) -> (f64, f64) {
    (rho_cons, xi_mom / rho_cons)
}

/// Converts primitive to conservative variables
#[inline(always)]
pub fn prim_to_cons(rho_prim: f64, xi_vel: f64) -> (f64, f64) {
    (rho_prim, xi_vel * rho_prim)
}

/// Checks vars for inconsistency, like mass density
#[inline(always)]
pub fn validate<const E: usize, const S: usize>(vars: &Variables<E, S>, j_rho: usize) -> Result<()> {
    ensure!(
        vars.prim.row(j_rho).fold(true, |acc, x| acc && x > &0.0),
        "Mass density must be positive! Got: {}",
        vars.prim.row(j_rho)
    );
    Ok(())
}
