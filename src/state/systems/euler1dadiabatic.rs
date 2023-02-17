// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [Euler1DAdiabatic] struct, which is an implementer for [Physics] for 1-dimensional
//! adiabatic Euler equations

use crate::{state::Physics, variables::Variables};
use color_eyre::{eyre::ensure, Result};

const E: usize = 3;

/// Implenter for [Physics] that manipulates [State](crate::state::State) for 1-dimensional
/// adiabatic Euler equations.
///
/// # Variables
///
/// The equations are sorted such that they correspond to the primitive and conservative variables
/// such that:
///
/// equation index | primitive variable | conservative variable
/// ---            | ---                | ---
/// 0              | mass density       | mass density
/// 1              | xi velocity        | xi momentum
/// 2              | pressure           | inner energy
///
/// Keep in mind that the coordinates in [corries](crate) are generalised, and that it is a purely
/// 1-dimensional framework. `xi` is the only coordinate that can have more than one cell, whereas
/// the `eta` and `Phi` directions only ever have one cell.
///
/// [Euler1DAdiabatic] does not model velocities along the `eta` direction.
/// Accessors to any of these will return a vector of zeroes, whereas the corresponding equation
/// indexes will return [usize::MAX].
///
/// # Conversions
///
/// Let
///
/// * `up_j`: primitive variables at equation index `j`
/// * `uc_j`: conservative variables at equation index `j`
/// * `Fp_j`: physical flux at equation index `j`
/// * `gamma`: the adiabatic index
///
/// In combination with the table up to this means that, for example, `up_0` corresponds the
/// primitive mass density, and `uc_1` corresponds to the momentum along the xi direction.
///
/// Assuming we start with the primitive variables, the conservative variables are calculated with:
///
/// ```text
/// uc_0 = up_0
/// uc_1 = up_0 * up_1
/// uc_2 = up_2 / (gamma - 1.0) + 0.5 * up_0 * up_1 * up_1
/// ```
///
/// They are converted back with:
///
/// ```text
/// up_0 = uc_0
/// up_1 = uc_1 / uc_0
/// up_2 = (gamma - 1.0) * (uc_2 - 0.5 / uc_0 * uc_1 * uc_1)
/// ```
///
/// The physical flux is calculated with:
///
/// ```text
/// Fp_0 = uc_1
/// Fp_1 = up_1 * uc_1 + up_2
/// Fp_2 = (uc_2 + up_2) * up_1
/// ```
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Euler1DAdiabatic<const S: usize>;

unsafe impl<const S: usize> Send for Euler1DAdiabatic<S> {}
unsafe impl<const S: usize> Sync for Euler1DAdiabatic<S> {}

impl<const S: usize> Physics<E, S> for Euler1DAdiabatic<S> {
    const NUM_EQ: usize = E;
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
        for i in 0..S {
            (
                vars.flux[[Self::JRHO, i]],
                vars.flux[[Self::JXI, i]],
                vars.flux[[Self::JPRESSURE, i]],
            ) = calc_flux(
                vars.prim[[Self::JXI, i]],
                vars.cons[[Self::JXI, i]],
                vars.prim[[Self::JPRESSURE, i]],
                vars.cons[[Self::JPRESSURE, i]],
            );
        }
    }

    #[inline(always)]
    fn validate(vars: &Variables<E, S>) -> Result<()> {
        super::euler1disot::validate(vars, Self::JRHO)?;
        validate(vars, Self::JPRESSURE)
    }
}

/// Updates physical flux
#[inline(always)]
pub fn calc_flux(xi_vel: f64, xi_mom: f64, pressure: f64, energy: f64) -> (f64, f64, f64) {
    (xi_mom, xi_vel * xi_mom + pressure, (energy + pressure) * xi_vel)
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
    check_positive_arrayd!(vars.prim.row(j_pressure));
    Ok(())
}
