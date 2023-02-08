// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [Physics] trait which defines the behavior embedded in objects

use color_eyre::{eyre::bail, Result};
use ndarray::{ArrayView1, Zip};

use crate::{boundaryconditions::BoundaryCondition, variables::Variables, Mesh};

/// Trait for Physics objects
pub trait Physics<const E: usize, const S: usize> {
    /// Whether the implementer of this trait is adiabatic
    const IS_ADIABATIC: bool;

    /// Whether the implementer of this trait is isothermal
    const IS_ISOTHERMAL: bool = !Self::IS_ADIABATIC;

    /// The row index for the mass density in Variables::prim or Variables::cons; might be set to
    /// usize::MAX if the implementer does not carry the density in either of those arrays.
    const JRHO: usize;

    /// The row index for the xi velocity in Variables::prim or Variables::cons; might be set to
    /// usize::MAX if the implementer does not carry the density in either of those arrays.
    const JXI: usize;

    /// The row index for the eta velocity in Variables::prim or Variables::cons; might be set to
    /// usize::MAX if the implementer does not carry the density in either of those arrays.
    const JETA: usize;

    /// The row index for the pressure in Variables::prim or Variables::cons; might be set to
    /// usize::MAX if the implementer does not carry the density in either of those arrays.
    const JPRESSURE: usize;

    /// Construct a new Physics object
    fn new() -> Self;

    /// todo
    fn rho_prim(vars: &Variables<E, S>) -> ArrayView1<f64> {
        if Self::JRHO < usize::MAX {
            vars.prim.row(Self::JRHO)
        } else {
            vars.zero_vec()
        }
    }

    /// todo
    fn rho_cons(vars: &Variables<E, S>) -> ArrayView1<f64> {
        if Self::JRHO < usize::MAX {
            vars.cons.row(Self::JRHO)
        } else {
            vars.zero_vec()
        }
    }

    /// todo
    fn xi_vel(vars: &Variables<E, S>) -> ArrayView1<f64> {
        if Self::JXI < usize::MAX {
            vars.prim.row(Self::JXI)
        } else {
            vars.zero_vec()
        }
    }

    /// todo
    fn xi_mom(vars: &Variables<E, S>) -> ArrayView1<f64> {
        if Self::JXI < usize::MAX {
            vars.cons.row(Self::JXI)
        } else {
            vars.zero_vec()
        }
    }

    /// todo
    fn eta_vel(vars: &Variables<E, S>) -> ArrayView1<f64> {
        if Self::JETA < usize::MAX {
            vars.prim.row(Self::JETA)
        } else {
            vars.zero_vec()
        }
    }

    /// todo
    fn eta_mom(vars: &Variables<E, S>) -> ArrayView1<f64> {
        if Self::JETA < usize::MAX {
            vars.cons.row(Self::JETA)
        } else {
            vars.zero_vec()
        }
    }

    /// todo
    fn pressure(vars: &Variables<E, S>) -> ArrayView1<f64> {
        if Self::JPRESSURE < usize::MAX {
            vars.prim.row(Self::JPRESSURE)
        } else {
            vars.zero_vec()
        }
    }

    /// todo
    fn energy(vars: &Variables<E, S>) -> ArrayView1<f64> {
        if Self::JPRESSURE < usize::MAX {
            vars.cons.row(Self::JPRESSURE)
        } else {
            vars.zero_vec()
        }
    }

    /// Updates primitive variables in vars.prim
    fn update_prim(vars: &mut Variables<E, S>);

    /// Updates conservative variables vars.cons
    fn update_cons(vars: &mut Variables<E, S>);

    /// Update values not part of the primitive and conservative variables, i.e. speed of sound,
    /// eigen values, and maybe others
    fn update_derived_variables(vars: &mut Variables<E, S>) {
        if Self::IS_ADIABATIC {
            Self::update_c_sound(vars);
        }
        Self::update_eigen_vals(vars);
    }

    /// Updates the speed of sound in vars.c_sound
    fn update_c_sound(vars: &mut Variables<E, S>) {
        if Self::IS_ADIABATIC {
            // PERF: This was benchmarked with the following alternatives:
            // - raw index loop (~10% faster than the alternatives)
            // - ndarray's `assign`
            // - ndarray's `Zip`
            for i in 0..S {
                vars.c_sound[i] = (vars.gamma * Self::pressure(vars)[i] / Self::rho_prim(vars)[i]).sqrt();
            }
        }
        // isothermal case is a no-op
    }

    /// Updates the eigen values in vars.eigen_vals
    fn update_eigen_vals(vars: &mut Variables<E, S>) {
        Self::update_eigen_vals_min(vars);
        for j in 1..E - 1 {
            vars.eigen_vals.row_mut(j).assign(&vars.prim.row(Self::JXI));
        }
        Self::update_eigen_vals_max(vars);
    }

    /// Updates the minimal eigen values, i.e. the row in vars.eigen_vals representing the smallest
    /// eigen value per cell
    fn update_eigen_vals_min(vars: &mut Variables<E, S>) {
        vars.eigen_vals
            .row_mut(0)
            .assign(&(&vars.prim.row(Self::JXI) - &vars.c_sound));
    }

    /// Updates the maximal eigen values, i.e. the row in vars.eigen_vals representing the largest
    /// eigen value per cell
    fn update_eigen_vals_max(vars: &mut Variables<E, S>) {
        vars.eigen_vals
            .row_mut(E - 1)
            .assign(&(&vars.prim.row(Self::JXI) + &vars.c_sound));
    }

    /// Updates everything in vars, assuming vars.prim is already up-to-date
    fn update_vars_from_prim(
        vars: &mut Variables<E, S>,
        boundary_west: &mut Box<dyn BoundaryCondition<E, S>>,
        boundary_east: &mut Box<dyn BoundaryCondition<E, S>>,
        mesh: &Mesh<S>,
    ) {
        boundary_west.apply(vars, mesh);
        boundary_east.apply(vars, mesh);
        Self::update_cons(vars);
        Self::update_derived_variables(vars);
    }

    /// Updates everything in vars, assuming vars.cons is already up-to-date
    fn update_vars_from_cons(
        vars: &mut Variables<E, S>,
        boundary_west: &mut Box<dyn BoundaryCondition<E, S>>,
        boundary_east: &mut Box<dyn BoundaryCondition<E, S>>,
        mesh: &Mesh<S>,
    ) {
        Self::update_prim(vars);
        Self::update_vars_from_prim(vars, boundary_west, boundary_east, mesh);
    }

    /// Updates the physical flux in vars.flux
    fn update_flux(vars: &mut Variables<E, S>);

    /// Calculate the CFL limited time step width
    fn calc_dt_cfl(eigen_max: &ArrayView1<f64>, c_cfl: f64, mesh: &Mesh<S>) -> Result<f64> {
        let dt = c_cfl
            / Zip::from(eigen_max)
                .and(&mesh.cell_width_inv)
                .fold(0.0f64, |acc, eigenval, cw| acc.max(f64::abs(eigenval * cw)));
        if !dt.is_finite() {
            bail!("dt_cfl turned non-finite! Got dt_cfl = {}", dt);
        }
        Ok(dt)
    }

    /// Validates the contents of vars, making sure it does not have non-finite numbers for example
    fn validate(vars: &Variables<E, S>) -> Result<()>;
}
