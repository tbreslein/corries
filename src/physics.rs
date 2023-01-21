// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

use color_eyre::{eyre::bail, Result};
use ndarray::{ArrayView1, Zip};

use crate::{boundaryconditions::BoundaryCondition, config::physicsconfig::PhysicsConfig, Mesh};

pub mod systems;
pub mod variables;

pub use systems::*;

/// Trait for Physics objects
pub trait Physics {
    /// Number of equations in the system
    const E: usize;

    /// Number of grid cells in the mesh
    const S: usize;

    /// The exact type of the generic Variables container
    type Vars;

    /// construct a new trait object
    fn new(physics_config: &PhysicsConfig) -> Self;

    /// Returns a reference to the variables at the center of the mesh's cells
    fn cent<'a>(&self) -> &'a Self::Vars;

    /// Returns a reference to the variables at the west border of the mesh's cells
    fn west<'a>(&self) -> &'a Self::Vars;

    /// Returns a reference to the variables at the east border of the mesh's cells
    fn east<'a>(&self) -> &'a Self::Vars;

    /// Whether the physics system is adiabatic or not
    fn is_adiabatic(&self) -> bool;

    /// Update primitive variables
    fn update_prim(&mut self);

    /// Update conservative variables
    fn update_cons(&mut self);

    /// Update values not part of the primitive and conservative variables, i.e. speed of sound,
    /// eigen values, and maybe others
    fn update_derived_values(&mut self);

    fn update_flux_cent(&mut self);

    /// Calculate the CFL limited time step width
    fn calc_dt_cfl<const S: usize>(&self, eigen_max: &ArrayView1<f64>, c_cfl: f64, mesh: &Mesh<S>) -> Result<f64> {
        let dt = c_cfl
            / Zip::from(eigen_max)
                .and(&mesh.cell_width_inv)
                .fold(0.0f64, |acc, eigenval, cw| acc.max(f64::abs(eigenval * cw)));
        if !dt.is_finite() {
            bail!("dt_cfl turned non-finite! Got dt_cfl = {}", dt);
        }
        return Ok(dt);
    }
}

/// Update everything assuming that the conservative variables are up-to-date
pub fn update_everything_from_cons<P: Physics, const S: usize>(
    u: &mut P,
    boundary_west: &mut Box<dyn BoundaryCondition<P, S>>,
    boundary_east: &mut Box<dyn BoundaryCondition<P, S>>,
    mesh: &Mesh<S>,
) {
    u.update_prim();
    update_everything_from_prim(u, boundary_west, boundary_east, mesh);
}

/// Update everything assuming that the primitive variables are up-to-date
pub fn update_everything_from_prim<P: Physics, const S: usize>(
    u: &mut P,
    boundary_west: &mut Box<dyn BoundaryCondition<P, S>>,
    boundary_east: &mut Box<dyn BoundaryCondition<P, S>>,
    mesh: &Mesh<S>,
) {
    boundary_west.apply(u, mesh);
    boundary_east.apply(u, mesh);
    u.update_cons();
    u.update_derived_values();
}

pub fn calc_dt_cfl_generic<const S: usize>(eigen_max: ArrayView1<f64>, c_cfl: f64, mesh: &Mesh<S>) -> Result<f64> {
    let dt = c_cfl
        / Zip::from(eigen_max)
            .and(&mesh.cell_width_inv)
            .fold(0.0f64, |acc, eigenval, cw| acc.max(f64::abs(eigenval * cw)));
    if !dt.is_finite() {
        bail!("dt_cfl turned non-finite! Got dt_cfl = {}", dt);
    }
    return Ok(dt);
}

#[cfg(test)]
mod tests {
    use crate::UnitsMode;

    use super::*;
    use proptest::prelude::*;
    const S: usize = 2;
    const PHYSICS_CONFIG: PhysicsConfig = PhysicsConfig {
        adiabatic_index: 1.4,
        units_mode: UnitsMode::SI,
    };

    mod euler1dadiabatic {
        use super::*;
        use approx::assert_relative_eq;
        proptest! {
            #[test]
            fn conversion(p0 in 0.1f64..10.0, p1 in -10.0f64..10.0, p2 in 0.1f64..10.0) {
                // converting to cons and back to prim should be idempotent
                let mut u0 = Euler1DAdiabatic::<S>::new(&PHYSICS_CONFIG);
                u0.cent().prim.row_mut(0).fill(p0);
                u0.cent().prim.row_mut(1).fill(p1);
                u0.cent().prim.row_mut(2).fill(p2);
                let mut u = Euler1DAdiabatic::<S>::new(&PHYSICS_CONFIG);
                u.cent().prim.row_mut(0).fill(p0);
                u.cent().prim.row_mut(1).fill(p1);
                u.cent().prim.row_mut(2).fill(p2);
                u.update_cons();
                u.update_prim();
                assert_relative_eq!(u.cent().prim.row(0), u0.cent().prim.row(0), max_relative = 1.0e-12);
                assert_relative_eq!(u.cent().prim.row(1), u0.cent().prim.row(1), max_relative = 1.0e-12);
                assert_relative_eq!(u.cent().prim.row(2), u0.cent().prim.row(2), max_relative = 1.0e-8);
            }
        }
    }
    mod euler1disot {
        use super::*;
        use approx::assert_relative_eq;
        proptest! {
            #[test]
            fn conversion(p0 in 0.1f64..100_000.0, p1 in -100_000.0f64..100_000.0) {
                // converting to cons and back to prim should be idempotent
                let mut u0 = Euler1DIsot::<S>::new(&PHYSICS_CONFIG);
                u0.cent().prim.row_mut(0).fill(p0);
                u0.cent().prim.row_mut(1).fill(p1);
                let mut u = Euler1DIsot::<S>::new(&PHYSICS_CONFIG);
                u.cent().prim.row_mut(0).fill(p0);
                u.cent().prim.row_mut(1).fill(p1);
                u.update_cons();
                u.update_prim();
                assert_relative_eq!(u.cent().prim, u0.cent().prim, max_relative = 1.0e-12);
            }
        }
    }
}
