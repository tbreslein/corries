// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

use color_eyre::{
    eyre::{bail, ensure},
    Result,
};
use ndarray::{ArrayView1, ArrayView2, Zip};

use crate::{
    boundaryconditions::BoundaryCondition,
    config::physicsconfig::PhysicsConfig,
    data::{Data, DataName, StructAssociation},
    Collectable, Mesh,
};

pub mod systems;
pub mod variables;

pub use systems::*;

/// Trait for Physics objects
pub trait Physics {
    /// construct a new trait object
    fn new(physics_config: &PhysicsConfig) -> Self;

    /// Whether the physics system is adiabatic or not
    fn is_adiabatic(&self) -> bool;

    /// Return copy of the primitive variable in row j at column i
    fn prim_entry(&self, j: usize, i: usize) -> f64;

    /// Return a view of the row j of the primitive variables
    fn prim_row(&self, j: usize) -> ArrayView1<f64>;

    /// Return a view of primitive variables
    fn prim(&self) -> ArrayView2<f64>;

    /// Return copy of the conservative variable in row j at column i
    fn cons_entry(&self, j: usize, i: usize) -> f64;

    /// Return a view of the row j of the conservative variables
    fn cons_row(&self, j: usize) -> ArrayView1<f64>;

    /// Return a view of conservative variables
    fn cons(&self) -> ArrayView2<f64>;

    /// Return a view of eigen value matrix
    fn eigen_vals(&self) -> ArrayView2<f64>;

    /// Return a view of the vector of minimal eigen values
    fn eigen_min(&self) -> ArrayView1<f64>;

    /// Return a view of the vector of maximal eigen values
    fn eigen_max(&self) -> ArrayView1<f64>;

    /// Return a view of speed of sound vector
    fn c_sound(&self) -> ArrayView1<f64>;

    /// Return copy of the physical flux in row j at column i
    fn flux_entry(&self, j: usize, i: usize) -> f64;

    /// Return a view of the row j of the physical flux
    fn flux_row(&self, j: usize) -> ArrayView1<f64>;

    /// Return a view of the physical flux
    fn flux(&self) -> ArrayView2<f64>;

    /// Update primitive variables
    fn update_prim(&mut self);

    /// Update conservative variables
    fn update_cons(&mut self);

    /// Update values not part of the primitive and conservative variables, i.e. speed of sound,
    /// eigen values, physical flux and maybe others
    fn update_derived_values(&mut self);

    /// Assign rhs to primitive variables at row j and column i
    fn assign_prim_entry(&mut self, j: usize, i: usize, rhs: f64);

    /// Assign rhs to primitive variables
    fn assign_prim(&mut self, rhs: &ArrayView2<f64>);

    /// Assign rhs to primitive variables at row j and column i
    fn assign_cons_entry(&mut self, j: usize, i: usize, rhs: f64);

    /// Assign rhs to conservative variables
    fn assign_cons(&mut self, rhs: &ArrayView2<f64>);

    fn assign_c_sound(&mut self, rhs: &ArrayView1<f64>);

    /// Assign the object rhs to this one
    fn assign(&mut self, rhs: &Self) {
        self.assign_prim(&rhs.prim());
        self.assign_cons(&rhs.cons());
        self.assign_c_sound(&rhs.c_sound());
    }

    /// Calculate the CFL limited time step width
    fn calc_dt_cfl<const S: usize>(&self, c_cfl: f64, mesh: &Mesh<S>) -> Result<f64> {
        let dt = c_cfl
            / Zip::from(self.eigen_max())
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

pub fn collect_data<P: Physics + Collectable>(u: &P, data: &mut Data, mesh_offset: usize) -> Result<()> {
    match (data.association, data.name) {
        (StructAssociation::Physics, DataName::Prim(j)) => u.write_vector(&u.prim_row(j), data, mesh_offset),
        (StructAssociation::Physics, DataName::Cons(j)) => u.write_vector(&u.cons_row(j), data, mesh_offset),
        (StructAssociation::Physics, DataName::CSound) => u.write_vector(&u.c_sound(), data, mesh_offset),
        (StructAssociation::Physics, x) => bail!("Tried associating {:?} with Physics!", x),
        (StructAssociation::Mesh, x) | (StructAssociation::TimeStep, x) => {
            bail!("name.association() for {:?} returned {:?}", x, data.association)
        },
    }?;
    return Ok(());
}

#[inline(always)]
pub fn validate<P: Physics>(u: &P) -> Result<()> {
    ensure!(
        u.prim().fold(true, |acc, x| acc && x.is_finite()),
        "Physics::prim must be finite! Got: {}",
        u.prim()
    );
    ensure!(
        u.cons().fold(true, |acc, x| acc && x.is_finite()),
        "Physics::cons must be finite! Got: {}",
        u.cons()
    );
    return Ok(());
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
                u0.prim.row_mut(0).fill(p0);
                u0.prim.row_mut(1).fill(p1);
                u0.prim.row_mut(2).fill(p2);
                let mut u = Euler1DAdiabatic::<S>::new(&PHYSICS_CONFIG);
                u.prim.row_mut(0).fill(p0);
                u.prim.row_mut(1).fill(p1);
                u.prim.row_mut(2).fill(p2);
                println!("before:");
                dbg!(&u.prim);
                dbg!(&u.cons);
                u.update_cons();
                println!("after update_cons:");
                dbg!(&u.prim);
                dbg!(&u.cons);
                u.update_prim();
                println!("after update_prim:");
                dbg!(&u.prim);
                dbg!(&u.cons);
                assert_relative_eq!(u.prim.row(0), u0.prim.row(0), max_relative = 1.0e-12);
                assert_relative_eq!(u.prim.row(1), u0.prim.row(1), max_relative = 1.0e-12);
                assert_relative_eq!(u.prim.row(2), u0.prim.row(2), max_relative = 1.0e-8);
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
                u0.prim.row_mut(0).fill(p0);
                u0.prim.row_mut(1).fill(p1);
                let mut u = Euler1DIsot::<S>::new(&PHYSICS_CONFIG);
                u.prim.row_mut(0).fill(p0);
                u.prim.row_mut(1).fill(p1);
                u.update_cons();
                u.update_prim();
                assert_relative_eq!(u.prim, u0.prim, max_relative = 1.0e-12);
            }
        }
    }
}
