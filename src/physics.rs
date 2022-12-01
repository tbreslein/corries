// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

use ndarray::{Array2, ArrayView1, ArrayView2};

use crate::config::physicsconfig::PhysicsConfig;

pub mod systems;

pub use systems::*;

/// Trait for Physics objects
pub trait Physics {
    /// construct a new trait object
    fn new(physics_config: &PhysicsConfig) -> Self;

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
    fn assign_prim(&mut self, rhs: &Array2<f64>);

    /// Assign rhs to primitive variables at row j and column i
    fn assign_cons_entry(&mut self, j: usize, i: usize, rhs: f64);

    /// Assign rhs to conservative variables
    fn assign_cons(&mut self, rhs: &Array2<f64>);
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;
    const S: usize = 10;
    const PHYSICS_CONFIG: PhysicsConfig = PhysicsConfig { adiabatic_index: 1.4 };

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
                u.update_cons();
                u.update_prim();
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
