// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

use ndarray::{Array2, ArrayView1, ArrayView2};

mod systems;

/// Trait for Physics objects
pub trait Physics {
    /// Return copy of the primitive variable in equation j at index i
    fn prim_entry(&self, j: usize, i: usize) -> f64;

    /// Return a view of the equation j of the primitive variables
    fn prim_row(&self, j: usize) -> ArrayView1<f64>;

    /// Return a view of primitive variables
    fn prim(&self) -> ArrayView2<f64>;

    /// Assign rhs to primitive variables
    fn assign_prim(&mut self, rhs: &Array2<f64>);

    /// Return copy of the conservative variable in equation j at index i
    fn cons_entry(&self, j: usize, i: usize) -> f64;

    /// Return a view of the equation j of the conservative variables
    fn cons_row(&self, j: usize) -> ArrayView1<f64>;

    /// Return a view of conservative variables
    fn cons(&self) -> ArrayView2<f64>;

    /// Assign rhs to conservative variables
    fn assign_cons(&mut self, rhs: &Array2<f64>);

    /// Update primitive variables
    fn update_prim(&mut self);

    /// Update conservative variables
    fn update_cons(&mut self);
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;
    const S: usize = 10;

    mod euler1disot {
        use crate::physics::systems::euler1disot::Euler1DIsot;

        use super::*;
        use approx::assert_relative_eq;
        proptest! {
            #[test]
            fn conversion(p0 in 0.1f64..100_000.0, p1 in -100_000.0f64..100_000.0) {
                let mut u0 = Euler1DIsot::<S>::new();
                u0.prim.row_mut(0).fill(p0);
                u0.prim.row_mut(1).fill(p1);
                let mut u = Euler1DIsot::<S>::new();
                u.prim.row_mut(0).fill(p0);
                u.prim.row_mut(1).fill(p1);

                // converting to cons and back to prim should be idempotent
                u.update_cons();
                u.update_prim();
                assert_relative_eq!(u.prim, u0.prim, max_relative = 1.0e-12);
            }
        }
    }
}
