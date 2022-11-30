// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [Physics] struct that handles the variables and physical state of the simulation.

use ndarray::{Array1, Array2, ArrayView1, ArrayView2};

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

/// Test struct for using a trait for Physics
#[derive(Debug)]
pub struct Euler1DIsot<const S: usize> {
    /// Primitive variables
    pub prim: Array2<f64>,

    /// Conservative variables
    pub cons: Array2<f64>,

    /// Speed of sound
    pub c_sound: Array1<f64>,

    /// Eigen values
    pub eigen_vals: Array2<f64>,
}

impl<const S: usize> Physics for Euler1DIsot<S> {
    fn prim_entry(&self, j: usize, i: usize) -> f64 { return self.prim[[j, i]]; }
    fn prim_row(&self, j: usize) -> ArrayView1<f64> { return self.prim.row(j); }
    fn prim(&self) -> ArrayView2<f64> { return self.prim.view(); }
    fn assign_prim(&mut self, rhs: &Array2<f64>) { self.prim.assign(rhs); }
    fn cons_entry(&self, j: usize, i: usize) -> f64 { return self.cons[[j, i]]; }
    fn cons_row(&self, j: usize) -> ArrayView1<f64> { return self.cons.row(j); }
    fn cons(&self) -> ArrayView2<f64> { return self.cons.view(); }
    fn assign_cons(&mut self, rhs: &Array2<f64>) { self.cons.assign(rhs); }

    fn update_prim(&mut self) {
        self.prim.row_mut(0).assign(&self.cons.row(0));
        self.prim.row_mut(1).assign(&(&self.cons.row(1) / &self.cons.row(0)));
    }
    fn update_cons(&mut self) {
        self.cons.row_mut(0).assign(&self.prim.row(0));
        self.cons.row_mut(1).assign(&(&self.prim.row(1) * &self.prim.row(0)));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;
    const S: usize = 10;

    mod euler1disot {
        use super::*;
        use approx::assert_relative_eq;
        proptest! {
            #[test]
            fn conversion(p0 in 0.1f64..100_000.0, p1 in -100_000.0f64..100_000.0) {
                let mut u0: Euler1DIsot<S> = Euler1DIsot {prim: Array2::zeros((2, S)), cons: Array2::zeros((2, S))};
                u0.prim.row_mut(0).fill(p0);
                u0.prim.row_mut(1).fill(p1);
                let mut u: Euler1DIsot<S> = Euler1DIsot {prim: Array2::zeros((2, S)), cons: Array2::zeros((2, S))};
                u.prim.row_mut(0).fill(p0);
                u.prim.row_mut(1).fill(p1);
                // let mut u0: Physics<S, EQ> = Physics::new(&physicsconf);
                // u0.prim.row_mut(0).fill(p0);
                // u0.prim.row_mut(1).fill(p1);
                // let mut u: Physics<S, EQ> = Physics::new(&physicsconf);
                // u.prim.row_mut(0).fill(p0);
                // u.prim.row_mut(1).fill(p1);

                // converting to cons and back to prim should be idempotent
                u.update_cons();
                u.update_prim();
                assert_relative_eq!(u.prim, u0.prim, max_relative = 1.0e-12);
            }
        }
    }
}
