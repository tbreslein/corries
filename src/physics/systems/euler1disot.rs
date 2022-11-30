// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

use ndarray::{Array1, Array2, ArrayView1, ArrayView2};

use crate::physics::Physics;

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

impl<const S: usize> Euler1DIsot<S> {
    /// Constructs a new [Euler1DIsot] object
    pub fn new() -> Self {
        return Self {
            prim: Array2::zeros((2, S)),
            cons: Array2::zeros((2, S)),
            c_sound: Array1::zeros(S),
            eigen_vals: Array2::zeros((2, S)),
        };
    }
}

impl<const S: usize> Physics for Euler1DIsot<S> {
    fn prim_entry(&self, j: usize, i: usize) -> f64 {
        return self.prim[[j, i]];
    }
    fn prim_row(&self, j: usize) -> ArrayView1<f64> {
        return self.prim.row(j);
    }
    fn prim(&self) -> ArrayView2<f64> {
        return self.prim.view();
    }
    fn assign_prim(&mut self, rhs: &Array2<f64>) {
        self.prim.assign(rhs);
    }
    fn cons_entry(&self, j: usize, i: usize) -> f64 {
        return self.cons[[j, i]];
    }
    fn cons_row(&self, j: usize) -> ArrayView1<f64> {
        return self.cons.row(j);
    }
    fn cons(&self) -> ArrayView2<f64> {
        return self.cons.view();
    }
    fn assign_cons(&mut self, rhs: &Array2<f64>) {
        self.cons.assign(rhs);
    }

    fn update_prim(&mut self) {
        self.prim.row_mut(0).assign(&self.cons.row(0));
        self.prim.row_mut(1).assign(&(&self.cons.row(1) / &self.cons.row(0)));
    }
    fn update_cons(&mut self) {
        self.cons.row_mut(0).assign(&self.prim.row(0));
        self.cons.row_mut(1).assign(&(&self.prim.row(1) * &self.prim.row(0)));
    }
}
