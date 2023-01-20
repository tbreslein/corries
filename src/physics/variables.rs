// Copyright (c) 2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

use ndarray::{Array1, Array2, ArrayView1, ArrayView2};

use crate::PhysicsConfig;

/// Contains the variables of a physics system
#[derive(Debug)]
pub struct Variables<const E: usize, const S: usize> {
    /// Primitive variables
    pub prim: Array2<f64>,

    /// Conservative variables
    pub cons: Array2<f64>,

    /// Speed of sound
    pub c_sound: Array1<f64>,

    /// Eigen values
    pub eigen_vals: Array2<f64>,

    /// Physical flux
    pub flux: Array2<f64>,

    /// Adiabatic index
    pub gamma: f64,
}

impl<const E: usize, const S: usize> Variables<E, S> {
    pub fn new(physics_config: &PhysicsConfig) -> Self {
        return Variables {
            prim: Array2::zeros((E, S)),
            cons: Array2::zeros((E, S)),
            c_sound: Array1::zeros(S),
            eigen_vals: Array2::zeros((E, S)),
            flux: Array2::zeros((E, S)),
            gamma: physics_config.adiabatic_index,
        };
    }

    /// Return copy of the primitive variable in row j at column i
    #[inline(always)]
    pub fn prim_entry(&self, j: usize, i: usize) -> f64 {
        return self.prim[[j, i]];
    }

    /// Return a view of the row j of the primitive variables
    #[inline(always)]
    pub fn prim_row(&self, j: usize) -> ArrayView1<f64> {
        return self.prim.row(j);
    }

    /// Return a view of primitive variables
    #[inline(always)]
    pub fn prim(&self) -> ArrayView2<f64> {
        return self.prim.view();
    }

    /// Return copy of the conservative variable in row j at column i
    #[inline(always)]
    pub fn cons_entry(&self, j: usize, i: usize) -> f64 {
        return self.cons[[j, i]];
    }

    /// Return a view of the row j of the conservative variables
    #[inline(always)]
    pub fn cons_row(&self, j: usize) -> ArrayView1<f64> {
        return self.cons.row(j);
    }

    /// Return a view of conservative variables
    #[inline(always)]
    pub fn cons(&self) -> ArrayView2<f64> {
        return self.cons.view();
    }

    /// Return a view of eigen value matrix
    #[inline(always)]
    pub fn eigen_vals(&self) -> ArrayView2<f64> {
        return self.eigen_vals.view();
    }

    /// Return a view of the vector of minimal eigen values
    #[inline(always)]
    pub fn eigen_min(&self) -> ArrayView1<f64> {
        return self.eigen_vals.row(0);
    }

    /// Return a view of the vector of maximal eigen values
    #[inline(always)]
    pub fn eigen_max(&self) -> ArrayView1<f64> {
        return self.eigen_vals.row(E-1);
    }

    /// Return a view of speed of sound vector
    #[inline(always)]
    pub fn c_sound(&self) -> ArrayView1<f64> {
        return self.c_sound.view();
    }

    /// Return copy of the physical flux in row j at column i
    #[inline(always)]
    pub fn flux_entry(&self, j: usize, i: usize) -> f64 {
        return self.flux[[j, i]];
    }

    /// Return a view of the row j of the physical flux
    #[inline(always)]
    pub fn flux_row(&self, j: usize) -> ArrayView1<f64> {
        return self.flux.row(j);
    }

    /// Return a view of the physical flux
    #[inline(always)]
    pub fn flux(&self) -> ArrayView2<f64> {
        return self.flux.view();
    }
}
