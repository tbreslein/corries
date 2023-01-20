// Copyright (c) 2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

use ndarray::{Array1, Array2};

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
}

impl<const E: usize, const S: usize> Variables<E, S> {
    pub fn new() -> Self {
        return Variables {
            prim: Array2::zeros((E, S)),
            cons: Array2::zeros((E, S)),
            c_sound: Array1::zeros(S),
            eigen_vals: Array2::zeros((E, S)),
            flux: Array2::zeros((E, S)),
        };
    }

    pub fn prim_entry(&self, j: usize, i: usize) -> f64;
}

// /// Return copy of the primitive variable in row j at column i
// fn prim_entry(&self, j: usize, i: usize) -> f64;
//
// /// Return a view of the row j of the primitive variables
// fn prim_row(&self, j: usize) -> ArrayView1<f64>;
//
// /// Return a view of primitive variables
// fn prim(&self) -> ArrayView2<f64>;
//
// /// Return copy of the conservative variable in row j at column i
// fn cons_entry(&self, j: usize, i: usize) -> f64;
//
// /// Return a view of the row j of the conservative variables
// fn cons_row(&self, j: usize) -> ArrayView1<f64>;
//
// /// Return a view of conservative variables
// fn cons(&self) -> ArrayView2<f64>;
//
// /// Return a view of eigen value matrix
// fn eigen_vals(&self) -> ArrayView2<f64>;
//
// /// Return a view of the vector of minimal eigen values
// fn eigen_min(&self) -> ArrayView1<f64>;
//
// /// Return a view of the vector of maximal eigen values
// fn eigen_max(&self) -> ArrayView1<f64>;
//
// /// Return a view of speed of sound vector
// fn c_sound(&self) -> ArrayView1<f64>;
//
// /// Return copy of the physical flux in row j at column i
// fn flux_entry(&self, j: usize, i: usize) -> f64;
//
// /// Return a view of the row j of the physical flux
// fn flux_row(&self, j: usize) -> ArrayView1<f64>;
//
// /// Return a view of the physical flux
// fn flux(&self) -> ArrayView2<f64>;
//
