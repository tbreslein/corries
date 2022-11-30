// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

use ndarray::{Array1, Array2, ArrayView1, ArrayView2, Zip, ArrayViewMut2};

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
    #[inline(always)]
    fn prim_entry(&self, j: usize, i: usize) -> f64 {
        return self.prim[[j, i]];
    }

    #[inline(always)]
    fn prim_row(&self, j: usize) -> ArrayView1<f64> {
        return self.prim.row(j);
    }

    #[inline(always)]
    fn prim(&self) -> ArrayView2<f64> {
        return self.prim.view();
    }

    #[inline(always)]
    fn assign_prim(&mut self, rhs: &Array2<f64>) {
        self.prim.assign(rhs);
    }

    #[inline(always)]
    fn cons_entry(&self, j: usize, i: usize) -> f64 {
        return self.cons[[j, i]];
    }

    #[inline(always)]
    fn cons_row(&self, j: usize) -> ArrayView1<f64> {
        return self.cons.row(j);
    }

    #[inline(always)]
    fn cons(&self) -> ArrayView2<f64> {
        return self.cons.view();
    }

    #[inline(always)]
    fn assign_cons(&mut self, rhs: &Array2<f64>) {
        self.cons.assign(rhs);
    }

    fn update_prim(&mut self) {
        cons_to_prim(self.prim.view_mut(), self.cons.view());
    }
    fn update_cons(&mut self) {
        prim_to_cons(self.cons.view_mut(), self.prim.view());
    }
}

#[inline(always)]
pub fn cons_to_prim(mut prim: ArrayViewMut2<f64>, cons: ArrayView2<f64>) {
        Zip::from(prim.row_mut(0))
            .and(cons.row(0))
            .for_each(|rho_prim, &rho_cons| {
                *rho_prim = rho_cons;
            });

        Zip::from(prim.row_mut(1))
            .and(cons.row(0))
            .and(cons.row(1))
            .for_each(|xi_vel, &xi_mom, &rho_cons| {
                *xi_vel = xi_mom / rho_cons;
            });
}

#[inline(always)]
pub fn prim_to_cons(mut cons: ArrayViewMut2<f64>, prim: ArrayView2<f64>) {
        Zip::from(cons.row_mut(0))
            .and(prim.row(0))
            .for_each(|rho_cons, &rho_prim| {
                *rho_cons = rho_prim;
            });

        Zip::from(cons.row_mut(1))
            .and(prim.row(0))
            .and(prim.row(1))
            .for_each(|xi_mom, &xi_vel, &rho_cons| {
                *xi_mom = xi_vel / rho_cons;
            });
}
