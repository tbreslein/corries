// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [Hll] struct.

use ndarray::{Array2, Array1, s};

use crate::{physics::Physics, mesh::Mesh};

use super::NumFlux;
use super::calc_dflux_xi_generic;

/// Handles calculating numerical flux using the HLL scheme
pub struct Hll<const S: usize, const EQ: usize> {
    /// Physical flux
    flux_phys: Array2<f64>,

    /// Left-side characteristics
    sl: Array1<f64>,

    /// Right-side characteristics
    sr: Array1<f64>,

    /// Helper array containting `1.0 / (sr - sl)`
    inv_sr_minus_sl: Array1<f64>,

    /// Helper array containting `sr * sl`
    sr_times_sl: Array1<f64>,

    /// Numerical flux
    flux_num: Array2<f64>,
}

impl<const S: usize, const EQ: usize> Hll<S, EQ> {
    /// Constructs a new [Hll] object
    pub fn new() -> Self {
        return Self {flux_phys: Array2::zeros((EQ, S)), sl: Array1::zeros(S), sr: Array1::zeros(S), inv_sr_minus_sl: Array1::zeros(S), sr_times_sl: Array1::zeros(S), flux_num: Array2::zeros((EQ, S))};
    }
}

impl<const S: usize, const EQ: usize> NumFlux<S, EQ> for Hll<S, EQ> {
    fn calc_dflux_dxi(&mut self, dflux_dxi: &mut Array2<f64>, u: &mut Physics<S, EQ>, mesh: &Mesh<S>) {
        u.calc_physical_flux_euler1d_adiabatic(&mut self.flux_phys);
        u.update_eigen_vals();
        let slice = s![mesh.ixi_in-1..=mesh.ixi_out];
        let slice_p1 = s![mesh.ixi_in..=mesh.ixi_out+1];

        for i in mesh.ixi_in-1..=mesh.ixi_out {
            self.sl[i] = 0.0f64.min(u.eigen_vals[[0, i]].min(u.eigen_vals[[0, i+1]]));
            self.sr[i] = 0.0f64.max(u.eigen_vals[[EQ-1, i]].max(u.eigen_vals[[EQ-1, i+1]]));
        }

        // TODO: Benchmark this version against azip, Zip and raw for loops
        self.inv_sr_minus_sl.slice_mut(slice).assign(&(1.0 / (&self.sr.slice(slice) - &self.sl.slice(slice))));
        self.sr_times_sl.slice_mut(slice).assign(&(&self.sr.slice(slice) * &self.sl.slice(slice)));

        // TODO: Benchmark this version against azip, Zip and raw for loops
        for j in 0..EQ {
            // TODO: Benchmark whether these explicit views create overhead
            let mut flux_num_j = self.flux_num.row_mut(j);
            let flux_phys_j = self.flux_phys.row(j);
            let uc_j = u.cons.row(j);
            let a = self.inv_sr_minus_sl.slice(slice);
            let b = self.sr_times_sl.slice(slice);
            let sl = self.sl.slice(slice);
            let sr = self.sr.slice(slice);

            flux_num_j.slice_mut(slice).assign(&(&a * (&sr * &flux_phys_j.slice(slice) - &sl * &flux_phys_j.slice(slice_p1) + &b * (&uc_j.slice(slice_p1) - &uc_j.slice(slice)))));
        }
        calc_dflux_xi_generic(dflux_dxi, &self.flux_num, &mesh);
    }
}
