// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [NumFlux] trait that identifies structs that can calculate numerical flux, and the
//! [init_numflux] function that handles constructing [NumFlux] objects.

use color_eyre::Result;
use ndarray::{s, Array2};

use crate::{
    config::numericsconfig::{NumFluxMode, NumericsConfig},
    mesh::Mesh,
    physics::Physics,
};

use self::hll::Hll;

mod hll;

/// Trait for structs that can calculate numerical flux
pub trait NumFlux<const S: usize, const EQ: usize> {
    /// Calculates the numerical derivative along the xi direction
    fn calc_dflux_dxi(&mut self, dflux_dxi: &mut Array2<f64>, u: &mut Physics<S, EQ>, mesh: &Mesh<S>) -> Result<()>;
}

/// Constructs an `impl Numflux<S, EQ>`.
pub fn init_numflux<const S: usize, const EQ: usize>(numericsconf: &NumericsConfig) -> impl NumFlux<S, EQ> {
    return match numericsconf.numflux_mode {
        NumFluxMode::Hll => Hll::new(),
    };
}

/// Generic function to calculate the derivative of the numerical flux along the xi direction.
///
/// The purpose of this function is that, after a [NumFlux] object has already calculated the
/// numerical flux, calculating the derivative of that flux is the same for most numerical flux
/// schemes.
///
/// # Arguments
///
/// * `dflux_dxi` - numerical flux derivative
/// * `flux_num` - numerical flux
/// * `mesh` - [Mesh] object
fn calc_dflux_xi_generic<const S: usize>(dflux_dxi: &mut Array2<f64>, flux_num: &Array2<f64>, mesh: &Mesh<S>) {
    for j in 0..dflux_dxi.nrows() {
        let s = s![mesh.ixi_in..=mesh.ixi_out];
        let s_m1 = s![mesh.ixi_in - 1..=mesh.ixi_out - 1];
        dflux_dxi
            .row_mut(j)
            .slice_mut(s)
            .assign(&(&mesh.deta_dphi_d_volume.slice(s) * (&flux_num.row(j).slice(s) - &flux_num.row(j).slice(s_m1))));
    }
}
