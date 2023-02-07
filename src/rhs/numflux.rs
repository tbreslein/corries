// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [NumFlux] trait that identifies structs that can calculate numerical flux, and the
//! [init_numflux] function that handles constructing [NumFlux] objects.

use crate::{mesh::Mesh, state::Physics, NumFluxConfig, State};
use color_eyre::Result;
use ndarray::{s, Array2};

pub mod hll;
pub mod kt;
pub use self::{hll::Hll, kt::Kt};

/// Trait for structs that can calculate numerical flux
pub trait NumFlux<const E: usize, const S: usize> {
    /// construct a new trait object
    fn new(numflux_config: &NumFluxConfig, mesh: &Mesh<S>) -> Result<Self>
    where
        Self: Sized;

    /// Calculates the numerical derivative along the xi direction
    fn calc_dflux_dxi<P: Physics<E, S>>(
        &mut self,
        dflux_dxi: &mut Array2<f64>,
        u: &mut State<P, E, S>,
        mesh: &Mesh<S>,
    ) -> Result<()>;
}

/// Constructs an `impl Numflux<S, EQ>`.
pub fn init_numflux<N: NumFlux<E, S>, const E: usize, const S: usize>(
    numflux_config: &NumFluxConfig,
    mesh: &Mesh<S>,
) -> Result<N> {
    N::new(numflux_config, mesh)
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
fn calc_dflux_xi_generic<const E: usize, const S: usize>(
    dflux_dxi: &mut Array2<f64>,
    flux_num: &Array2<f64>,
    mesh: &Mesh<S>,
) {
    // NOTE: This was benchmarked against using par_azip calls instead of regular assign calls. The
    // performance differences were negligable.
    let s = s![mesh.ixi_in..=mesh.ixi_out];
    let sm1 = s![mesh.ixi_in - 1..=mesh.ixi_out - 1];
    for j in 0..E {
        dflux_dxi
            .row_mut(j)
            .slice_mut(s)
            .assign(&(&mesh.deta_dphi_d_volume.slice(s) * (&flux_num.row(j).slice(s) - &flux_num.row(j).slice(sm1))));
    }
}
