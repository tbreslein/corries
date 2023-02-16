// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [NumFlux] trait that identifies structs that can calculate numerical flux.

use crate::{mesh::Mesh, state::Physics, NumFluxConfig, State};
use color_eyre::Result;
use ndarray::{s, Array2};

pub mod hll;
pub mod kt;
pub use self::{hll::Hll, kt::Kt};

/// Trait for structs that can calculate numerical flux
pub trait NumFlux<const E: usize, const S: usize> {
    /// Construct a new [NumFlux] trait object wrapped in a [`color_eyre::Result`]
    ///
    /// # Arguments
    ///
    /// * `numflux_config` - configures the [NumFlux] object
    /// * `mesh` - the [Mesh] this simulation runs on
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    ///
    /// const S: usize = 100;
    /// set_Physics_and_E!(Euler1DIsot);
    /// let mesh: Mesh<S> = Mesh::<S>::new(&MeshConfig::default_riemann_test()).unwrap();
    ///
    /// // construct an Hll numerical flux solver
    /// let hll = Hll::<E, S>::new(&NumFluxConfig::Hll, &mesh);
    ///
    /// // construct a Kt numerical flux solver
    /// let kt = Kt::<E, S>::new(&NumFluxConfig::Kt { limiter_mode: LimiterMode::VanLeer }, &mesh);
    /// ```
    fn new(numflux_config: &NumFluxConfig, mesh: &Mesh<S>) -> Result<Self>
    where
        Self: Sized;

    /// Calculates the numerical derivative along the xi direction.
    ///
    /// # Arguments
    ///
    /// * `dflux_dxi` - derivative of the numerical flux over xi
    /// * `u` - current [Physics] state of the simulation
    /// * `mesh` - the [Mesh] this simulation runs on
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use corries::prelude::*;
    /// use ndarray::Array2;
    ///
    /// const S: usize = 100;
    /// set_Physics_and_E!(Euler1DIsot);
    /// type N = Hll<E,S>;
    ///
    /// // build a configuratoin for riemann tests
    /// let t_end = 0.5;
    /// let folder_name = "results";
    /// let file_name = "noh";
    /// let config = CorriesConfig::default_riemann_test::<N, E, S>(t_end, folder_name, file_name);
    ///
    /// let mesh: Mesh<S> = Mesh::<S>::new(&config.mesh_config).unwrap();
    /// let mut u: State<P, E, S> = State::new(&config.physics_config);
    ///
    /// // construct an Hll numerical flux solver
    /// let mut hll = N::new(&NumFluxConfig::Hll, &mesh).unwrap();
    ///
    /// /* assume that we apply some initial conditions to u */
    ///
    /// // this will store the numerical flux derivative
    /// let mut dflux_dxi = Array2::<f64>::zeros((E, S));
    ///
    /// hll.calc_dflux_dxi(&mut dflux_dxi, &mut u, &mesh).unwrap();
    /// ```
    fn calc_dflux_dxi<P: Physics<E, S>>(
        &mut self,
        dflux_dxi: &mut Array2<f64>,
        u: &mut State<P, E, S>,
        mesh: &Mesh<S>,
    ) -> Result<()>;
}

/// Generic function to calculate the derivative of the numerical flux along the xi direction.
///
/// The purpose of this function is that, after a [NumFlux] object has already calculated the
/// numerical flux, calculating the derivative of that flux is the same for most numerical flux
/// schemes. This function abstracts this derivation.
///
/// # Arguments
///
/// * `dflux_dxi` - numerical flux derivative
/// * `flux_num` - numerical flux
/// * `mesh` - [Mesh] object the simulation runs on
fn calc_dflux_xi_generic<const E: usize, const S: usize>(
    dflux_dxi: &mut Array2<f64>,
    flux_num: &Array2<f64>,
    mesh: &Mesh<S>,
) {
    // PERF: This was benchmarked against using par_azip calls instead of regular assign calls. The
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
