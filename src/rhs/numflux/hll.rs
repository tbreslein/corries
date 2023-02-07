// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [Hll] struct.

use color_eyre::eyre::{bail, ensure, Context};
use color_eyre::Result;
use ndarray::{par_azip, s, Array1, Array2};

use crate::errorhandling::Validation;
use crate::{mesh::Mesh, state::Physics};
use crate::{NumFluxConfig, State};

use super::calc_dflux_xi_generic;
use super::NumFlux;

/// Handles calculating numerical flux using the HLL scheme with 0-order reconstruction
pub struct Hll<const E: usize, const S: usize> {
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

impl<const E: usize, const S: usize> NumFlux<E, S> for Hll<E, S> {
    fn new(numflux_config: &NumFluxConfig, _: &Mesh<S>) -> Result<Self> {
        match numflux_config {
            NumFluxConfig::Hll => Ok(Self {
                sl: Array1::zeros(S),
                sr: Array1::zeros(S),
                inv_sr_minus_sl: Array1::zeros(S),
                sr_times_sl: Array1::zeros(S),
                flux_num: Array2::zeros((E, S)),
            }),
            _ => bail!("Tried constructing Hll, but numflux_config does not contain NumFluxConfig::Hll!"),
        }
    }

    fn calc_dflux_dxi<P: Physics<E, S>>(
        &mut self,
        dflux_dxi: &mut Array2<f64>,
        u: &mut State<P, E, S>,
        mesh: &Mesh<S>,
    ) -> Result<()> {
        // NOTE: Assumes that u.eigen_vals are already up to date
        u.update_flux();

        let uflux = &u.cent.flux;
        let ucons = &u.cent.cons;
        let eigen_min = &u.cent.eigen_min();
        let eigen_max = &u.cent.eigen_max();
        let s = s![(mesh.ixi_in - 1)..=mesh.ixi_out];
        let sp1 = s![mesh.ixi_in..=(mesh.ixi_out + 1)];

        // PERF: This was benchmarked in comparison with running these 4 calls in blocks of two
        // parallel calls using rayon::join. The operations themselves, even on a 500 cell mesh,
        // were fast enough to outpace the overhead of the joins. Calling the parallel zip as 4
        // sequential calls was the clear benchmark winner on two seperate machines.
        par_azip!((
                sl in &mut self.sl.slice_mut(s), &ev1 in &eigen_min.slice(s), &ev2 in &eigen_min.slice(sp1))
                *sl = 0.0f64.min(ev1.min(ev2)));
        par_azip!((
                sr in &mut self.sr.slice_mut(s), &ev1 in &eigen_max.slice(s), &ev2 in &eigen_max.slice(sp1))
                *sr = 0.0f64.max(ev1.max(ev2)));
        par_azip!((
                a in &mut self.inv_sr_minus_sl.slice_mut(s), &sr in &self.sr.slice(s), &sl in &self.sl.slice(s))
                *a = 1.0 / (sr - sl));
        par_azip!((
                b in &mut self.sr_times_sl.slice_mut(s), &sr in &self.sr.slice(s), &sl in &self.sl.slice(s))
                *b = sr * sl);

        // PERF: This loop was benchmarked against using a parallel iterator on self.flux_num with
        // the benchmark using 3 equations and a mesh of 500 cells. The raw loop won over the
        // parallel iterator in that case, probably due to parallelisation in the .assign(), since
        // that uses parallel BLAS under the hood.
        // Adding a parallel iter on top of that creates too much overhead to be useful.
        for j in 0..E {
            self.flux_num.row_mut(j).slice_mut(s).assign(
                &(&self.inv_sr_minus_sl.slice(s)
                    * (&self.sr.slice(s) * &uflux.row(j).slice(s) - &self.sl.slice(s) * &uflux.row(j).slice(sp1)
                        + &self.sr_times_sl.slice(s) * (&ucons.row(j).slice(sp1) - &ucons.row(j).slice(s)))),
            );
        }

        calc_dflux_xi_generic::<E, S>(dflux_dxi, &self.flux_num, mesh);
        if cfg!(feature = "validation") {
            self.validate()
                .context("Calling Hll::validate in Hll::calc_dflux_dxi")?;
        }
        Ok(())
    }
}

impl<const E: usize, const S: usize> Validation for Hll<E, S> {
    fn validate(&self) -> Result<()> {
        ensure!(
            self.sl.fold(true, |acc, x| acc && x.is_finite()),
            "Hll::sl must be finite! Got: {}",
            self.sl
        );
        ensure!(
            self.sr.fold(true, |acc, x| acc && x.is_finite()),
            "Hll::sr must be finite! Got: {}",
            self.sr
        );
        ensure!(
            self.inv_sr_minus_sl.fold(true, |acc, x| acc && x.is_finite()),
            "Hll::inv_sr_minus_sl must be finite! Got:\ninv_sr_minus_sl = {}\nsr = {}\nsl = {}",
            self.inv_sr_minus_sl,
            self.sr,
            self.sl
        );
        ensure!(
            self.flux_num.fold(true, |acc, x| acc && x.is_finite()),
            "Hll::flux_num must be finite! Got: {}",
            self.flux_num
        );
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::prelude::*;

    use super::*;
    use approx::assert_relative_eq;
    use ndarray::Array2;
    const S: usize = 8;
    const EQ: usize = 2;
    const MESHCONFIG: MeshConfig = MeshConfig {
        mode: MeshMode::Cartesian,
        xi_in: 2.0,
        xi_out: 3.0,
        ratio_disk: 1.0,
    };
    const PHYSICSCONFIG: PhysicsConfig = PhysicsConfig {
        units_mode: UnitsMode::SI,
        adiabatic_index: 1.4,
    };

    fn init_noh<P: Physics<E, S>, const E: usize, const S: usize>(u: &mut State<P, E, S>) {
        let breakpoint_index = (S as f64 * 0.5) as usize;
        let mut prim = Array2::zeros((E, S));
        prim.fill(0.0);
        for i in 0..breakpoint_index {
            prim[[0, i]] = 1.0;
            prim[[1, i]] = 1.0;
        }
        for i in breakpoint_index..S {
            prim[[0, i]] = 1.0;
            prim[[1, i]] = -1.0;
        }
        if u.is_adiabatic() {
            prim.row_mut(E - 1).fill(1.0E-5)
        } else {
            let c_sound = Array1::ones(S);
            u.cent.c_sound.assign(&c_sound.view());
        }
        u.cent.prim.assign(&prim.view());
    }

    set_Physics_and_E!(Euler1DIsot);

    #[test]
    fn hll_test() {
        // Also test this on log meshes
        let mesh: Mesh<S> = Mesh::new(&MESHCONFIG).unwrap();
        let mut u = State::<P, E, S>::new(&PHYSICSCONFIG);
        init_noh(&mut u);
        u.update_cons();
        u.update_derived_variables();
        let mut hll: Hll<E, S> = Hll::new(&NumFluxConfig::Hll, &mesh).unwrap();

        let dflux_dxi_prim_expect = Array2::from_shape_vec(
            (EQ, S),
            vec![
                0.0, 0.0, 0.0, -4.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, -8.0, 0.0, 0.0, 0.0,
            ],
        )
        .unwrap();

        let mut dflux_dxi = Array2::zeros((EQ, S));
        hll.calc_dflux_dxi(&mut dflux_dxi, &mut u, &mesh).unwrap();
        assert_relative_eq!(dflux_dxi, dflux_dxi_prim_expect, max_relative = 1.0e-12);
    }
}
