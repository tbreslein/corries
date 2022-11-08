// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [Hll] struct.

use color_eyre::eyre::{ensure, Context};
use color_eyre::Result;
use ndarray::{par_azip, s, Array1, Array2};

use crate::errorhandling::Validation;
use crate::{mesh::Mesh, physics::Physics};

use super::calc_dflux_xi_generic;
use super::NumFlux;

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
        return Self {
            flux_phys: Array2::zeros((EQ, S)),
            sl: Array1::zeros(S),
            sr: Array1::zeros(S),
            inv_sr_minus_sl: Array1::zeros(S),
            sr_times_sl: Array1::zeros(S),
            flux_num: Array2::zeros((EQ, S)),
        };
    }
}

impl<const S: usize, const EQ: usize> NumFlux<S, EQ> for Hll<S, EQ> {
    fn calc_dflux_dxi(&mut self, dflux_dxi: &mut Array2<f64>, u: &mut Physics<S, EQ>, mesh: &Mesh<S>) -> Result<()> {
        // NOTE: Assumes that u.eigen_vals are already up to date
        // NOTE: This function was benchmarked using many different combinations of raw loops,
        // azips, par_azips, and Zip::par_for_each. This ended up being the most performant
        // version. The lession we can take from this is that, that the BLAS parallelisation used
        // in the .assign method is already good enough, and that parallelisation over the
        // equations creates to much overhead. This function in particular ran about 7.5% slower in
        // my notebook when replacing the for j loop with a Zip::par_for_each.
        // My guess is that the for loop probably got fully unrolled, since the compiler knows the
        // length of the loop at compile time, and that lead to a far greater performance boost
        // compared to what the Zip::par_for_each had to offer.
        // The assignments for sl, sr, inv_sr_minus_sl, and sr_times_sl on the other hand performed
        // just well in both version: simple .assign calls and par_azip calls. I chose to stick
        // with the par_azip version because I think it's easier to read, and might(!) scale
        // better.
        u.calc_physical_flux(&mut self.flux_phys);
        let s = s![(mesh.ixi_in - 1)..=mesh.ixi_out];
        let sp1 = s![mesh.ixi_in..=(mesh.ixi_out + 1)];

        par_azip!((sl in &mut self.sl.slice_mut(s), &ev1 in &u.eigen_vals.row(0).slice(s), &ev2 in &u.eigen_vals.row(0).slice(sp1))
            *sl = 0.0f64.min(ev1.min(ev2))
        );
        par_azip!((sr in &mut self.sr.slice_mut(s), &ev1 in &u.eigen_vals.row(EQ-1).slice(s), &ev2 in &u.eigen_vals.row(EQ-1).slice(sp1))
            *sr = 0.0f64.max(ev1.max(ev2))
        );
        par_azip!((a in &mut self.inv_sr_minus_sl.slice_mut(s), &sr in &self.sr.slice(s), &sl in &self.sl.slice(s))
            *a = 1.0 / (sr - sl)
        );
        par_azip!((b in &mut self.sr_times_sl.slice_mut(s), &sr in &self.sr.slice(s), &sl in &self.sl.slice(s))
            *b = sr * sl
        );

        for j in 0..EQ {
            let mut flux_num_j = self.flux_num.row_mut(j);
            let flux_phys_j = self.flux_phys.row(j);
            let uc_j = u.cons.row(j);
            flux_num_j.slice_mut(s).assign(
                &(&self.inv_sr_minus_sl.slice(s)
                    * (&self.sr.slice(s) * &flux_phys_j.slice(s) - &self.sl.slice(s) * &flux_phys_j.slice(sp1)
                        + &self.sr_times_sl.slice(s) * (&uc_j.slice(sp1) - &uc_j.slice(s)))),
            );
        }

        calc_dflux_xi_generic::<S, EQ>(dflux_dxi, &self.flux_num, mesh);
        self.validate()
            .context("Calling Hll::validate in Hll::calc_dflux_dxi")?;
        return Ok(());
    }
}

impl<const S: usize, const EQ: usize> Validation for Hll<S, EQ> {
    fn validate(&self) -> Result<()> {
        ensure!(
            self.flux_phys.fold(true, |acc, x| acc && x.is_finite()),
            "Hll::flux_phys must be finite! Got: {}",
            self.flux_phys
        );
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
        return Ok(());
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        config::{
            meshconfig::{MeshConfig, MeshMode},
            physicsconfig::{PhysicsConfig, PhysicsMode},
        },
        units::UnitsMode,
    };

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
        mode: PhysicsMode::Euler1DIsot,
        units_mode: UnitsMode::SI,
        adiabatic_index: 1.4,
    };

    fn init_noh<const S: usize, const EQ: usize>(u: &mut Physics<S, EQ>) {
        let breakpoint_index = (S as f64 * 0.5) as usize;
        u.prim.fill(0.0);
        u.cons.fill(0.0);
        for i in 0..breakpoint_index {
            u.prim[[u.jdensity, i]] = 1.0;
            u.prim[[u.jxivelocity, i]] = 1.0;
        }
        for i in breakpoint_index..S {
            u.prim[[u.jdensity, i]] = 1.0;
            u.prim[[u.jxivelocity, i]] = -1.0;
        }
        u.c_sound.fill(1.0);
        return;
    }

    #[test]
    fn hll_test() {
        // Also test this on log meshes
        let mesh: Mesh<S> = Mesh::new(&MESHCONFIG).unwrap();
        let mut u: Physics<S, EQ> = Physics::new(&PHYSICSCONFIG);
        init_noh(&mut u);
        u.update_cons();
        u.update_derived_variables();
        let mut hll: Hll<S, EQ> = Hll::new();

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
