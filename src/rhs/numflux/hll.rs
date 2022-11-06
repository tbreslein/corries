// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [Hll] struct.

use color_eyre::eyre::{ensure, Context};
use color_eyre::Result;
use ndarray::{s, Array1, Array2};

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
        u.calc_physical_flux_euler1d_adiabatic(&mut self.flux_phys);
        let slice = s![mesh.ixi_in - 1..=mesh.ixi_out];
        let slice_p1 = s![mesh.ixi_in..=mesh.ixi_out + 1];

        for i in mesh.ixi_in - 1..=mesh.ixi_out {
            self.sl[i] = 0.0f64.min(u.eigen_vals[[0, i]].min(u.eigen_vals[[0, i + 1]]));
            self.sr[i] = 0.0f64.max(u.eigen_vals[[EQ - 1, i]].max(u.eigen_vals[[EQ - 1, i + 1]]));
        }

        // TODO: Benchmark this version against azip, Zip and raw for loops
        self.inv_sr_minus_sl
            .slice_mut(slice)
            .assign(&(1.0 / (&self.sr.slice(slice) - &self.sl.slice(slice))));
        self.sr_times_sl
            .slice_mut(slice)
            .assign(&(&self.sr.slice(slice) * &self.sl.slice(slice)));

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

            flux_num_j.slice_mut(slice).assign(
                &(&a * (&sr * &flux_phys_j.slice(slice) - &sl * &flux_phys_j.slice(slice_p1)
                    + &b * (&uc_j.slice(slice_p1) - &uc_j.slice(slice)))),
            );
        }
        calc_dflux_xi_generic(dflux_dxi, &self.flux_num, mesh);

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
    const S: usize = 6;
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
        c_sound_0: 1.0,
    };

fn init_noh<const S: usize, const EQ: usize>(u: &mut Physics<S, EQ>) {
    let breakpoint_index = ((S - 1) as f64 * 0.5) as usize;
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
    return;
}

    #[test]
    fn hll_test() {
        // Also test this on log meshes
        let mesh: Mesh<S> = Mesh::new(&MESHCONFIG).unwrap();
        let mut u: Physics<S, EQ> = Physics::new(&PHYSICSCONFIG);
        init_noh(&mut u);
        let mut hll: Hll<S, EQ> = Hll::new();

        let u_prim_expect = Array2::from_shape_vec(
            (EQ, S),
            vec![0.0, 0.0, -4.0, -4.0, 0.0, 0.0, 0.0, 0.0, 8.0, -8.0, 0.0, 0.0],
        )
        .unwrap();

        let mut dflux_dxi = Array2::zeros((EQ, S));
        hll.calc_dflux_dxi(&mut dflux_dxi, &mut u, &mesh).unwrap();
        assert_relative_eq!(u.prim, u_prim_expect, max_relative = 1.0e-12);
    }
}
