// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [Kt] struct.

use color_eyre::{eyre::{bail, ensure, Context}, Result};
use ndarray::{par_azip, s, Array1, Array2};
use crate::{
    directions::Direction,
    errorhandling::Validation,
    mesh::Mesh,
    state::Physics,
    LimiterMode,
    NumFluxConfig,
    State
};
use super::{calc_dflux_xi_generic, NumFlux};

macro_rules! min {
    ($a:expr, $b:expr) => {
        if $a < $b {
            $a
        } else {
            $b
        }
    };
}
macro_rules! max {
    ($a:expr, $b:expr) => {
        if $a > $b {
            $a
        } else {
            $b
        }
    };
}

/// Handles calculating numerical flux using the Kurganov-Tadmor scheme
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Kt<const E: usize, const S: usize> {
    /// Left-side characteristics
    a_plus: Array1<f64>,

    /// Right-side characteristics
    a_minus: Array1<f64>,

    /// Helper array containting `1.0 / (sr - sl)`
    b: Array1<f64>,

    /// Helper array containting `sr * sl`
    c: Array1<f64>,

    /// Numerical flux
    flux_num: Array2<f64>,

    /// The type of limiter to use during reconstruction
    limiter_mode: LimiterMode,

    /// inverse of the xi coordinate differential
    inv_dxi: f64,

    /// monocent parameter, only used when limiter_mode == LimiterMode == Monocent; taken from that
    /// enums payload
    theta: f64,

    /// distances between cell centres and the west facing cell borders
    dist_west: Array1<f64>,

    /// distances between cell centres and the east facing cell borders
    dist_east: Array1<f64>,
}

unsafe impl<const E: usize, const S: usize> Send for Kt<E, S> {}
unsafe impl<const E: usize, const S: usize> Sync for Kt<E, S> {}

impl<const E: usize, const S: usize> NumFlux<E, S> for Kt<E, S> {
    fn new(numflux_config: &NumFluxConfig, mesh: &Mesh<S>) -> Result<Self> {
        match numflux_config {
            NumFluxConfig::Kt { limiter_mode } => Ok(Self {
                a_plus: Array1::zeros(S),
                a_minus: Array1::zeros(S),
                b: Array1::zeros(S),
                c: Array1::zeros(S),
                flux_num: Array2::zeros((E, S)),
                limiter_mode: *limiter_mode,
                inv_dxi: 1.0 / mesh.dxi,
                theta: match limiter_mode {
                    LimiterMode::Monocent(x) => *x,
                    _ => 0.0,
                },
                dist_west: Array1::from_shape_fn(S, |i| mesh.xi_west[i] - mesh.xi_cent[i]),
                dist_east: Array1::from_shape_fn(S, |i| mesh.xi_east[i] - mesh.xi_cent[i]),
            }),
            _ => bail!("Tried constructing Kt, but numflux_config does not contain NumFluxConfig::Kt!"),
        }
    }

    fn calc_dflux_dxi<P: Physics<E, S>>(
        &mut self,
        dflux_dxi: &mut Array2<f64>,
        u: &mut State<P, E, S>,
        mesh: &Mesh<S>,
    ) -> Result<()> {
        // NOTE: Assumes that u.cons is already up to date
        self.reconstruct(u);

        u.update_prim_d::<{ Direction::West as u8 }>();
        u.update_prim_d::<{ Direction::East as u8 }>();
        u.update_c_sound_d::<{ Direction::West as u8 }>();
        u.update_c_sound_d::<{ Direction::East as u8 }>();
        u.update_eigen_vals_min_d::<{ Direction::West as u8 }>();
        u.update_eigen_vals_min_d::<{ Direction::East as u8 }>();
        u.update_eigen_vals_max_d::<{ Direction::West as u8 }>();
        u.update_eigen_vals_max_d::<{ Direction::East as u8 }>();
        u.update_flux_d::<{ Direction::West as u8 }>();
        u.update_flux_d::<{ Direction::East as u8 }>();

        let eigen_min_west = &u.west.eigen_min();
        let eigen_min_east = &u.east.eigen_min();
        let eigen_max_west = &u.west.eigen_max();
        let eigen_max_east = &u.east.eigen_max();

        let uflux_west = &u.west.flux;
        let uflux_east = &u.east.flux;
        let ucons_west = &u.west.cons;
        let ucons_east = &u.east.cons;

        let s = s![(mesh.ixi_in - 1)..=mesh.ixi_out];
        let sp1 = s![mesh.ixi_in..=(mesh.ixi_out + 1)];

        par_azip!((
                a_plus in &mut self.a_plus.slice_mut(s),
                &ev_max_west in &eigen_max_west.slice(sp1),
                &ev_max_east in &eigen_max_east.slice(s))
            *a_plus = max!(0.0, max!(ev_max_west, ev_max_east)));
        par_azip!((
                a_minus in &mut self.a_minus.slice_mut(s),
                &ev_min_west in &eigen_min_west.slice(sp1),
                &ev_min_east in &eigen_min_east.slice(s))
            *a_minus = min!(0.0, min!(ev_min_west, ev_min_east)));

        par_azip!((
                b in &mut self.b.slice_mut(s),
                &y in &mesh.d_area_xi_deta_dphi_east.slice(s),
                &ap in &self.a_plus.slice(s),
                &am in &self.a_minus.slice(s))
            *b = y / (ap - am));
        par_azip!((
                c in &mut self.c.slice_mut(s),
                &ap in &self.a_plus.slice(s),
                &am in &self.a_minus.slice(s))
            *c = ap * am);

        for j in 0..E {
            self.flux_num.row_mut(j).slice_mut(s).assign(
                &(&self.b.slice(s)
                    * (&self.a_plus.slice(s) * &uflux_east.row(j).slice(s)
                        - &self.a_minus.slice(s) * &uflux_west.row(j).slice(sp1)
                        + &self.c.slice(s) * (&ucons_west.row(j).slice(sp1) - &ucons_east.row(j).slice(s)))),
            );
        }

        calc_dflux_xi_generic::<E, S>(dflux_dxi, &self.flux_num, mesh);
        if cfg!(feature = "validation") {
            self.validate().context("Calling Kt::validate in Kt::calc_dflux_dxi")?;
        }
        Ok(())
    }
}

fn signum(a: f64) -> i32 {
    if a > 0.0 {
        1
    } else {
        -1
    }
}

impl<const E: usize, const S: usize> Kt<E, S> {
    fn reconstruct<P: Physics<E, S>>(&self, u: &mut State<P, E, S>) {
        let slope_fn = match self.limiter_mode {
            LimiterMode::NoLimiter => |a: f64, b: f64, _: f64, _: f64| 0.5 * (a + b),
            LimiterMode::MinMod => |a: f64, b: f64, _: f64, _: f64| {
                if signum(a) * signum(b) > 0 {
                    signum(a) as f64 * min!(a.abs(), b.abs())
                } else {
                    0.0
                }
            },
            LimiterMode::Superbee => |a: f64, b: f64, _: f64, _: f64| {
                if a * b > 0.0 {
                    signum(a) as f64 * min!(min!(a.abs(), b.abs()), 0.5 * max!(a.abs(), b.abs()))
                } else {
                    0.0
                }
            },
            LimiterMode::Monocent(_) => |a: f64, b: f64, c: f64, p: f64| {
                if signum(a) * signum(b) > 0 && signum(b) * signum(c) > 0 {
                    signum(a) as f64 * min!((p * a).abs(), min!((p * b).abs(), c.abs()))
                } else {
                    0.0
                }
            },
            LimiterMode::VanLeer => |a: f64, b: f64, _: f64, _: f64| {
                let abs_a = a.abs();
                let abs_b = b.abs();
                (a * abs_b + b * abs_a) / (abs_a + abs_b + f64::MIN)
            },
        };
        for j in 0..E {
            for i in 1..S - 1 {
                let slope = self.inv_dxi
                    * slope_fn(
                        u.cent.cons[[j, i]] - u.cent.cons[[j, i - 1]],
                        u.cent.cons[[j, i + 1]] - u.cent.cons[[j, i]],
                        0.5 * (u.cent.cons[[j, i + 1]] - u.cent.cons[[j, i - 1]]),
                        self.theta,
                    );
                u.west.cons[[j, i]] = u.cent.cons[[j, i]] + slope * self.dist_west[i];
                u.east.cons[[j, i]] = u.cent.cons[[j, i]] + slope * self.dist_east[i];
            }
        }
    }
}

impl<const E: usize, const S: usize> Validation for Kt<E, S> {
    fn validate(&self) -> Result<()> {
        ensure!(
            self.a_plus.fold(true, |acc, x| acc && x.is_finite()),
            "Kt::a_plus must be finite! Got: {}",
            self.a_plus
        );
        ensure!(
            self.a_minus.fold(true, |acc, x| acc && x.is_finite()),
            "Kt::a_minus must be finite! Got: {}",
            self.a_minus
        );
        ensure!(
            self.b.fold(true, |acc, x| acc && x.is_finite()),
            "Kt::b = dA / (deta * dPhi * (a_plus - a_minus)) must be finite! Got:\nb = {}\na_plus = {}\na_minus = {}",
            self.b,
            self.a_plus,
            self.a_minus
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
    fn kt_test() {
        // Also test this on log meshes
        let mesh: Mesh<S> = Mesh::new(&MESHCONFIG).unwrap();
        let mut u = State::<P, E, S>::new(&PHYSICSCONFIG);
        init_noh(&mut u);
        u.update_cons();
        u.update_derived_variables();
        u.init_west_east();
        let mut kt = Kt::new(
            &NumFluxConfig::Kt {
                limiter_mode: LimiterMode::Monocent(1.2),
            },
            &mesh,
        )
        .unwrap();

        let dflux_dxi_prim_expect = Array2::from_shape_vec(
            (EQ, S),
            vec![
                0.0, 0.0, 0.0, -4.0, -4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, -8.0, 0.0, 0.0, 0.0,
            ],
        )
        .unwrap();

        let mut dflux_dxi = Array2::zeros((EQ, S));
        kt.calc_dflux_dxi(&mut dflux_dxi, &mut u, &mesh).unwrap();
        assert_relative_eq!(dflux_dxi, dflux_dxi_prim_expect, max_relative = 1.0e-12);
    }
}
