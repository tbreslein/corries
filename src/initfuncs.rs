// Copyright (c) 2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports commonly used initial condition function that are most importantly used in the
//! integration tests. These functions all can be passed to
//! [CorriesConfig::init_corries()](crate::config::CorriesConfig::init_corries()).

use crate::{Euler1DAdiabatic, Euler1DIsot, Mesh, NumFlux, Physics, Solver, State, TimeSolver};
use color_eyre::{eyre::bail, Result};
use std::any::TypeId;

/// The form an initial condition function needs to have.
///
/// This includes passing a [State], a [Solver] as well as a [Mesh] to the function, and the
/// function should return a [color_eyre::Result<()>].
pub type InitFn<P, N, T, const E: usize, const S: usize> =
    fn(&mut State<P, E, S>, &mut Solver<P, N, T, E, S>, &Mesh<S>) -> Result<()>;

/// Sets up the initial conditions for the Noh test.
///
/// This sets up:
///
/// * mass density of 1.0 everywhere
/// * xi velocity of 1.0 on the left hand side of the shocktube
/// * xi velocity of -1.0 on the right hand side of the shocktube
/// * pressure of 1.0E-5 in case of adiabatic physics, or speed of sound of 1.0 otherwise
///
/// ```
/// use corries::prelude::*;
///
/// const S: usize = 100;
/// set_Physics_and_E!(Euler1DAdiabatic);
/// type N = Kt<E, S>;
/// type T = RungeKuttaFehlberg<P, E, S>;
///
/// let (u, _, _, _) = CorriesConfig::default_riemann_test::<N, E, S>(
///     0.5,
///     "results/integrationtests/noh_euler1d_adiabatic",
///     "noh_euler1d_adiabatic",
/// )
/// .init_corries::<P, N, T, E, S>(corries::initfuncs::init_noh).unwrap();
///
///
/// use approx::assert_relative_eq;
/// use ndarray::Array2;
///
/// let breakpoint_index = (S as f64 * 0.5) as usize;
///
/// let mut u0_prim = Array2::zeros((E, S)) ;
/// u0_prim.row_mut(P::JRHO).fill(1.0);
/// u0_prim.row_mut(P::JPRESSURE).fill(1.0E-5);
/// for i in 0..breakpoint_index {
///     u0_prim[[P::JXI, i]] = 1.0;
/// }
/// for i in breakpoint_index..S {
///     u0_prim[[P::JXI, i]] = -1.0;
/// }
///
/// assert_relative_eq!(u.cent.prim, u0_prim);
/// ```
pub fn init_noh<P, N, T, const E: usize, const S: usize>(
    u: &mut State<P, E, S>,
    _: &mut Solver<P, N, T, E, S>,
    _: &Mesh<S>,
) -> Result<()>
where
    P: Physics<E, S> + 'static,
    N: NumFlux<E, S>,
    T: TimeSolver<P, E, S>,
{
    // TODO:
    // I can probably turn this into a static assert by attaching a constant identifier to
    // my P type that I can then static assert here
    if TypeId::of::<P>() != TypeId::of::<Euler1DIsot<S>>() && TypeId::of::<P>() != TypeId::of::<Euler1DAdiabatic<S>>() {
        bail!("init_noh cannot run when the Physics type is set to: {}!", P::name())
    };
    let breakpoint_index = (S as f64 * 0.5) as usize;
    u.cent.prim.row_mut(P::JRHO).fill(1.0);

    for i in 0..breakpoint_index {
        u.cent.prim[[P::JXI, i]] = 1.0;
    }
    for i in breakpoint_index..S {
        u.cent.prim[[P::JXI, i]] = -1.0;
    }

    if u.is_adiabatic() {
        u.cent.prim.row_mut(P::JPRESSURE).fill(1.0E-5);
    } else {
        u.cent.c_sound.fill(1.0);
    }

    Ok(())
}

/// Sets up the initial conditions for the Sod test.
///
/// This sets up:
///
/// * mass density and pressure of 1.0 on the left hand side of the shocktube
/// * mass density of 0.125 and pressure of 0.1 on the right hand side of the shocktube
/// * all velocities set to 0.0
///
/// ```
/// use corries::prelude::*;
///
/// const S: usize = 100;
/// set_Physics_and_E!(Euler1DAdiabatic);
/// type N = Kt<E, S>;
/// type T = RungeKuttaFehlberg<P, E, S>;
///
/// let (u, _, _, _) = CorriesConfig::default_riemann_test::<N, E, S>(
///     0.5,
///     "results/integrationtests/noh_euler1d_adiabatic",
///     "noh_euler1d_adiabatic",
/// )
/// .init_corries::<P, N, T, E, S>(corries::initfuncs::init_sod).unwrap();
///
///
/// use approx::assert_relative_eq;
/// use ndarray::Array2;
///
/// let breakpoint_index = (S as f64 * 0.5) as usize;
///
/// let mut u0_prim = Array2::zeros((E, S)) ;
/// for i in 0..breakpoint_index {
///     u0_prim[[P::JRHO, i]] = 1.0;
///     u0_prim[[P::JPRESSURE, i]] = 1.0;
/// }
/// for i in breakpoint_index..S {
///     u0_prim[[P::JRHO, i]] = 0.125;
///     u0_prim[[P::JPRESSURE, i]] = 0.1;
/// }
///
/// assert_relative_eq!(u.cent.prim, u0_prim);
/// ```
pub fn init_sod<P, N, T, const E: usize, const S: usize>(
    u: &mut State<P, E, S>,
    _: &mut Solver<P, N, T, E, S>,
    _: &Mesh<S>,
) -> Result<()>
where
    P: Physics<E, S> + 'static,
    N: NumFlux<E, S>,
    T: TimeSolver<P, E, S>,
{
    // TODO:
    // I can probably turn this into a static assert by attaching a constant identifier to
    // my P type that I can then static assert here
    if TypeId::of::<P>() != TypeId::of::<Euler1DAdiabatic<S>>() {
        bail!("init_sod cannot run when the Physics type is set to: {}!", P::name())
    };

    let breakpoint_index = (S as f64 * 0.5) as usize;
    u.cent.prim.fill(0.0);
    for i in 0..breakpoint_index {
        u.cent.prim[[P::JRHO, i]] = 1.0;
        u.cent.prim[[P::JPRESSURE, i]] = 1.0;
    }
    for i in breakpoint_index..S {
        u.cent.prim[[P::JRHO, i]] = 0.125;
        u.cent.prim[[P::JPRESSURE, i]] = 0.1;
    }
    Ok(())
}
