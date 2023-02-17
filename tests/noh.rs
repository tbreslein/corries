// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::{eyre::Context, Result};
use corries::prelude::*;
const S: usize = 100;

fn init<P, N, T, const E: usize, const S: usize>(
    u: &mut State<P, E, S>,
    _: &mut Solver<P, N, T, E, S>,
    _: &Mesh<S>,
) -> Result<()>
where
    P: Physics<E, S>,
    N: NumFlux<E, S>,
    T: TimeSolver<P, E, S>,
{
    let breakpoint_index = (S as f64 * 0.5) as usize;
    u.cent.prim.fill(1.0);
    for i in breakpoint_index..S {
        u.cent.prim[[1, i]] = -1.0;
    }
    if u.is_adiabatic() {
        u.cent.prim.row_mut(E - 1).fill(1.0E-5);
    } else {
        u.cent.c_sound.fill(1.0);
    }
    Ok(())
}

#[test]
fn noh_euler1d_adiabatic() -> Result<()> {
    set_Physics_and_E!(Euler1DAdiabatic);
    type N = Hll<E, S>;
    type T = RungeKuttaFehlberg<P, E, S>;
    CorriesConfig::default_riemann_test::<N, E, S>(
        0.5,
        "results/integrationtests/noh_euler1d_adiabatic",
        "noh_euler1d_adiabatic",
    )
    .init_corries::<P, N, T, E, S>(init)
    .context("While calling CorriesConfig::init_corries")?
    .run_corries()
}

#[test]
fn noh_euler1d_isot() -> Result<()> {
    set_Physics_and_E!(Euler1DIsot);
    type N = Hll<E, S>;
    type T = RungeKuttaFehlberg<P, E, S>;
    CorriesConfig::default_riemann_test::<N, E, S>(
        0.5,
        "results/integrationtests/noh_euler1d_isothermal",
        "noh_euler1d_isothermal",
    )
    .init_corries::<P, N, T, E, S>(init)
    .context("While calling CorriesConfig::init_corries")?
    .run_corries()
}
