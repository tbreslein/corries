// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::{eyre::Context, Result};
use corries::prelude::*;
use ndarray::{Array1, Array2};
const S: usize = 100;

fn init_noh<P: Physics<E, S>, const E: usize, const S: usize>(u: &mut State<P, E, S>) {
    let breakpoint_index = (S as f64 * 0.5) as usize;
    let mut prim = Array2::ones((E, S));
    for i in breakpoint_index..S {
        prim[[1, i]] = -1.0;
    }
    if u.is_adiabatic() {
        prim.row_mut(E - 1).fill(1.0E-5);
    } else {
        u.cent.c_sound.assign(&Array1::ones(S).view());
    }
    u.cent.prim.assign(&prim.view());
    return;
}

#[test]
fn noh_euler1d_adiabatic() -> Result<()> {
    set_Physics_and_E!(Euler1DAdiabatic);
    type N = Hll<E, S>;
    type T = RungeKuttaFehlberg<P, E, S>;
    let t_end = 0.5;

    let config = CorriesConfig::default_riemann_test::<N, E, S>(
        t_end,
        "results/integrationtests/noh_euler1d_adiabatic",
        "noh_euler1d_adiabatic",
    );
    let (mut u, mut rhs, mut time, mesh, mut writer) = init_corries::<P, N, T, E, S>(&config).unwrap();

    init_noh::<P, E, S>(&mut u);
    u.update_vars_from_prim(&mut rhs.boundary_west, &mut rhs.boundary_east, &mesh);

    run_corries::<P, N, T, E, S>(&mut u, &mut rhs, &mut time, &mesh, &mut writer)
        .context("Calling run_loop in noh test")?;
    return Ok(());
}

#[test]
fn noh_euler1d_isot() -> Result<()> {
    set_Physics_and_E!(Euler1DIsot);
    type N = Hll<E, S>;
    type T = RungeKuttaFehlberg<P, E, S>;
    let t_end = 0.5;

    let config = CorriesConfig::default_riemann_test::<N, E, S>(
        t_end,
        "results/integrationtests/noh_euler1d_isothermal",
        "noh_euler1d_isothermal",
    );
    let (mut u, mut rhs, mut time, mesh, mut writer) = init_corries::<P, N, T, E, S>(&config).unwrap();

    init_noh::<P, E, S>(&mut u);
    u.update_vars_from_prim(&mut rhs.boundary_west, &mut rhs.boundary_east, &mesh);

    run_corries::<P, N, T, E, S>(&mut u, &mut rhs, &mut time, &mesh, &mut writer)
        .context("Calling run_loop in noh test")?;
    return Ok(());
}
