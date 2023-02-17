// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::Result;
use corries::prelude::*;
use ndarray::{Array1, Array2};

const S: usize = 500;

fn init<P: Physics<E, S>, N: NumFlux<E, S>, const E: usize, const S: usize>(
    u: &mut State<P, E, S>,
    rhs: &mut Rhs<N, E, S>,
    mesh: &Mesh<S>,
) {
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
    u.update_vars_from_prim(&mut rhs.boundary_west, &mut rhs.boundary_east, &mesh);
    u.init_west_east();
    return;
}

fn main() -> Result<()> {
    set_Physics_and_E!(Euler1DAdiabatic);
    type N = Kt<E, S>;
    type T = RungeKuttaFehlberg<P, E, S>;

    let config = CorriesConfig::default_riemann_test::<N, E, S>(0.5, "results/examples/bench_noh", "bench_noh");
    let (mut u, mut rhs, mut time, mesh, mut writer) = init_corries::<P, N, T, E, S>(&config).unwrap();

    init(&mut u, &mut rhs, &mesh);
    run_corries::<P, N, T, E, S>(&mut u, &mut rhs, &mut time, &mesh, &mut writer)
}
