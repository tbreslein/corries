// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::Result;
use corries::prelude::*;

const S: usize = 500;

fn init<P: Physics<E, S>, N: NumFlux<E, S>, T: TimeSolver<P, E, S>, const E: usize, const S: usize>(
    u: &mut State<P, E, S>,
    _: &mut Solver<P, N, T, E, S>,
    _: &Mesh<S>,
) -> Result<()> {
    let breakpoint_index = (S as f64 * 0.5) as usize;
    u.cent.prim.fill(1.0);
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

fn main() -> Result<()> {
    set_Physics_and_E!(Euler1DAdiabatic);
    type N = Kt<E, S>;
    type T = RungeKuttaFehlberg<P, E, S>;

    CorriesConfig::default_riemann_test::<N, E, S>(0.5, "results/examples/bench_noh", "bench_noh")
        .init_corries::<P, N, T, E, S>(init)?
        .run_corries()
}
