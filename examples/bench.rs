// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::Result;
use corries::{initfuncs::init_noh, prelude::*};

const S: usize = 500;

fn main() -> Result<()> {
    set_Physics_and_E!(Euler1DAdiabatic);
    type N = Kt<E, S>;
    type T = RungeKuttaFehlberg<P, E, S>;

    CorriesConfig::default_riemann_test::<N, E, S>(0.5, "results/examples/bench_noh", "bench_noh")
        .init_corries::<P, N, T, E, S>(init_noh)?
        .run_corries()
}
