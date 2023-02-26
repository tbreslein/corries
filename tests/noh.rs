// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::{eyre::Context, Result};
use corries::{initfuncs::init_noh, prelude::*};
const S: usize = 100;

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
    .init_corries::<P, N, T, E, S>(init_noh)
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
    .init_corries::<P, N, T, E, S>(init_noh)
    .context("While calling CorriesConfig::init_corries")?
    .run_corries()
}
