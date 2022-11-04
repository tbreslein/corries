// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

use color_eyre::{eyre::Context, Result};

use crate::{
    config::{CorriesConfig, InitialConditions},
    physics::Physics,
};

use self::noh::init_noh;

mod noh;

pub fn apply_initial_conditions<const S: usize, const EQ: usize>(
    config: &CorriesConfig,
    u: &mut Physics<S, EQ>,
) -> Result<()> {
    return match config.initial_conditions {
        InitialConditions::Noh => init_noh(u).context("Applying Noh conditions"),
    };
}
