// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::{eyre::Context, Result};
use config::{Config, CorriesConfig};
use mesh::Mesh;

pub mod config;
#[macro_use]
mod errorhandling;
mod mesh;

/// Runs a Corries simulation.
///
/// # Arguments
///
/// * `config` - The `CorriesConfig` the simulation is based on
pub fn run_sim(config: CorriesConfig) -> Result<()> {
    config.validate().context("Validating config")?;
    let _mesh = Mesh::new(&config.meshconf)?;
    return Ok(());
}
