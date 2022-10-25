// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::{Result, eyre::Context};
use config::{Config, CorriesConfig};
use mesh::Mesh;

pub mod config;
mod mesh;

/// Runs a Corries simulation.
///
/// # Arguments
///
/// * `config` - The `CorriesConfig` the simulation is based on
pub fn run_sim(config: CorriesConfig) -> Result<()> {
    config.validate().context("Validating config")?;
    let _mesh = Mesh::new(&config.meshconf);
    return Ok(());
}
