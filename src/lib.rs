// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::{eyre::Context, Result};
use config::CorriesConfig;
use mesh::Mesh;
use writer::Writer;

use crate::errorhandling::Validation;

pub mod config;
#[macro_use]
mod errorhandling;
mod mesh;
mod writer;

/// Runs a Corries simulation.
///
/// # Arguments
///
/// * `config` - The `CorriesConfig` the simulation is based on
pub fn run_sim(config: CorriesConfig) -> Result<()> {
    config.validate().context("Validating config")?;
    let mesh = Mesh::new(&config.meshconf).context("Constructing Mesh")?;
    let _writer = Writer::new(&config, &mesh);
    return Ok(());
}
