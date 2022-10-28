// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::{eyre::Context, Result};
use config::CorriesConfig;
use mesh::Mesh;
use writer::Writer;

use crate::{errorhandling::Validation, rhs::Rhs};

pub mod config;
#[macro_use]
mod errorhandling;
mod mesh;
mod rhs;
mod writer;

/// Runs a Corries simulation.
///
/// # Arguments
///
/// * `config` - The `CorriesConfig` the simulation is based on
pub fn run_sim(config: CorriesConfig) -> Result<()> {
    config.validate().context("Validating config")?;
    let mesh = Mesh::new(&config.meshconf).context("Constructing Mesh")?;
    let _rhs = Rhs::new();

    // TEMP:
    let output_count_max = 2;
    let mut writer = Writer::new(&config, &mesh, output_count_max)?;

    // first output
    writer.update_data_matrices(&mesh)?;
    writer.write_metadata(&config)?;
    println!("{}", mesh.xi_cent);
    writer.write_output()?;
    return Ok(());
}
