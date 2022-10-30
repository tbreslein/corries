// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::{eyre::Context, Result};
use config::CorriesConfig;
use mesh::Mesh;
use physics::init_physics;
use writer::Writer;

use crate::{config::outputconfig::StreamMode, errorhandling::Validation, rhs::Rhs};

pub mod config;
#[macro_use]
mod errorhandling;
mod mesh;
mod physics;
mod rhs;
pub mod units;
mod writer;

/// Runs a Corries simulation.
///
/// # Arguments
///
/// * `config` - The `CorriesConfig` the simulation is based on
pub fn run_sim(config: CorriesConfig) -> Result<()> {
    config.validate().context("Validating config")?;
    let mesh = Mesh::new(&config.meshconf).context("Constructing Mesh")?;
    let _u = init_physics(&config.physicsconf, &mesh).context("Constructing Physics")?;
    let _rhs = Rhs::new();

    // TEMP:
    let output_count_max = 2;
    let mut writer = Writer::new(&config, &mesh, output_count_max)?;

    // first output
    if writer.outputs.iter().any(|o| o.stream_mode == StreamMode::Stdout) {
        print_banner();
    }
    writer.update_data_matrices(&mesh)?;
    writer.write_metadata(&config)?;
    writer.write_output()?;
    return Ok(());
}

const VERSION: &str = env!("CARGO_PKG_VERSION");
fn print_banner() {
    println!("# ********************************************");
    println!("# *** Corries - corrosive Riemann solver ");
    println!("# *** ");
    println!("# *** Version: {}", VERSION);
    println!("# *** Copyright (c) 2022");
    println!("# *** Author: tbreslein <github.com/tbreslein>");
    println!("# *** License: MIT");
    println!("# ********************************************");
}
