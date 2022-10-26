// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::eyre::{ensure, Context};
use color_eyre::Result;

use self::meshconfig::MeshConfig;
use self::outputconfig::OutputConfig;
use crate::errorhandling::Validation;

pub mod meshconfig;
pub mod outputconfig;

/// Struct that carries the full configuration info for a simulation.
///
/// This struct is used in the beginning of a run to initialise all the runtime-objects that are
/// used throughout the simulation.
#[derive(Debug)]
pub struct CorriesConfig {
    /// Name of the simulation; may not be empty
    pub name: String,

    /// Config for Mesh objects
    pub meshconf: MeshConfig,

    /// Config for Writer objects
    pub writerconf: Vec<OutputConfig>,
}

impl Validation for CorriesConfig {
    fn validate(&self) -> Result<()> {
        ensure!(!self.name.is_empty(), "name must not be empty!");
        self.meshconf
            .validate()
            .context("Validating config.meshconf")?;
        for outputconf in self.writerconf.iter() {
            outputconf.validate().context("Validating config.writerconf")?;
        }
        return Ok(());
    }
}
