use color_eyre::Result;
use color_eyre::eyre::{ensure, Context};

use self::meshconfig::MeshConfig;

pub mod meshconfig;

#[derive(Debug)]
pub struct CorriesConfig {
    pub name: String,
    pub meshconf: MeshConfig,
}

pub trait Config {
    fn validate(&self) -> Result<()>;
}

impl Config for CorriesConfig {
    fn validate(&self) -> Result<()> {
        ensure!(!self.name.is_empty(), "name must not be empty!");
        self.meshconf.validate().context("Validating config.meshconf")?;
        return Ok(());
    }
}
