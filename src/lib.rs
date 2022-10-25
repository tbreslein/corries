use color_eyre::{Result, eyre::Context};
use config::{Config, CorriesConfig};

pub mod config;

pub fn run_sim(config: CorriesConfig) -> Result<()> {
    config.validate().context("Validating config")?;
    return Ok(());
}
