use color_eyre::eyre::ensure;
use color_eyre::Result;

use super::Config;

#[derive(Debug)]
pub enum MeshMode {
    Cartesian,
}

#[derive(Debug)]
pub struct MeshConfig {
    pub mode: MeshMode,
    pub n_comp: u32,
    pub n_gc: u32,
    pub xi_in: f64,
    pub xi_out: f64,
    pub ratio: f64,
}

impl Config for MeshConfig {
    fn validate(&self) -> Result<()> {
        ensure!(
            self.xi_in < self.xi_out,
            "This must hold: xi_in < xi_out! Got xi_in = {} ; xi_out = {}",
            self.xi_in,
            self.xi_out
        );
        ensure!(
            self.ratio > 0.0 && self.ratio <= 1.0,
            "This must hold: 0.0 < ratio <= 1.0! Got ratio = {}",
            self.ratio
        );
        return Ok(());
    }
}
