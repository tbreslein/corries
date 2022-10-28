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
            outputconf
                .validate()
                .context("Validating config.writerconf")?;
        }
        return Ok(());
    }
}

impl CorriesConfig {
    pub fn metadata_dump(&self) -> String {
        let mut s = "".to_string();
        s += "### Corries Configuration\n";
        s += &format!("# name: {}\n", self.name);
        s += "# meshconf:\n";
        s += &format!("#     mode: {:?}\n", self.meshconf.mode);
        s += &format!("#     n_comp: {}\n", self.meshconf.n_comp);
        s += &format!("#     n_gc: {}\n", self.meshconf.n_gc);
        s += &format!("#     xi_in: {}\n", self.meshconf.xi_in);
        s += &format!("#     xi_out: {}\n", self.meshconf.xi_out);
        s += &format!("#     ratio_disk: {}\n", self.meshconf.ratio_disk);
        s += "# writerconf:\n";
        for (i, outputconf) in self.writerconf.iter().enumerate() {
            s += &format!("#   outputconf[{}]\n", i);
            s += &format!("#     stream_mode: {:?}\n", outputconf.stream_mode);
            s += &format!("#     formatter_mode: {:?}\n", outputconf.formatter_mode);
            s += &format!(
                "#     string_conversion_mode: {:?}\n",
                outputconf.string_conversion_mode
            );
            s += &format!("#     folder_name: {}\n", outputconf.folder_name);
            s += &format!("#     file_name: {}\n", outputconf.file_name);
            s += &format!("#     precision: {}\n", outputconf.precision);
            s += &format!(
                "#     should_print_ghostcells: {}\n",
                outputconf.should_print_ghostcells
            );
            s += &format!(
                "#     should_print_metadata: {}\n",
                outputconf.should_print_metadata
            );
            s += "#     data_names: [\n";
            s += &format!("#       {:?}\n", outputconf.data_names);
            s += "#     ]\n";
        }
        s += "###\n";
        return s;
    }
}
