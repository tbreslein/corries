// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [CorriesConfig] structs and its nested structs for configuring Corries simulations.

use color_eyre::eyre::{ensure, Context};
use color_eyre::Result;

use self::meshconfig::MeshConfig;
use self::numericsconfig::NumericsConfig;
use self::outputconfig::OutputConfig;
use self::physicsconfig::PhysicsConfig;
use crate::errorhandling::Validation;

pub mod meshconfig;
pub mod numericsconfig;
pub mod outputconfig;
pub mod physicsconfig;

/// Struct that carries the full configuration info for a simulation.
///
/// This struct is used in the beginning of a run to initialise all the runtime-objects that are
/// used throughout the simulation.
#[derive(Debug)]
pub struct CorriesConfig {
    /// Name of the simulation; may not be empty
    pub name: String,

    /// Whether to print the Corries banner to stdout
    pub print_banner: bool,

    /// Config for Mesh objects
    pub meshconf: MeshConfig,

    /// Config for Physics objects
    pub physicsconf: PhysicsConfig,

    /// Config for everything related to numerics
    pub numericsconf: NumericsConfig,

    /// Config for Writer objects
    pub writerconf: Vec<OutputConfig>,
}

impl Validation for CorriesConfig {
    fn validate(&self) -> Result<()> {
        ensure!(!self.name.is_empty(), "name must not be empty!");
        self.meshconf.validate().context("Validating config.meshconf")?;
        for outputconf in self.writerconf.iter() {
            outputconf.validate().context("Validating config.writerconf")?;
        }
        return Ok(());
    }
}

impl CorriesConfig {
    /// Writes metadata, mostly just the contents of this struct and its nested structs, into a
    /// `String`.
    pub fn metadata_dump<const MESH_COMP_SIZE: usize>(&self) -> String {
        let mut s = "".to_string();
        s += "### Corries Configuration\n";
        s += &format!("# name: {}\n", self.name);

        s += "# meshconf:\n";
        s += &format!("#     mode: {:?}\n", self.meshconf.mode);
        s += &format!("#     n_comp: {}\n", MESH_COMP_SIZE);
        // s += &format!("#     n_gc: {}\n", NGC);
        s += &format!("#     xi_in: {}\n", self.meshconf.xi_in);
        s += &format!("#     xi_out: {}\n", self.meshconf.xi_out);
        s += &format!("#     ratio_disk: {}\n", self.meshconf.ratio_disk);

        s += "# physicsconf:\n";
        s += &format!("#     mode: {:?}\n", self.physicsconf.mode);
        s += &format!("#     units_mode: {:?}\n", self.physicsconf.units_mode);
        s += &format!("#     adiabatic_index: {}\n", self.physicsconf.adiabatic_index);
        s += &format!("#     c_sound_0: {}\n", self.physicsconf.c_sound_0);

        s += "# numericsconf:\n";
        s += &format!("#     numflux_mode: {:?}\n", self.numericsconf.numflux_mode);

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
            s += &format!(
                "#     should_clear_out_folder: {}\n",
                outputconf.should_clear_out_folder
            );
            s += &format!("#     file_name: {}\n", outputconf.file_name);
            s += &format!("#     precision: {}\n", outputconf.precision);
            s += &format!(
                "#     should_print_ghostcells: {}\n",
                outputconf.should_print_ghostcells
            );
            s += &format!("#     should_print_metadata: {}\n", outputconf.should_print_metadata);
            s += "#     data_names: [\n";
            s += &format!("#       {:?}\n", outputconf.data_names);
            s += "#     ]\n";
        }
        s += "###\n";
        return s;
    }
}
