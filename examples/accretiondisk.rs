// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::Result;
use corries::config;
use corries::config::meshconfig::{MeshConfig, MeshMode};
use corries::config::outputconfig::{DataName, FormatterMode, OutputConfig, StreamMode, ToStringConversionMode};
use corries::config::physicsconfig::{PhysicsConfig, PhysicsMode, UnitsMode};
use corries::run_sim;

fn main() -> Result<()> {
    run_sim(config::CorriesConfig {
        name: "accretiondisk".to_string(),
        meshconf: MeshConfig {
            mode: MeshMode::Cartesian,
            n_comp: 10,
            n_gc: 2,
            xi_in: 0.1,
            xi_out: 10.0,
            ratio_disk: 1.0,
        },
        physicsconf: PhysicsConfig {
            mode: PhysicsMode::Euler1DIsot,
            adiabatic_index: 2.0 / 3.0,
            c_sound_0: 1.0,
        },
        writerconf: vec![
            OutputConfig {
                stream_mode: StreamMode::Stdout,
                formatter_mode: FormatterMode::TSV,
                string_conversion_mode: ToStringConversionMode::Scalar,
                folder_name: "".to_string(),
                should_clear_out_folder: false,
                file_name: "".to_string(),
                precision: 3,
                should_print_ghostcells: false,
                should_print_metadata: true,
                data_names: vec![DataName::NComp, DataName::NAll],
            },
            OutputConfig {
                stream_mode: StreamMode::File,
                formatter_mode: FormatterMode::TSV,
                string_conversion_mode: ToStringConversionMode::Vector,
                folder_name: "results/accretiondisk".to_string(),
                should_clear_out_folder: true,
                file_name: "accretiondisk".to_string(),
                precision: 7,
                should_print_ghostcells: true,
                should_print_metadata: false,
                data_names: vec![DataName::NComp, DataName::NAll, DataName::XiCent],
            },
        ],
    })
}
