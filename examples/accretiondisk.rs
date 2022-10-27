// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::Result;
use corries::config;
use corries::config::meshconfig::MeshMode;
use corries::config::outputconfig::{DataName, FormatterMode, StreamMode, ToStringConversionMode};
use corries::run_sim;

fn main() -> Result<()> {
    run_sim(config::CorriesConfig {
        name: "accretiondisk".to_string(),
        meshconf: config::meshconfig::MeshConfig {
            mode: MeshMode::Cartesian,
            n_comp: 10,
            n_gc: 2,
            xi_in: 0.1,
            xi_out: 100.0,
            ratio_disk: 1.0,
        },
        writerconf: vec![
            config::outputconfig::OutputConfig {
                stream_mode: StreamMode::Stdout,
                formatter_mode: FormatterMode::CSV,
                string_conversion_mode: ToStringConversionMode::Vector,
                folder_name: "".to_string(),
                file_name: "".to_string(),
                precision: 3,
                should_print_ghostcells: false,
                should_print_metadata: true,
                data: vec![DataName::NComp, DataName::XiCent],
            },
            // config::outputconfig::OutputConfig {
            //     stream_mode: StreamMode::File,
            //     formatter_mode: FormatterMode::CSV,
            //     string_conversion_mode: ToStringConversionMode::Vector,
            //     folder_name: "results".to_string(),
            //     file_name: "accretiondisk".to_string(),
            //     precision: 7,
            //     should_print_ghostcells: true,
            //     should_print_metadata: false,
            //     data: vec![DataName::XiCent],
            // },
        ],
    })
}
