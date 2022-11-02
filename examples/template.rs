// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::Result;
use corries::config::meshconfig::{MeshConfig, MeshMode};
use corries::config::numericsconfig::{NumFluxMode, NumericsConfig};
use corries::config::outputconfig::{DataName, FormatterMode, OutputConfig, StreamMode, ToStringConversionMode};
use corries::config::physicsconfig::{PhysicsConfig, PhysicsMode};
use corries::run_sim;
use corries::units::UnitsMode;
use corries::{config, get_n_equations};

fn main() -> Result<()> {
    const SIZE: usize = 14;
    const N_EQUATIONS: usize = get_n_equations(PhysicsMode::Euler1DIsot);

    run_sim::<SIZE, N_EQUATIONS>(config::CorriesConfig {
        name: "accretiondisk".to_string(),
        print_banner: true,
        meshconf: MeshConfig {
            mode: MeshMode::Cartesian,
            xi_in: 0.1,
            xi_out: 10.0,
            ratio_disk: 1.0,
        },
        physicsconf: PhysicsConfig {
            mode: PhysicsMode::Euler1DIsot,
            units_mode: UnitsMode::SI,
            adiabatic_index: 2.0 / 3.0,
            c_sound_0: 1.0,
        },
        numericsconf: NumericsConfig {
            numflux_mode: NumFluxMode::Hll,
        },
        writerconf: vec![
            OutputConfig {
                stream_mode: StreamMode::Stdout,
                formatter_mode: FormatterMode::TSV,
                string_conversion_mode: ToStringConversionMode::Vector,
                folder_name: "".to_string(),
                should_clear_out_folder: false,
                file_name: "".to_string(),
                precision: 3,
                should_print_ghostcells: false,
                should_print_metadata: true,
                data_names: vec![DataName::Prim(0), DataName::Cons(1)],
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
                data_names: vec![DataName::XiCent],
            },
        ],
    })
}
