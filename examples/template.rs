// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::Result;
use corries::config::meshconfig::{MeshConfig, MeshMode};
use corries::config::numericsconfig::{NumFluxMode, NumericsConfig, RkfConfig, TimeIntegrationConfig};
use corries::config::outputconfig::{DataName, FormatterMode, OutputConfig, StreamMode, ToStringConversionMode};
use corries::config::physicsconfig::{PhysicsConfig, PhysicsMode};
use corries::config::{BoundaryMode, CustomBoundaryMode, PhysicsVariable};
use corries::units::UnitsMode;
use corries::{config, get_n_equations};
use corries::{init_sim, run_loop};

fn main() -> Result<()> {
    const SIZE: usize = 14;
    const N_EQUATIONS: usize = get_n_equations(PhysicsMode::Euler1DIsot);

    let (mut u, mut rhs, mut timeintegration, mesh, mut writer) =
        init_sim::<SIZE, N_EQUATIONS>(&config::CorriesConfig {
            print_banner: true,
            meshconfig: MeshConfig {
                mode: MeshMode::Cartesian,
                xi_in: 0.1,
                xi_out: 10.0,
                ratio_disk: 1.0,
            },
            physicsconfig: PhysicsConfig {
                mode: PhysicsMode::Euler1DIsot,
                units_mode: UnitsMode::SI,
                adiabatic_index: 1.4,
                c_sound_0: 1.0,
            },
            boundary_condition_west: BoundaryMode::Custom(vec![
                (PhysicsVariable::Density, CustomBoundaryMode::NoGradients),
                (PhysicsVariable::XiVelocity, CustomBoundaryMode::NoGradients),
            ]),
            boundary_condition_east: BoundaryMode::Custom(vec![
                (PhysicsVariable::Density, CustomBoundaryMode::NoGradients),
                (PhysicsVariable::XiVelocity, CustomBoundaryMode::NoGradients),
            ]),
            numericsconfig: NumericsConfig {
                numflux_mode: NumFluxMode::Hll,
                time_integration_config: TimeIntegrationConfig::Rkf(RkfConfig {
                    rkf_mode: config::numericsconfig::RKFMode::RK4,
                    asc: false,
                    asc_relative_tolerance: 0.001,
                    asc_absolute_tolerance: 0.1,
                    asc_timestep_friction: 0.08,
                }),
                iter_max: usize::MAX - 2,
                t0: 0.0,
                t_end: 10.0,
                dt_min: 1.0e-12,
                dt_max: f64::MAX,
                dt_cfl_param: 0.1,
            },
            output_counter_max: 10,
            writerconfig: vec![
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
        })?;

    return run_loop(&mut u, &mut rhs, &mut timeintegration, &mesh, &mut writer);
}
