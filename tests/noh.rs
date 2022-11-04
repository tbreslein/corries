// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use corries::config::meshconfig::{MeshConfig, MeshMode};
use corries::config::numericsconfig::{NumFluxMode, NumericsConfig, RkfConfig, TimeIntegrationConfig};
use corries::config::outputconfig::{DataName, FormatterMode, OutputConfig, StreamMode, ToStringConversionMode};
use corries::config::physicsconfig::{PhysicsConfig, PhysicsMode};
use corries::config::{BoundaryMode, CustomBoundaryMode, InitialConditions, PhysicsVariable};
use corries::run_sim;
use corries::units::UnitsMode;
use corries::{config, get_n_equations};

const SIZE: usize = 1000;

fn get_config(mode: PhysicsMode) -> config::CorriesConfig {
    let boundary_conditions_west = match mode {
        PhysicsMode::Euler1DAdiabatic => BoundaryMode::Custom(vec![
            (PhysicsVariable::Density, CustomBoundaryMode::NoGradients),
            (PhysicsVariable::XiVelocity, CustomBoundaryMode::NoGradients),
            (PhysicsVariable::Pressure, CustomBoundaryMode::NoGradients),
        ]),
        PhysicsMode::Euler1DIsot => BoundaryMode::Custom(vec![
            (PhysicsVariable::Density, CustomBoundaryMode::NoGradients),
            (PhysicsVariable::XiVelocity, CustomBoundaryMode::NoGradients),
        ]),
        PhysicsMode::Euler2DIsot => BoundaryMode::Custom(vec![
            (PhysicsVariable::Density, CustomBoundaryMode::NoGradients),
            (PhysicsVariable::XiVelocity, CustomBoundaryMode::NoGradients),
            (PhysicsVariable::EtaVelocity, CustomBoundaryMode::NoGradients),
        ]),
    };

    let boundary_conditions_east = match mode {
        PhysicsMode::Euler1DAdiabatic => BoundaryMode::Custom(vec![
            (PhysicsVariable::Density, CustomBoundaryMode::NoGradients),
            (PhysicsVariable::XiVelocity, CustomBoundaryMode::NoGradients),
            (PhysicsVariable::Pressure, CustomBoundaryMode::NoGradients),
        ]),
        PhysicsMode::Euler1DIsot => BoundaryMode::Custom(vec![
            (PhysicsVariable::Density, CustomBoundaryMode::NoGradients),
            (PhysicsVariable::XiVelocity, CustomBoundaryMode::NoGradients),
        ]),
        PhysicsMode::Euler2DIsot => BoundaryMode::Custom(vec![
            (PhysicsVariable::Density, CustomBoundaryMode::NoGradients),
            (PhysicsVariable::XiVelocity, CustomBoundaryMode::NoGradients),
            (PhysicsVariable::EtaVelocity, CustomBoundaryMode::NoGradients),
        ]),
    };

    let file_name = "noh_".to_owned() + match mode {
        PhysicsMode::Euler1DAdiabatic => "euler1d_adiabatic",
        PhysicsMode::Euler1DIsot => "euler1d_isot",
        PhysicsMode::Euler2DIsot => "euler2d_isot",
    };
    let folder_name = "results/integrationtests/noh_".to_owned() + &file_name;

    return config::CorriesConfig {
        print_banner: true,
        meshconfig: MeshConfig {
            mode: MeshMode::Cartesian,
            xi_in: 1.0,
            xi_out: 2.0,
            ratio_disk: 1.0,
        },
        physicsconfig: PhysicsConfig {
            mode,
            units_mode: UnitsMode::SI,
            adiabatic_index: 1.4,
            c_sound_0: 1.0,
        },
        initial_conditions: InitialConditions::Noh,
        boundary_condition_west: boundary_conditions_west,
        boundary_condition_east: boundary_conditions_east,
        numericsconfig: NumericsConfig {
            numflux_mode: NumFluxMode::Hll,
            time_integration_config: TimeIntegrationConfig::Rkf(RkfConfig {
                rkf_mode: config::numericsconfig::RKFMode::SSPRK5,
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
            dt_cfl_param: 0.4,
        },
        output_counter_max: 1,
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
                formatter_mode: FormatterMode::CSV,
                string_conversion_mode: ToStringConversionMode::Vector,
                folder_name,
                should_clear_out_folder: true,
                file_name,
                precision: 7,
                should_print_ghostcells: true,
                should_print_metadata: false,
                data_names: vec![DataName::XiCent],
            },
            ],
    };
}

#[test]
fn noh_euler1d_adiabatic() {
    const PHYSICS_MODE: PhysicsMode = PhysicsMode::Euler1DAdiabatic;
    const N_EQUATIONS: usize = get_n_equations(PHYSICS_MODE);
    assert!(run_sim::<SIZE, N_EQUATIONS>(get_config(PHYSICS_MODE)).is_ok());
}

#[test]
fn noh_euler1d_isot() {
    const PHYSICS_MODE: PhysicsMode = PhysicsMode::Euler1DIsot;
    const N_EQUATIONS: usize = get_n_equations(PHYSICS_MODE);
    assert!(run_sim::<SIZE, N_EQUATIONS>(get_config(PHYSICS_MODE)).is_ok());
}

#[test]
fn noh_euler2d_isot() {
    const PHYSICS_MODE: PhysicsMode = PhysicsMode::Euler2DIsot;
    const N_EQUATIONS: usize = get_n_equations(PHYSICS_MODE);
    assert!(run_sim::<SIZE, N_EQUATIONS>(get_config(PHYSICS_MODE)).is_ok());
}
