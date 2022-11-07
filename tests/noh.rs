// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::eyre::Context;
use color_eyre::Result;
use corries::config::meshconfig::{MeshConfig, MeshMode};
use corries::config::numericsconfig::{NumFluxMode, NumericsConfig, RkfConfig, TimeIntegrationConfig};
use corries::config::outputconfig::{DataName, FormatterMode, OutputConfig, StreamMode, ToStringConversionMode};
use corries::config::physicsconfig::{PhysicsConfig, PhysicsMode};
use corries::config::{BoundaryMode, CustomBoundaryMode, PhysicsVariable};
use corries::physics::Physics;
use corries::units::UnitsMode;
use corries::{config, get_n_equations};
use corries::{init_sim, run_loop};

const SIZE: usize = 100;

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

    let file_name = "noh_".to_owned()
        + match mode {
            PhysicsMode::Euler1DAdiabatic => "euler1d_adiabatic",
            PhysicsMode::Euler1DIsot => "euler1d_isot",
            PhysicsMode::Euler2DIsot => "euler2d_isot",
        };
    let folder_name = "results/integrationtests/noh_".to_owned() + &file_name;

    let t_end = match mode {
        PhysicsMode::Euler1DAdiabatic => 1.25,
        PhysicsMode::Euler1DIsot | PhysicsMode::Euler2DIsot => 0.5,
    };

    return config::CorriesConfig {
        print_banner: false,
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
        },
        boundary_condition_west: boundary_conditions_west,
        boundary_condition_east: boundary_conditions_east,
        numericsconfig: NumericsConfig {
            numflux_mode: NumFluxMode::Hll,
            time_integration_config: TimeIntegrationConfig::Rkf(RkfConfig {
                rkf_mode: config::numericsconfig::RKFMode::RK1,
                asc: false,
                asc_relative_tolerance: 0.001,
                asc_absolute_tolerance: 0.1,
                asc_timestep_friction: 0.08,
            }),
            iter_max: usize::MAX - 2,
            t0: 0.0,
            t_end,
            dt_min: 1.0e-12,
            dt_max: f64::MAX,
            dt_cfl_param: 0.4,
        },
        output_counter_max: 1,
        writerconfig: vec![
            OutputConfig {
                stream_mode: StreamMode::Stdout,
                formatter_mode: FormatterMode::TSV,
                string_conversion_mode: ToStringConversionMode::Scalar,
                folder_name: "".to_string(),
                should_clear_out_folder: false,
                file_name: "".to_string(),
                precision: 3,
                should_print_ghostcells: false,
                should_print_metadata: false,
                data_names: vec![DataName::Iter, DataName::T, DataName::Dt, DataName::DtKind],
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
                data_names: vec![DataName::XiCent, DataName::Prim(0), DataName::Prim(1)],
            },
        ],
    };
}

fn init_noh<const S: usize, const EQ: usize>(u: &mut Physics<S, EQ>) {
    let breakpoint_index = (S as f64 * 0.5) as usize;
    u.prim.fill(0.0);
    u.cons.fill(0.0);
    for i in 0..breakpoint_index {
        u.prim[[u.jdensity, i]] = 1.0;
        u.prim[[u.jxivelocity, i]] = 1.0;
    }
    for i in breakpoint_index..S {
        u.prim[[u.jdensity, i]] = 1.0;
        u.prim[[u.jxivelocity, i]] = -1.0;
    }
    if u.is_adiabatic {
        u.prim.row_mut(u.jpressure).fill(1.0e-5)
    }
    u.c_sound.fill(1.0);
    return;
}

// #[test]
// fn noh_euler1d_adiabatic() {
//     const PHYSICS_MODE: PhysicsMode = PhysicsMode::Euler1DAdiabatic;
//     const N_EQUATIONS: usize = get_n_equations(PHYSICS_MODE);
//     let (mut u, mut rhs, mut timeintegration, mesh, mut writer) = init_sim::<SIZE, N_EQUATIONS>(&get_config(PHYSICS_MODE)).unwrap();
//     init_noh(&mut u);
//     assert!(run_loop(&mut u, &mut rhs, &mut timeintegration, &mesh, &mut writer).is_ok());
// }

#[test]
fn noh_euler1d_isot() -> Result<()> {
    const PHYSICS_MODE: PhysicsMode = PhysicsMode::Euler1DIsot;
    const N_EQUATIONS: usize = get_n_equations(PHYSICS_MODE);
    let (mut u, mut rhs, mut timeintegration, mesh, mut writer) =
        init_sim::<SIZE, N_EQUATIONS>(&get_config(PHYSICS_MODE)).context("Calling init_sim in noh test")?;
    init_noh(&mut u);
    u.update_everything_from_prim(&mut rhs.boundary_conditions, &mesh)
        .context("Calling u.update_everything_from_prim in noh test")?;
    run_loop(&mut u, &mut rhs, &mut timeintegration, &mesh, &mut writer).context("Calling run_loop in noh test")?;
    return Ok(());
}

// #[test]
// fn noh_euler2d_isot() {
//     const PHYSICS_MODE: PhysicsMode = PhysicsMode::Euler2DIsot;
//     const N_EQUATIONS: usize = get_n_equations(PHYSICS_MODE);
//     let (mut u, mut rhs, mut timeintegration, mesh, mut writer) = init_sim::<SIZE, N_EQUATIONS>(&get_config(PHYSICS_MODE)).unwrap();
//     init_noh(&mut u);
//     assert!(run_loop(&mut u, &mut rhs, &mut timeintegration, &mesh, &mut writer).is_ok());
// }
