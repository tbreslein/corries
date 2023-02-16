// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::{eyre::Context, Result};
use corries::prelude::*;
use ndarray::Array2;
const S: usize = 100;

fn get_config(numflux_config: NumFluxConfig) -> CorriesConfig {
    let boundary_conditions_west = BoundaryMode::Custom(vec![
        (0, CustomBoundaryMode::NoGradients),
        (1, CustomBoundaryMode::NoGradients),
        (2, CustomBoundaryMode::NoGradients),
    ]);

    let boundary_conditions_east = BoundaryMode::Custom(vec![
        (0, CustomBoundaryMode::NoGradients),
        (1, CustomBoundaryMode::NoGradients),
        (2, CustomBoundaryMode::NoGradients),
    ]);

    let file_name = "sod_".to_owned()
        + match numflux_config {
            NumFluxConfig::Hll => "hll",
            NumFluxConfig::Kt { limiter_mode: _ } => "kt",
        };
    let folder_name = "results/integrationtests/".to_owned() + &file_name;
    let data_names_vector = vec![
        DataName::XiCent,
        DataName::T,
        DataName::Prim(0),
        DataName::Prim(1),
        DataName::Prim(2),
    ];

    return CorriesConfig {
        print_banner: false,
        mesh_config: MeshConfig {
            mode: MeshMode::Cartesian,
            xi_in: 1.0,
            xi_out: 2.0,
            ratio_disk: 1.0,
        },
        physics_config: PhysicsConfig {
            units_mode: UnitsMode::SI,
            adiabatic_index: 1.4,
        },
        boundary_condition_west: boundary_conditions_west,
        boundary_condition_east: boundary_conditions_east,
        numerics_config: NumericsConfig {
            numflux_config,
            time_integration_config: TimeIntegrationConfig::Rkf(RkfConfig {
                rkf_mode: RKFMode::SSPRK5,
                asc: false,
                asc_relative_tolerance: 0.001,
                asc_absolute_tolerance: 0.001,
                asc_timestep_friction: 0.08,
            }),
            iter_max: usize::MAX - 2,
            t0: 0.0,
            t_end: 0.25,
            dt_min: 1.0e-12,
            dt_max: f64::MAX,
            dt_cfl_param: 0.4,
        },
        output_counter_max: 1,
        writer_config: vec![
            OutputConfig {
                stream_mode: StreamMode::Stdout,
                formatting_mode: FormattingMode::TSV,
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
                formatting_mode: FormattingMode::CSV,
                string_conversion_mode: ToStringConversionMode::Vector,
                folder_name,
                should_clear_out_folder: true,
                file_name,
                precision: 7,
                should_print_ghostcells: true,
                should_print_metadata: false,
                data_names: data_names_vector,
            },
        ],
    };
}

fn init<P: Physics<E, S>, const E: usize, const S: usize>(u: &mut State<P, E, S>) {
    let breakpoint_index = (S as f64 * 0.5) as usize;
    let mut prim = Array2::zeros((E, S));
    for i in 0..breakpoint_index {
        prim[[0, i]] = 1.0;
        prim[[E - 1, i]] = 1.0;
    }
    for i in breakpoint_index..S {
        prim[[0, i]] = 0.125;
        prim[[E - 1, i]] = 0.1;
    }
    u.cent.prim.assign(&prim.view());
    return;
}

#[test]
fn sod_hll() -> Result<()> {
    set_Physics_and_E!(Euler1DAdiabatic);
    type N = Hll<E, S>;
    type T = RungeKuttaFehlberg<P, E, S>;

    let config = get_config(NumFluxConfig::Hll);
    let mesh: Mesh<S> = Mesh::new(&config.mesh_config).context("Constructing Mesh")?;
    let mut u = State::<P, E, S>::new(&config.physics_config);
    let mut rhs: Rhs<N, E, S> = Rhs::<N, E, S>::new(&config, &mesh)?;
    let mut time: Time<P, T, E, S> = Time::new(&config)?;
    let mut writer = Writer::new::<S>(&config, &mesh)?;

    init::<P, E, S>(&mut u);
    u.update_vars_from_prim(&mut rhs.boundary_west, &mut rhs.boundary_east, &mesh);

    run_corries::<P, N, T, E, S>(&mut u, &mut rhs, &mut time, &mesh, &mut writer)
        .context("Calling run_loop in noh test")?;
    return Ok(());
}

#[test]
fn sod_kt() -> Result<()> {
    set_Physics_and_E!(Euler1DAdiabatic);
    type N = Kt<E, S>;
    type T = RungeKuttaFehlberg<P, E, S>;

    let config = get_config(NumFluxConfig::Kt {
        limiter_mode: LimiterMode::Monocent(1.2),
    });
    let mesh: Mesh<S> = Mesh::new(&config.mesh_config).context("Constructing Mesh")?;
    let mut u = State::<P, E, S>::new(&config.physics_config);
    let mut rhs: Rhs<N, E, S> = Rhs::<N, E, S>::new(&config, &mesh)?;
    let mut time: Time<P, T, E, S> = Time::new(&config)?;
    let mut writer = Writer::new::<S>(&config, &mesh)?;

    init::<P, E, S>(&mut u);
    u.update_vars_from_prim(&mut rhs.boundary_west, &mut rhs.boundary_east, &mesh);
    u.init_west_east();

    run_corries::<P, N, T, E, S>(&mut u, &mut rhs, &mut time, &mesh, &mut writer)
        .context("Calling run_loop in noh test")?;
    return Ok(());
}
