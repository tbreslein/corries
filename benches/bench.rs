// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use corries::prelude::*;
use criterion::{criterion_group, criterion_main, Criterion};
use ndarray::{Array1, Array2};

const S: usize = 500;

const EULER1D_ADIABATIC: usize = 0;
const EULER1D_ISOT: usize = 1;
const EULER2D_ISOT: usize = 2;

fn get_config(mode: usize) -> CorriesConfig {
    let boundary_conditions_west = match mode {
        EULER1D_ADIABATIC => BoundaryMode::Custom(vec![
            (0, CustomBoundaryMode::NoGradients),
            (1, CustomBoundaryMode::NoGradients),
            (2, CustomBoundaryMode::NoGradients),
        ]),
        EULER1D_ISOT => BoundaryMode::Custom(vec![
            (0, CustomBoundaryMode::NoGradients),
            (1, CustomBoundaryMode::NoGradients),
        ]),
        EULER2D_ISOT => BoundaryMode::Custom(vec![
            (0, CustomBoundaryMode::NoGradients),
            (1, CustomBoundaryMode::NoGradients),
            (2, CustomBoundaryMode::NoGradients),
        ]),
        _ => BoundaryMode::Custom(vec![]),
    };

    let boundary_conditions_east = match mode {
        EULER1D_ADIABATIC => BoundaryMode::Custom(vec![
            (0, CustomBoundaryMode::NoGradients),
            (1, CustomBoundaryMode::NoGradients),
            (2, CustomBoundaryMode::NoGradients),
        ]),
        EULER1D_ISOT => BoundaryMode::Custom(vec![
            (0, CustomBoundaryMode::NoGradients),
            (1, CustomBoundaryMode::NoGradients),
        ]),
        EULER2D_ISOT => BoundaryMode::Custom(vec![
            (0, CustomBoundaryMode::NoGradients),
            (1, CustomBoundaryMode::NoGradients),
            (2, CustomBoundaryMode::NoGradients),
        ]),
        _ => BoundaryMode::Custom(vec![]),
    };

    let file_name = "noh_".to_owned()
        + match mode {
            EULER1D_ADIABATIC => "euler1d_adiabatic",
            EULER1D_ISOT => "euler1d_isot",
            EULER2D_ISOT => "euler2d_isot",
            _ => "",
        };
    let folder_name = "results/integrationtests/".to_owned() + &file_name;
    let data_names_vector = match mode {
        EULER1D_ADIABATIC => vec![
            DataName::XiCent,
            DataName::T,
            DataName::Prim(0),
            DataName::Prim(1),
            DataName::Prim(2),
        ],
        EULER1D_ISOT => vec![DataName::XiCent, DataName::T, DataName::Prim(0), DataName::Prim(1)],
        EULER2D_ISOT => vec![
            DataName::XiCent,
            DataName::T,
            DataName::Prim(0),
            DataName::Prim(1),
            DataName::Prim(2),
        ],
        _ => vec![],
    };

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
            adiabatic_index: 5.0 / 3.0,
        },
        boundary_condition_west: boundary_conditions_west,
        boundary_condition_east: boundary_conditions_east,
        numerics_config: NumericsConfig {
            numflux_config: NumFluxConfig::Kt {
                limiter_mode: LimiterMode::VanLeer,
            },
            time_integration_config: TimeIntegrationConfig::Rkf(RkfConfig {
                rkf_mode: RKFMode::SSPRK5,
                asc: false,
                asc_relative_tolerance: 0.001,
                asc_absolute_tolerance: 0.001,
                asc_timestep_friction: 0.08,
            }),
            iter_max: usize::MAX - 2,
            t0: 0.0,
            t_end: 0.5,
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

fn init_noh<P: Physics<E, S>, const E: usize, const S: usize>(u: &mut State<P, E, S>) {
    let breakpoint_index = (S as f64 * 0.5) as usize;
    let mut prim = Array2::zeros((E, S));
    prim.fill(0.0);
    for i in 0..breakpoint_index {
        prim[[0, i]] = 1.0;
        prim[[1, i]] = 1.0;
    }
    for i in breakpoint_index..S {
        prim[[0, i]] = 1.0;
        prim[[1, i]] = -1.0;
    }
    if u.is_adiabatic() {
        prim.row_mut(E - 1).fill(1.0E-5)
    } else {
        let c_sound = Array1::ones(S);
        u.cent.c_sound.assign(&c_sound.view());
    }
    u.cent.prim.assign(&prim.view());
    return;
}

pub fn noh_run(c: &mut Criterion) {
    let mut group = c.benchmark_group("noh_run");

    set_Physics_and_E!(Euler1DAdiabatic);
    let config = get_config(EULER1D_ADIABATIC);
    type N = Hll<E, S>;
    type T = RungeKuttaFehlberg<P, E, S>;

    let (mut u, mut rhs, mut time, mesh, mut writer) = init_corries::<P, N, T, E, S>(&config).unwrap();
    init_noh::<P, E, S>(&mut u);
    u.update_cent_from_prim(&mut rhs.boundary_west, &mut rhs.boundary_east, &mesh);

    group.bench_function("noh test run", |b| {
        b.iter(|| {
            run_corries::<P, N, T, E, S>(&mut u, &mut rhs, &mut time, &mesh, &mut writer).unwrap();
        })
    });
    group.finish();
}

pub fn euler1d_isot_conversions(c: &mut Criterion) {
    const PHYSICS_CONFIG: PhysicsConfig = PhysicsConfig {
        adiabatic_index: 1.4,
        units_mode: UnitsMode::SI,
    };
    let mut group = c.benchmark_group("physics_conversions");
    group.sample_size(1000);

    let mut u = State::<Euler1DIsot<S>, 2, S>::new(&PHYSICS_CONFIG);
    u.cent.prim.row_mut(0).fill(2.0);
    u.cent.prim.row_mut(1).fill(3.0);
    group.bench_function("from prim to cons and back; euler 1d isot", |b| {
        b.iter(|| {
            u.update_cons_cent();
            u.update_prim_cent();
        })
    });

    let mut u = State::<Euler1DAdiabatic<S>, 3, S>::new(&PHYSICS_CONFIG);
    u.cent.prim.row_mut(0).fill(2.0);
    u.cent.prim.row_mut(1).fill(3.0);
    u.cent.prim.row_mut(2).fill(4.0);
    group.bench_function("from prim to cons and back; euler 1d adiabatic", |b| {
        b.iter(|| {
            u.update_cons_cent();
            u.update_prim_cent();
        })
    });

    group.finish();
}

criterion_group!(benches, euler1d_isot_conversions, noh_run);
criterion_main!(benches);
