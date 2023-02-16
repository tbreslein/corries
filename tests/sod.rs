// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::{eyre::Context, Result};
use corries::prelude::*;
use ndarray::Array2;
const S: usize = 100;

fn get_config<N: NumFlux<E,S> + 'static, const E: usize>(folder_name: &str, file_name: &str) -> CorriesConfig {
    return CorriesConfig {
        print_banner: false,
        mesh_config: MeshConfig::default_riemann_test(),
        physics_config: PhysicsConfig {
            units_mode: UnitsMode::SI,
            adiabatic_index: 1.4,
        },
        boundary_condition_west: BoundaryMode::NoGradients,
        boundary_condition_east: BoundaryMode::NoGradients,
        numerics_config: NumericsConfig::default_riemann_test::<N,E,S>(0.25),
        output_counter_max: 1,
        writer_config: vec![
            OutputConfig::default_stdout(),
            OutputConfig::default_file(folder_name, file_name, E),
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

    let config = get_config::<N,E>("results/integrationtests/sod_hll", "sod_hll");
    let (mut u, mut rhs, mut time, mesh, mut writer) = init_corries::<P, N, T, E, S>(&config).unwrap();

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

    let config = get_config::<N,E>("results/integrationtests/sod_kt", "sod_kt");
    let (mut u, mut rhs, mut time, mesh, mut writer) = init_corries::<P, N, T, E, S>(&config).unwrap();

    init::<P, E, S>(&mut u);
    u.update_vars_from_prim(&mut rhs.boundary_west, &mut rhs.boundary_east, &mesh);
    u.init_west_east();

    run_corries::<P, N, T, E, S>(&mut u, &mut rhs, &mut time, &mesh, &mut writer)
        .context("Calling run_loop in noh test")?;
    return Ok(());
}
