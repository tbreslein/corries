// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::{eyre::Context, Result};
use corries::prelude::*;
const S: usize = 100;

fn get_config<N: NumFlux<E, S> + 'static, const E: usize>(folder_name: &str, file_name: &str) -> CorriesConfig {
    return CorriesConfig {
        print_banner: false,
        mesh_config: MeshConfig::default_riemann_test(),
        physics_config: PhysicsConfig {
            units_mode: UnitsMode::SI,
            adiabatic_index: 1.4,
        },
        boundary_condition_west: BoundaryMode::NoGradients,
        boundary_condition_east: BoundaryMode::NoGradients,
        numerics_config: NumericsConfig::default_riemann_test::<N, E, S>(0.25),
        output_counter_max: 1,
        writer_config: vec![
            OutputConfig::default_stdout(),
            OutputConfig::default_file(folder_name, file_name, E),
        ],
    };
}

fn init<P, N, T, const E: usize, const S: usize>(
    u: &mut State<P, E, S>,
    _: &mut Solver<P, N, T, E, S>,
    _: &Mesh<S>,
) -> Result<()>
where
    P: Physics<E, S>,
    N: NumFlux<E, S>,
    T: TimeSolver<P, E, S>,
{
    let breakpoint_index = (S as f64 * 0.5) as usize;
    u.cent.prim.fill(0.0);
    for i in 0..breakpoint_index {
        u.cent.prim[[P::JRHO, i]] = 1.0;
        u.cent.prim[[P::JPRESSURE, i]] = 1.0;
    }
    for i in breakpoint_index..S {
        u.cent.prim[[P::JRHO, i]] = 0.125;
        u.cent.prim[[P::JPRESSURE, i]] = 0.1;
    }
    Ok(())
}

#[test]
fn sod_hll() -> Result<()> {
    set_Physics_and_E!(Euler1DAdiabatic);
    type N = Hll<E, S>;
    type T = RungeKuttaFehlberg<P, E, S>;

    get_config::<N, E>("results/integrationtests/sod_hll", "sod_hll")
        .init_corries::<P, N, T, E, S>(init)
        .context("While calling CorriesConfig::init_corries")?
        .run_corries()
}

#[test]
fn sod_kt() -> Result<()> {
    set_Physics_and_E!(Euler1DAdiabatic);
    type N = Kt<E, S>;
    type T = RungeKuttaFehlberg<P, E, S>;

    get_config::<N, E>("results/integrationtests/sod_kt", "sod_kt")
        .init_corries::<P, N, T, E, S>(init)
        .context("While calling CorriesConfig::init_corries")?
        .run_corries()
}
