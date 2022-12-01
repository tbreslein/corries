// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::{eyre::Context, Result};
use corries::prelude::*;

// =============
// CONFIGURATION
// =============

const S: usize = 10;
set_Physics_and_E!(Euler1DIsot);
type N = Hll;

fn main() -> Result<()> {
    let config: CorriesConfig = CorriesConfig {
        print_banner: true,
        mesh_config: MeshConfig {
            mode: MeshMode::Cartesian,
            xi_in: 1.0,
            xi_out: 2.0,
            ratio_disk: 1.0,
        },
        physics_config: PhysicsConfig { adiabatic_index: 1.4 },
        boundary_condition_west: BoundaryMode::Custom(vec![
            (0, CustomBoundaryMode::NoGradients),
            (1, CustomBoundaryMode::NoGradients),
        ]),
        boundary_condition_east: BoundaryMode::Custom(vec![
            (0, CustomBoundaryMode::NoGradients),
            (1, CustomBoundaryMode::NoGradients),
        ]),
        numerics_config: NumericsConfig {
            time_integration_config: TimeIntegrationConfig::Rkf(RkfConfig {
                rkf_mode: RKFMode::RK4,
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
    };

    let mesh: Mesh<S> = Mesh::new(&config.mesh_config).context("Constructing Mesh")?;
    let mut u: P = P::new(&config.physics_config);
    let mut rhs: Rhs<P, N, S> = Rhs::<P, N, S>::new::<E>(&config);
    run_corries(&mut u, &mut rhs, &mesh)?;
    return Ok(());
}
