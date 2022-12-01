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
    };

    let mesh: Mesh<S> = Mesh::new(&config.mesh_config).context("Constructing Mesh")?;
    let mut u: P = P::new(&config.physics_config);
    let mut rhs: Rhs<P, S> = Rhs::<P, S>::new::<E>(&config);
    run_corries(&mut u, &mut rhs, &mesh)?;
    return Ok(());
}
