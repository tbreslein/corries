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

const CONFIG: CorriesConfig = CorriesConfig {
    print_banner: true,
    mesh_config: MeshConfig {
        mode: MeshMode::Cartesian,
        xi_in: 1.0,
        xi_out: 2.0,
        ratio_disk: 1.0,
    },
    physics_config: PhysicsConfig { adiabatic_index: 1.4 }
};

fn main() -> Result<()> {
    let mesh: Mesh<S> = Mesh::new(&CONFIG.mesh_config).context("Constructing Mesh")?;
    let mut u: P = P::new(&CONFIG.physics_config);
    run_corries(&mut u, &mesh)?;
    return Ok(());
}
