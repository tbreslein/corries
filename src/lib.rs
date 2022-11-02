// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

#![warn(missing_docs)]

//! Corries - CORrosive RIEman Solver
//!
//! Library to run 1D-hydrodynamics simulations solved with Riemann solvers.
//!
//! Provides the [run_sim] and [get_n_equations] functions, as well as the [config::CorriesConfig] struct.
//! The [run_sim] function is used to start a simulation configured through its input arguments.

use color_eyre::{eyre::Context, Result};
use config::{physicsconfig::PhysicsMode, CorriesConfig};
use mesh::Mesh;
use physics::init_physics;
use rhs::numflux::init_numflux;
use writer::Writer;

use crate::{errorhandling::Validation, rhs::Rhs};

mod boundaryconditions;
pub mod config;
#[macro_use]
mod errorhandling;
mod mesh;
pub mod physics;
mod rhs;
pub mod units;
mod writer;

/// Runs a Corries simulation.
///
/// # Arguments
///
/// * `config` - The [CorriesConfig] the simulation is based on
pub fn run_sim<const S: usize, const EQ: usize>(config: CorriesConfig) -> Result<()> {
    config.validate().context("Validating config")?;
    let mesh: Mesh<S> = Mesh::new(&config.meshconf).context("Constructing Mesh")?;
    let mut u = init_physics::<S, EQ>(&config.physicsconf);
    let mut numflux = init_numflux(&config.numericsconf);
    let mut rhs: Rhs<S, EQ> = Rhs::new(&config, &u, &mut numflux);
    rhs.update_dflux_dxi(&mut u, &mesh);

    // // TEMP:
    let output_count_max = 2;
    let mut writer = Writer::new(&config, &mesh, output_count_max)?;
    //
    // // first output
    if config.print_banner {
        print_banner();
    }
    writer.update_data_matrices(&mesh, &u)?;
    writer.write_metadata::<S>(&config)?;
    writer.write_output()?;
    return Ok(());
}

/// Compile time function to return the number of equations according to the [PhysicsMode]
/// argument.
pub const fn get_n_equations(physics_mode: PhysicsMode) -> usize {
    return match physics_mode {
        PhysicsMode::Euler1DAdiabatic => 3,
        PhysicsMode::Euler1DIsot => 2,
        PhysicsMode::Euler2DIsot => 3,
    };
}

const VERSION: &str = env!("CARGO_PKG_VERSION");
fn print_banner() {
    println!("# ********************************************");
    println!("# *** Corries - corrosive Riemann solver ");
    println!("# *** ");
    println!("# *** Version: {}", VERSION);
    println!("# *** Copyright (c) 2022");
    println!("# *** Author: tbreslein <github.com/tbreslein>");
    println!("# *** License: MIT");
    println!("# ********************************************");
}
