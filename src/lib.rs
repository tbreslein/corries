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
use physics::Physics;
use rhs::numflux::init_numflux;
use timeintegration::TimeIntegration;
use writer::Writer;

use crate::{errorhandling::Validation, rhs::Rhs};

mod boundaryconditions;
pub mod config;
#[macro_use]
mod errorhandling;
mod mesh;
pub mod physics;
mod rhs;
mod timeintegration;
pub mod units;
mod writer;

/// Runs a Corries simulation.
///
/// # Arguments
///
/// * `config` - The [CorriesConfig] the simulation is based on
pub fn run_sim<const S: usize, const EQ: usize>(config: CorriesConfig) -> Result<()> {
    config.validate().context("Validating config")?;
    let mesh: Mesh<S> = Mesh::new(&config.meshconfig).context("Constructing Mesh")?;
    let mut u: Physics<S, EQ> = Physics::new(&config.physicsconfig);
    let mut numflux = init_numflux(&config.numericsconfig);
    let mut rhs: Rhs<S, EQ> = Rhs::new(&config, &u, &mut numflux);
    let mut timeintegration: TimeIntegration<S, EQ> = TimeIntegration::new(&config);

    let mut writer = Writer::new(&config, &mesh).context("Constructing Writer")?;

    // first output
    if config.print_banner {
        print_banner();
    }
    writer
        .write_metadata::<S>(&config)
        .context("Calling writer.write_metadata in run_sim")?;
    loop {
        if timeintegration.time.t >= timeintegration.time.t_next_output {
            timeintegration.time.t_next_output += timeintegration.time.dt_output;

            writer
                .update_data_matrices(&mesh, &u)
                .context("Calling writer.update_data_matrices in run_sim")?;
            // thread this call?
            writer
                .write_output()
                .context("Calling wirter.write_output in run_sim")?;
        }
        if timeintegration.time.t >= timeintegration.time.t_end {
            break;
        }

        // TODO: next solution
        if let err @ Err(_) = timeintegration.next_solution(&mut u, &mut rhs, &mesh) {
            writer
                .update_data_matrices(&mesh, &u)
                .context("Calling writer.update_data_matrices during the error dump in run_sim")?;
            writer
                .write_output()
                .context("Calling writer.write_output during the error dump in run_sim")?;
            return err;
        }
    }
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
