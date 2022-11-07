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

/// Initialises and returns the objects needed for running a corries simulation.
///
/// # Arguments
///
/// ## Generic arguments
///
/// * `S` - Size of the mesh, as in the number of cells
/// * `EQ` - The number of equations in the type of physics
///
/// ## Function arguments
///
/// * `config` - The [CorriesConfig] the simulation is based on
pub fn init_sim<const S: usize, const EQ: usize>(
    config: &CorriesConfig,
) -> Result<(Physics<S, EQ>, Rhs<S, EQ>, TimeIntegration<S, EQ>, Mesh<S>, Writer)> {
    config.validate().context("Validating config")?;
    let mesh: Mesh<S> = Mesh::new(&config.meshconfig).context("Constructing Mesh")?;
    let u: Physics<S, EQ> = Physics::new(&config.physicsconfig);
    let rhs: Rhs<S, EQ> = Rhs::new(config, &u);
    let timeintegration: TimeIntegration<S, EQ> =
        TimeIntegration::new(config, &u).context("Constructing TimeIntegration")?;
    let mut writer = Writer::new(config, &mesh).context("Constructing Writer")?;

    if config.print_banner {
        print_banner();
    }
    writer
        .write_metadata::<S>(config)
        .context("Calling writer.write_metadata in run_sim")?;
    return Ok((u, rhs, timeintegration, mesh, writer));
}

/// Runs the core loop of a corries simulation
///
/// # Arguments
///
/// * `u` - The [Physics] the simulation is based on
/// * `rhs` - The [Rhs] that solves the right-hand side
/// * `timeintegration` - The [TimeIntegration] that solves the time integration step
/// * `mesh` - The [Mesh] the simulation runs on
/// * `writer` - The [Writer] that writes output
pub fn run_loop<const S: usize, const EQ: usize>(
    u: &mut Physics<S, EQ>,
    rhs: &mut Rhs<S, EQ>,
    timeintegration: &mut TimeIntegration<S, EQ>,
    mesh: &Mesh<S>,
    writer: &mut Writer,
) -> Result<()> {
    loop {
        if timeintegration.time.t >= timeintegration.time.t_next_output {
            timeintegration.time.t_next_output += timeintegration.time.dt_output;

            writer
                .update_data_matrices(mesh, u, timeintegration)
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
        if let err @ Err(_) = timeintegration.next_solution(u, rhs, mesh) {
            writer
                .update_data_matrices(mesh, u, timeintegration)
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
