// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

// #![warn(missing_docs)]

//! Corries - CORrosive RIEman Solver
//!
//! Library to run 1D-hydrodynamics simulations solved with Riemann solvers.
//!
//! Most importantly exports the module [prelude].

use color_eyre::{eyre::Context, Result};

mod boundaryconditions;
pub mod config;
#[macro_use]
mod errorhandling;
#[macro_use]
pub mod macros;
pub mod mesh;
pub mod physics;
pub mod rhs;
pub mod time;
pub mod units;
pub mod writer;

/// Exports everything you need to run a corries simulation. This includes the following modules
///
/// - [config]
/// - [mesh]
/// - [physics]
/// - [rhs]
/// - [time]
/// - [units]
/// - [writer]
///
/// as well as the [set_Physics_and_E] macro, and the [run_corries] function
pub mod prelude {
    pub use crate::config::*;
    pub use crate::init_corries;
    pub use crate::mesh::*;
    pub use crate::physics::*;
    pub use crate::rhs::*;
    pub use crate::run_corries;
    pub use crate::set_Physics_and_E;
    pub use crate::time::*;
    pub use crate::units::*;
    pub use crate::writer::*;
}

use errorhandling::Validation;
pub use prelude::*;

/// Initialises all objects needed to run a corries simulation.
///
/// Apart from the `config` argument, the important bits that also help configuring the simulation
/// are the template Parameters. For example, the type you pass as the first template argument
/// determines the type of `Physics` used throughout the whole simulation!
pub fn init_corries<P, N, T, const E: usize, const S: usize>(
    config: &CorriesConfig,
) -> Result<(P, Rhs<P, N, S>, Time<P, T>, Mesh<S>, Writer)>
where
    P: Physics + Collectable + 'static,
    N: NumFlux,
    T: TimeSolver<P>,
{
    if cfg!(feature = "validation") {
        config.validate()?;
    }

    let mesh: Mesh<S> = Mesh::new(&config.mesh_config).context("Constructing Mesh")?;
    let u: P = P::new(&config.physics_config);
    let rhs: Rhs<P, N, S> = Rhs::<P, N, S>::new::<E>(config);
    let time: Time<P, T> = Time::new::<E, S>(config, &u)?;
    let mut writer = Writer::new::<S>(config, &mesh)?;

    if writer.print_banner {
        print_banner();
    }
    writer
        .write_metadata::<S>()
        .context("Calling writer.write_metadata in run_corries")?;
    return Ok((u, rhs, time, mesh, writer));
}

/// Runs a corries simulation.
///
/// This assumes that you already set your initial conditions in `u` and that it is fully
/// up-to-date.
/// During the run, [writer] will be writing output according to your specifications when you wrote
/// the [CorriesConfig] struct to configure everthing.
///
/// # Arguments
///
/// * `u` - Carries with the state of the physical simulation
/// * `rhs - Solves the right-hand side of the equations
/// * `time - Solves the time integration step and carries information about the time coordinate
/// and the time step
/// * `mesh` - The mesh the simulation runs on
/// * `writer` - Deals with writing output
pub fn run_corries<P: Physics + Collectable, N: NumFlux, T: TimeSolver<P>, const E: usize, const S: usize>(
    u: &mut P,
    rhs: &mut Rhs<P, N, S>,
    time: &mut Time<P, T>,
    mesh: &Mesh<S>,
    writer: &mut Writer,
) -> Result<()> {
    loop {
        if time.timestep.t >= time.timestep.t_next_output - time.timestep.dt_min {
            time.timestep.t_next_output += time.timestep.dt_output;
            writer
                .update_data(u, time, mesh)
                .context("Calling writer.update_data_matrices in run_corries")?;
            writer
                .write_output()
                .context("Calling wirter.write_output in run_corries")?;
        }
        if time.timestep.t + time.timestep.dt_min * 0.01 >= time.timestep.t_end {
            break;
        }

        if let err @ Err(_) = time.next_solution::<N, E, S>(u, rhs, mesh) {
            time.timestep.dt_kind = DtKind::ErrorDump;
            writer
                .update_data(u, time, mesh)
                .context("Calling writer.update_data_matrices during the error dump in run_corries")?;
            writer
                .write_output()
                .context("Calling writer.write_output during the error dump in run_corries")?;
            return err;
        }
    }
    return Ok(());
}

const VERSION: &str = env!("CARGO_PKG_VERSION");
fn print_banner() {
    println!("# ****************************************");
    println!("# Corries - corrosive Riemann solver ");
    println!("# ");
    println!("# Version: {}", VERSION);
    println!("# Copyright (c) 2022-2023");
    println!("# Author: tbreslein <github.com/tbreslein>");
    println!("# License: MIT");
    println!("# ****************************************");
}
