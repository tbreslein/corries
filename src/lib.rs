// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

// #![warn(missing_docs)]

//! Corries - CORrosive RIEman Solver
//!
//! Library to run 1D-hydrodynamics simulations solved with Riemann solvers.
//!
//! TODO

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

/// TODO
pub mod prelude {
    pub use crate::config::*;
    pub use crate::mesh::*;
    pub use crate::physics::*;
    pub use crate::rhs::*;
    pub use crate::run_corries;
    pub use crate::set_Physics_and_E;
    pub use crate::time::*;
    pub use crate::units::*;
    pub use crate::writer::*;
}

pub use prelude::*;

/// TODO
pub fn run_corries<P: Physics + Collectable, N: NumFlux, T: TimeSolver<P>, const E: usize, const S: usize>(
    u: &mut P,
    rhs: &mut Rhs<P, N, S>,
    time: &mut Time<P, T>,
    mesh: &Mesh<S>,
    writer: &mut Writer,
) -> Result<()> {
    if writer.print_banner {
        print_banner();
    }
    writer
        .write_metadata::<S>()
        .context("Calling writer.write_metadata in run_sim")?;
    loop {
        if time.timestep.t >= time.timestep.t_next_output - time.timestep.dt_min {
            time.timestep.t_next_output += time.timestep.dt_output;
            writer
                .update_data(u, time, mesh)
                .context("Calling writer.update_data_matrices in run_sim")?;
            writer
                .write_output()
                .context("Calling wirter.write_output in run_sim")?;
        }
        if time.timestep.t + time.timestep.dt_min * 0.01 >= time.timestep.t_end {
            break;
        }

        if let err @ Err(_) = time.next_solution::<N, E, S>(u, rhs, mesh) {
            time.timestep.dt_kind = DtKind::ErrorDump;
            writer
                .update_data(u, time, mesh)
                .context("Calling writer.update_data_matrices during the error dump in run_sim")?;
            writer
                .write_output()
                .context("Calling writer.write_output during the error dump in run_sim")?;
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
    println!("# Copyright (c) 2022");
    println!("# Author: tbreslein <github.com/tbreslein>");
    println!("# License: MIT");
    println!("# ****************************************");
}
