// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

// #![warn(missing_docs)]

//! Corries - CORrosive RIEman Solver
//!
//! Library to run 1D-hydrodynamics simulations solved with Riemann solvers.
//!
//! TODO

use color_eyre::Result;

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
// mod units;
// pub mod writer;

/// TODO
pub mod prelude {
    pub use crate::config::*;
    pub use crate::mesh::*;
    pub use crate::physics::*;
    pub use crate::rhs::*;
    pub use crate::run_corries;
    pub use crate::set_Physics_and_E;
    pub use crate::time::*;
}

pub use prelude::*;

/// TODO
pub fn run_corries<P: Physics, N: NumFlux, T: TimeSolver<P>, const S: usize>(
    _u: &mut P,
    _rhs: &mut Rhs<P, N, S>,
    _time: &mut Time<P, T>,
    _mesh: &Mesh<S>,
) -> Result<()> {
    print_banner();
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
