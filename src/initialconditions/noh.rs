// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

use color_eyre::{eyre::bail, Result};

use crate::{config::physicsconfig::PhysicsMode, physics::Physics};

pub fn init_noh<const S: usize, const EQ: usize>(u: &mut Physics<S, EQ>) -> Result<()> {
    let supported_physics = [
        PhysicsMode::Euler1DAdiabatic,
        PhysicsMode::Euler1DIsot,
        PhysicsMode::Euler2DIsot,
    ];
    if !supported_physics.contains(&u.mode) {
        bail!(
            "Noh initial conditions do not support your chosen PhysicsMode! Got mode = {:?}",
            u.mode
        );
    }

    let breakpoint_index = ((S - 1) as f64 * 0.5) as usize;
    u.prim.fill(0.0);
    u.cons.fill(0.0);
    for i in 0..breakpoint_index {
        u.prim[[u.jdensity, i]] = 1.0;
        u.prim[[u.jxivelocity, i]] = 1.0;
    }
    for i in breakpoint_index..S {
        u.prim[[u.jdensity, i]] = 1.0;
        u.prim[[u.jxivelocity, i]] = -1.0;
    }
    if u.is_adiabatic {
        u.prim.row_mut(u.jpressure).fill(1.0e-5)
    }
    return Ok(());
}
