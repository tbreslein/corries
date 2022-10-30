// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::Result;

use crate::{config::physicsconfig::{PhysicsConfig, PhysicsMode}, mesh::Mesh};

use self::systems::euler1disot::Euler1DIsot;

mod systems;
pub mod variables;

pub trait Physics {}

pub fn init_physics(physicsconf: &PhysicsConfig, mesh: &Mesh) -> Result<Box<dyn Physics>> {
    return match physicsconf.mode {
        PhysicsMode::Euler1DIsot => Euler1DIsot::new(physicsconf, mesh),
    };
}
