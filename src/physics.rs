// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use crate::config::physicsconfig::{PhysicsConfig, PhysicsMode};

use self::{systems::euler1disot::Euler1DIsot, variables::Variables};

mod systems;
pub mod variables;

pub trait Physics<const S: usize, const EQ: usize> {
    fn update_cons(&self, vars: &mut Variables<S, EQ>);
    fn update_prim(&self, vars: &mut Variables<S, EQ>);
}

pub fn init_physics<const S: usize, const EQ: usize>(physicsconf: &PhysicsConfig) -> Box<dyn Physics<S, EQ>> {
    return match physicsconf.mode {
        PhysicsMode::Euler1DIsot => Box::new(Euler1DIsot::new(physicsconf)),
    };
}
