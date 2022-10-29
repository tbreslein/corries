use color_eyre::Result;

use crate::config::physicsconfig::{PhysicsConfig, PhysicsMode};

use self::euler1disot::Euler1DIsot;

mod euler1disot;
pub mod variables;

pub trait Physics {}

pub fn init_physics(physicsconf: &PhysicsConfig) -> Result<Box<dyn Physics>> {
    return match physicsconf.mode {
        PhysicsMode::Euler1DIsot => Euler1DIsot::new(physicsconf),
    };
}
