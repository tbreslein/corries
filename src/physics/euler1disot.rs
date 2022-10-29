use color_eyre::Result;

use crate::config::physicsconfig::PhysicsConfig;

use super::{variables::Variables, Physics};

pub struct Euler1DIsot {
    pub vars: Variables,
}

impl Euler1DIsot {
    pub fn new(physicsconf: &PhysicsConfig) -> Result<Box<dyn Physics>> {
        return Ok(Box::new(Euler1DIsot {}));
    }
}

impl Physics for Euler1DIsot {}
