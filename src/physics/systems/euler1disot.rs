// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::Result;

use crate::{config::physicsconfig::PhysicsConfig, mesh::Mesh, units::Units, physics::{variables::Variables, Physics}};

pub struct Euler1DIsot {
    pub cent: Variables,
    pub units: Units,
}

impl Euler1DIsot {
    pub fn new(physicsconf: &PhysicsConfig, mesh: &Mesh) -> Result<Box<dyn Physics>> {
        let cent = Variables::new(physicsconf, mesh);
        let units = Units::new(physicsconf.units_mode);
        return Ok(Box::new(Euler1DIsot {cent, units}));
    }
}

impl Physics for Euler1DIsot {}
