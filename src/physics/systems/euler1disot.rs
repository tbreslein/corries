// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::Result;
use ndarray::Axis;

use crate::{
    config::physicsconfig::PhysicsConfig,
    mesh::Mesh,
    physics::{variables::Variables, Physics},
    units::Units,
};

pub struct Euler1DIsot {
    pub cent: Variables,
    pub units: Units,
}

impl Euler1DIsot {
    pub fn new(physicsconf: &PhysicsConfig, mesh: &Mesh) -> Result<Box<dyn Physics>> {
        let cent = Variables::new(physicsconf, mesh);
        let units = Units::new(physicsconf.units_mode);
        return Ok(Box::new(Euler1DIsot { cent, units }));
    }
}

impl Physics for Euler1DIsot {
    fn update_cons(&self, vars: &mut Variables) {
        vars.c.index_axis_mut(Axis(0), 0).assign(&vars.p.index_axis(Axis(0), 0));
        vars.c
            .index_axis_mut(Axis(0), 1)
            .assign(&(&vars.p.index_axis(Axis(0), 1) * &vars.p.index_axis(Axis(0), 0)));
    }

    fn update_prim(&self, vars: &mut Variables) {
        vars.p.index_axis_mut(Axis(0), 0).assign(&vars.c.index_axis(Axis(0), 0));
        vars.p
            .index_axis_mut(Axis(0), 1)
            .assign(&(&vars.c.index_axis(Axis(0), 1) / &vars.c.index_axis(Axis(0), 0)));
    }
}
