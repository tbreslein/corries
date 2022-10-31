// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use ndarray::Axis;

use crate::{
    config::physicsconfig::PhysicsConfig,
    physics::{variables::Variables, Physics},
    units::Units,
};

pub struct Euler1DIsot<const S: usize, const EQ: usize> {
    pub cent: Variables<S, EQ>,
    pub units: Units,
}

impl<const S: usize, const EQ: usize> Euler1DIsot<S, EQ> {
    pub fn new(physicsconf: &PhysicsConfig) -> Self {
        let cent = Variables::new(physicsconf);
        let units = Units::new(physicsconf.units_mode);
        return Euler1DIsot { cent, units };
    }
}

impl<const S: usize, const EQ: usize> Physics<S, EQ> for Euler1DIsot<S, EQ> {
    fn update_cons(&self, vars: &mut Variables<S, EQ>) {
        vars.c.index_axis_mut(Axis(0), 0).assign(&vars.p.index_axis(Axis(0), 0));
        vars.c
            .index_axis_mut(Axis(0), 1)
            .assign(&(&vars.p.index_axis(Axis(0), 1) * &vars.p.index_axis(Axis(0), 0)));
    }

    fn update_prim(&self, vars: &mut Variables<S, EQ>) {
        vars.p.index_axis_mut(Axis(0), 0).assign(&vars.c.index_axis(Axis(0), 0));
        vars.p
            .index_axis_mut(Axis(0), 1)
            .assign(&(&vars.c.index_axis(Axis(0), 1) / &vars.c.index_axis(Axis(0), 0)));
    }
}
