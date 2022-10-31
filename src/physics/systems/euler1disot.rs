// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use crate::{
    config::physicsconfig::PhysicsConfig,
    physics::{variables::Variables, Physics},
    units::Units,
};

pub struct Euler1DIsot<const S: usize> {
    pub cent: Variables<S, 2>,
    pub units: Units,
}

impl<const S: usize> Euler1DIsot<S> {
    pub fn new(physicsconf: &PhysicsConfig) -> Self {
        let cent = Variables::new(physicsconf);
        let units = Units::new(physicsconf.units_mode);
        return Euler1DIsot { cent, units };
    }
}

impl<const S: usize, const EQ: usize> Physics<S, EQ> for Euler1DIsot<S> {
    fn update_cons(&self, vars: &mut Variables<S, EQ>) {
        vars.c.row_mut(0).assign(&vars.p.row(0));
        vars.c.row_mut(1).assign(&(&vars.p.row(1) * &vars.p.row(0)));
    }

    fn update_prim(&self, vars: &mut Variables<S, EQ>) {
        vars.p.row_mut(0).assign(&vars.c.row(0));
        vars.p.row_mut(1).assign(&(&vars.c.row(1) / &vars.c.row(0)));
    }
}
