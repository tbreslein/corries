// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use ndarray::AssignElem;

use crate::{
    config::physicsconfig::PhysicsConfig,
    physics::{variables::Variables, Physics},
    units::Units,
};

pub struct Euler1DIsot<const S: usize> {
    pub vars: Variables<S, 2>,
    pub units: Units,
}

impl<const S: usize> Euler1DIsot<S> {
    pub fn new(physicsconf: &PhysicsConfig) -> Self {
        let vars = Variables::new(physicsconf);
        let units = Units::new(physicsconf.units_mode);
        return Euler1DIsot { vars, units };
    }
}

impl<const S: usize, const EQ: usize> Physics<S, EQ> for Euler1DIsot<S> {
    fn update_cons(&mut self) {
        self.vars.c.row_mut(0).assign(&self.vars.p.row(0));
        self.vars
            .c
            .row_mut(1)
            .assign(&(&self.vars.p.row(1) * &self.vars.p.row(0)));
    }

    fn update_prim(&mut self) {
        self.vars.p.row_mut(0).assign(&self.vars.c.row(0));
        self.vars
            .p
            .row_mut(1)
            .assign(&(&self.vars.c.row(1) / &self.vars.c.row(0)));
    }

    fn assign_prim(&mut self, i: usize, j: usize, x: f64) {
        self.vars.p[[i, j]] = x;
    }

    fn assign_cons(&mut self, i: usize, j: usize, x: f64) {
        self.vars.c[[i, j]] = x;
    }

    fn assign_prim_vec(&mut self, j: usize, x: &ndarray::Array1<f64>) {
        self.vars.p.row_mut(j).assign(x);
    }

    fn assign_cons_vec(&mut self, j: usize, x: &ndarray::Array1<f64>) {
        self.vars.c.row_mut(j).assign(x);
    }

    fn get_prim(&self, i: usize, j: usize) -> f64 {
        return self.vars.p[[i, j]];
    }

    fn get_cons(&self, i: usize, j: usize) -> f64 {
        return self.vars.c[[i, j]];
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ndarray::Array1;
    use proptest::prelude::*;
    const S: usize = 104;
    // TODO: These are unnecessary. I should just generate these in the proptest
    prop_compose! {
        fn arb_density()
                      (x in 0.1f64..100_000.0) -> Array1<f64> {
                          Array1::from_elem(S, x)
        }
    }
    prop_compose! {
        fn arb_vel()
                      (x in -100_000.0..100_000.0) -> Array1<f64> {
                          Array1::from_elem(S, x)
        }
    }

    proptest! {
        #[test]
        fn conversions_work_correctly(p0 in arb_density(), p1 in arb_vel()) {
            let physconf = PhysicsConfig { adiabatic_index: 0.5, c_sound_0: 1.0, mode: crate::config::physicsconfig::PhysicsMode::Euler1DIsot, units_mode: crate::units::UnitsMode::SI};
            let mut u: Box<dyn Physics<S, 2>> = Box::new(Euler1DIsot::new(&physconf));
            u.assign_prim_vec(0, &p0);
            u.assign_prim_vec(1, &p1);

            u.update_cons();
            u.update_prim();
            for i in 0..S {
                // assert_relative_eq!(p0[i], p0[i]);
                assert_relative_eq!(p0[i], u.get_prim(i, 0), epsilon = 10.0f64.powi(-12));
                assert_relative_eq!(p1[i], u.get_prim(i, 1), epsilon = 10.0f64.powi(-12));
            }
            // assert_relative_eq!(mesh.xi_in, mesh.xi_west[mesh.ixi_in], epsilon = 10.0f64.powi(-12));
        }
    }
}
