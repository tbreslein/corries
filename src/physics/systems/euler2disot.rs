// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use ndarray::{Array1, Array2};

use crate::{
    config::physicsconfig::{PhysicsConfig, PhysicsMode},
    physics::Physics,
    units::Units,
};

pub struct Euler2DIsot<const S: usize> {
    /// The type of physics equations we are solving
    pub mode: PhysicsMode,

    /// Primitive variables
    pub prim: Array2<f64>,

    /// Conservative variables
    pub cons: Array2<f64>,

    /// Speed of sound
    pub c_sound: Array1<f64>,

    /// Possible equation index for density
    pub jdensity: usize,

    /// Possible equation index for velocity in the xi direction
    pub jxivelocity: usize,

    /// Possible equation index for momentum in the xi direction
    pub jximomentum: usize,

    /// Possible equation index for velocity in the xi direction
    pub jetavelocity: usize,

    /// Possible equation index for momentum in the xi direction
    pub jetamomentum: usize,

    /// Adiabatic index
    pub adiabatic_index: f64,

    /// Whether this type of Physics / these variables are adiabatic
    pub is_adiabatic: bool,

    /// Whether this type of Physics / these variables are isothermal
    pub is_isothermal: bool,

    /// Helper struct for handling unit systems
    pub units: Units,
}

impl<const S: usize> Euler2DIsot<S> {
    pub fn new(physicsconf: &PhysicsConfig) -> Self {
        let mode = physicsconf.mode;

        let prim = Array2::zeros((3, S));
        let cons = Array2::zeros((3, S));
        let c_sound = Array1::zeros(S);

        let jdensity = 0;
        let jxivelocity = 1;
        let jximomentum = 1;
        let jetavelocity = 2;
        let jetamomentum = 2;

        let adiabatic_index = physicsconf.adiabatic_index;
        let is_adiabatic = false;
        let is_isothermal = !is_adiabatic;
        let units = Units::new(physicsconf.units_mode);
        return Self {
            mode,
            prim,
            cons,
            c_sound,
            jdensity,
            jxivelocity,
            jximomentum,
            jetavelocity,
            jetamomentum,
            adiabatic_index,
            is_adiabatic,
            is_isothermal,
            units,
        };
    }
}

impl<const S: usize, const EQ: usize> Physics<S, EQ> for Euler2DIsot<S> {
    fn assign_prim(&mut self, j: usize, i: usize, x: f64) {
        self.prim[[j, i]] = x;
    }

    fn assign_cons(&mut self, j: usize, i: usize, x: f64) {
        self.cons[[j, i]] = x;
    }

    fn assign_c_sound(&mut self, i: usize, x: f64) {
        self.c_sound[i] = x;
    }

    fn assign_prim_vec(&mut self, j: usize, x: &ndarray::Array1<f64>) {
        self.prim.row_mut(j).assign(x);
    }

    fn assign_cons_vec(&mut self, j: usize, x: &ndarray::Array1<f64>) {
        self.cons.row_mut(j).assign(x);
    }

    fn assign_c_sound_vec(&mut self, x: &Array1<f64>) {
        self.c_sound.assign(x);
    }

    fn get_prim(&self, j: usize, i: usize) -> f64 {
        return self.prim[[j, i]];
    }

    fn get_cons(&self, j: usize, i: usize) -> f64 {
        return self.cons[[j, i]];
    }

    fn get_c_sound(&self, i: usize) -> f64 {
        return self.c_sound[i];
    }

    fn update_cons(&mut self) {
        self.cons.row_mut(self.jdensity).assign(&self.prim.row(self.jdensity));
        self.cons
            .row_mut(self.jximomentum)
            .assign(&(&self.prim.row(self.jxivelocity) * &self.prim.row(self.jdensity)));
        self.cons
            .row_mut(self.jetamomentum)
            .assign(&(&self.prim.row(self.jetavelocity) * &self.prim.row(self.jdensity)));
    }

    fn update_prim(&mut self) {
        self.prim.row_mut(self.jdensity).assign(&self.cons.row(self.jdensity));
        self.prim
            .row_mut(self.jetavelocity)
            .assign(&(&self.cons.row(self.jetamomentum) / &self.cons.row(self.jdensity)));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use ndarray::Array1;
    use proptest::prelude::*;
    const S: usize = 10;
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
        fn conversions_work_correctly(p0 in arb_density(), p1 in arb_vel(), p2 in arb_vel()) {
            let physconf = PhysicsConfig { adiabatic_index: 0.5, c_sound_0: 1.0, mode: crate::config::physicsconfig::PhysicsMode::Euler1DIsot, units_mode: crate::units::UnitsMode::SI};
            let mut u: Box<dyn Physics<S, 3>> = Box::new(Euler2DIsot::new(&physconf));
            u.assign_prim_vec(0, &p0);
            u.assign_prim_vec(1, &p1);
            u.assign_prim_vec(2, &p2);

            u.update_cons();
            u.update_prim();
            for i in 0..S {
                assert_relative_eq!(p0[i], u.get_prim(0, i), epsilon = 10.0f64.powi(-12));
                assert_relative_eq!(p1[i], u.get_prim(1, i), epsilon = 10.0f64.powi(-12));
                assert_relative_eq!(p2[i], u.get_prim(2, i), epsilon = 10.0f64.powi(-12));
            }
        }
    }
}
