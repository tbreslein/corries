// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use ndarray::Array1;

use crate::config::physicsconfig::{PhysicsConfig, PhysicsMode};

use self::systems::euler1disot::Euler1DIsot;

mod systems;

pub trait Physics<const S: usize, const EQ: usize> {
    fn update_cons(&mut self);
    fn update_prim(&mut self);

    fn assign_prim(&mut self, i: usize, j: usize, x: f64);
    fn assign_cons(&mut self, i: usize, j: usize, x: f64);
    fn assign_c_sound(&mut self, i: usize, x: f64);

    fn assign_prim_vec(&mut self, j: usize, x: &Array1<f64>);
    fn assign_cons_vec(&mut self, j: usize, x: &Array1<f64>);
    fn assign_c_sound_vec(&mut self, x: &Array1<f64>);

    fn get_prim(&self, i: usize, j: usize) -> f64;
    fn get_cons(&self, i: usize, j: usize) -> f64;
    fn get_c_sound(&self, i: usize) -> f64;
}

pub fn init_physics<const S: usize, const EQ: usize>(physicsconf: &PhysicsConfig) -> Box<dyn Physics<S, EQ>> {
    return match physicsconf.mode {
        PhysicsMode::Euler1DIsot => Box::new(Euler1DIsot::new(physicsconf)),
    };
}
