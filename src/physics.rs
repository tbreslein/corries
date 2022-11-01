// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use ndarray::{Array1, Array2};

use crate::{config::physicsconfig::{PhysicsConfig, PhysicsMode}, units::Units};

use self::systems::{euler1disot::Euler1DIsot, euler2disot::Euler2DIsot};

mod systems;

pub struct PhysicsTest<const S: usize, const EQ: usize> {
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

    /// Adiabatic index
    pub adiabatic_index: f64,

    /// Whether this type of Physics / these variables are adiabatic
    pub is_adiabatic: bool,

    /// Whether this type of Physics / these variables are isothermal
    pub is_isothermal: bool,

    /// Helper struct for handling unit systems
    pub units: Units,
}

pub trait Physics<const S: usize, const EQ: usize> {
    fn update_cons(&mut self);
    fn update_prim(&mut self);

    fn assign_prim(&mut self, j: usize, i: usize, x: f64);
    fn assign_cons(&mut self, j: usize, i: usize, x: f64);
    fn assign_c_sound(&mut self, i: usize, x: f64);

    fn assign_prim_vec(&mut self, j: usize, x: &Array1<f64>);
    fn assign_cons_vec(&mut self, j: usize, x: &Array1<f64>);
    fn assign_c_sound_vec(&mut self, x: &Array1<f64>);

    fn get_prim(&self, j: usize, i: usize) -> f64;
    fn get_cons(&self, j: usize, i: usize) -> f64;
    fn get_c_sound(&self, i: usize) -> f64;
}

pub fn init_physics<const S: usize, const EQ: usize>(physicsconf: &PhysicsConfig) -> Box<dyn Physics<S, EQ>> {
    return match physicsconf.mode {
        PhysicsMode::Euler1DIsot => Box::new(Euler1DIsot::new(physicsconf)),
        PhysicsMode::Euler2DIsot => Box::new(Euler2DIsot::new(physicsconf)),
    };
}

pub fn init_physics_test<const S: usize, const EQ: usize>(physicsconf: &PhysicsConfig) -> PhysicsTest<S, EQ> {
    let mode = physicsconf.mode;

    let prim = Array2::zeros((EQ, S));
    let cons = Array2::zeros((EQ, S));
    let c_sound = Array1::zeros(S);

    let jdensity = 0;
    let jxivelocity = 1;
    let jximomentum = 1;

    let adiabatic_index = physicsconf.adiabatic_index;
    let is_adiabatic = false;
    let is_isothermal = !is_adiabatic;
    let units = Units::new(physicsconf.units_mode);
    return PhysicsTest { mode, prim, cons, c_sound, jdensity, jxivelocity, jximomentum, adiabatic_index, is_adiabatic, is_isothermal, units };
}
