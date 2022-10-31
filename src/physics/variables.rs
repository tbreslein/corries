// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use ndarray::{Array1, Array2};

use crate::config::physicsconfig::{PhysicsConfig, PhysicsMode};

pub struct Variables<const S: usize, const EQ: usize> {
    /// The type of physics equations we are solving
    pub mode: PhysicsMode,

    /// Number of equations
    pub n_equations: usize,

    /// Primitive variables
    pub p: Array2<f64>,

    /// Conservative variables
    pub c: Array2<f64>,

    /// Speed of sound
    pub c_sound: Array1<f64>,

    /// Possible equation index for density
    pub jdensity: Option<usize>,

    /// Possible equation index for velocity in the xi direction
    pub jxivelocity: Option<usize>,

    /// Possible equation index for momentum in the xi direction
    pub jximomentum: Option<usize>,

    /// Adiabatic index
    pub adiabatic_index: f64,

    /// Whether this type of Physics / these variables are adiabatic
    pub is_adiabatic: bool,

    /// Whether this type of Physics / these variables are isothermal
    pub is_isothermal: bool,
}

impl<const S: usize, const EQ: usize> Variables<S, EQ> {
    /// Builds a new `Variables` object.
    ///
    /// # Arguments
    ///
    /// * `physicsconf` - `PhysicsConfig` containing configuration to build `Physics` and `Variables` objects
    pub fn new(physicsconf: &PhysicsConfig) -> Variables<S, EQ> {
        let mode = physicsconf.mode;
        let n_equations = match mode {
            PhysicsMode::Euler1DIsot => 2,
        };

        // let mut v = Array1::zeros(S);
        let p = Array2::zeros((EQ, S));
        let c = Array2::zeros((EQ, S));
        let c_sound = Array1::zeros(S);

        let jdensity = match mode {
            PhysicsMode::Euler1DIsot => Some(0),
        };
        let jxivelocity = match mode {
            PhysicsMode::Euler1DIsot => Some(1),
        };
        let jximomentum = match mode {
            PhysicsMode::Euler1DIsot => Some(1),
        };

        let adiabatic_index = physicsconf.adiabatic_index;
        let is_adiabatic = match mode {
            PhysicsMode::Euler1DIsot => false,
        };
        let is_isothermal = !is_adiabatic;

        return Variables {
            mode,
            n_equations,
            p,
            c,
            c_sound,
            jdensity,
            jxivelocity,
            jximomentum,
            adiabatic_index,
            is_adiabatic,
            is_isothermal,
        };
    }

    // NOTE: When I start implementing Physics modes that do not support a value explicitly, I will
    // need match statements over the physics mode here.
    // pub fn get_density(&self, i: usize) -> Option<f64> {
    //     return match self.jdensity {
    //         Some(j) => Some(self.p[[i, j]]),
    //         None => None,
    //     };
    // }
    //
    // pub fn get_xi_velocity(&self, i: usize) -> Option<f64> {
    //     return match self.jxivelocity {
    //         Some(j) => Some(self.p[[i, j]]),
    //         None => None,
    //     };
    // }
    //
    // pub fn get_xi_momentum(&self, i: usize) -> Option<f64> {
    //     return match self.jximomentum {
    //         Some(j) => Some(self.p[[i, j]]),
    //         None => None,
    //     };
    // }
    //
    // pub fn get_density_unsafe(&self, i: usize) -> f64 {
    //     return self.p[[i, self.jdensity.unwrap()]];
    // }
    //
    // pub fn get_xi_velocity_unsafe(&self, i: usize) -> f64 {
    //     return self.p[[i, self.jxivelocity.unwrap()]];
    // }
    //
    // pub fn get_xi_momentum_unsafe(&self, i: usize) -> f64 {
    //     return self.c[[i, self.jximomentum.unwrap()]];
    // }
}
