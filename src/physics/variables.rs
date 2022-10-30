// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use ndarray::{ArrayD, IxDyn};

use crate::{config::physicsconfig::{PhysicsConfig, PhysicsMode}, mesh::Mesh};

pub struct Variables {
    /// The type of physics equations we are solving
    pub mode: PhysicsMode,

    /// Number of equations
    pub n_equations: usize,

    /// Primitive variables
    pub p: ArrayD<f64>,

    /// Conservative variables
    pub c: ArrayD<f64>,

    /// Speed of sound
    pub c_sound: ArrayD<f64>,

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

impl Variables {
    /// Builds a new `Variables` object.
    ///
    /// # Arguments
    ///
    /// * `physicsconf` - `PhysicsConfig` containing configuration to build `Physics` and `Variables` objects
    /// * `mesh` - provides information about the mesh in this simulation
    pub fn new(physicsconf: &PhysicsConfig, mesh: &Mesh) -> Variables {
        let mode = physicsconf.mode;
        let n_equations = match mode {
            PhysicsMode::Euler1DIsot => 2,
        };

        let p = ArrayD::from_elem(IxDyn(&[mesh.n_all, n_equations]), 0.0);
        let c = ArrayD::from_elem(IxDyn(&[mesh.n_all, n_equations]), 0.0);
        let c_sound = ArrayD::from_elem(IxDyn(&[mesh.n_all]), 0.0);

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

        return Variables { mode, n_equations, p, c, c_sound, jdensity, jxivelocity, jximomentum, adiabatic_index, is_adiabatic, is_isothermal };
    }
}
