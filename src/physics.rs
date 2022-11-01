// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [Physics] struct that handles the variables and physical state of the simulation.

use ndarray::{Array1, Array2};

use crate::{
    config::physicsconfig::{PhysicsConfig, PhysicsMode},
    units::Units,
};

mod systems;

/// Struct that governs the variables and state for a system of differential equations
#[derive(Debug)]
pub struct Physics<const S: usize, const EQ: usize> {
    /// The type of physics equations we are solving
    pub mode: PhysicsMode,

    /// Primitive variables
    pub prim: Array2<f64>,

    /// Conservative variables
    pub cons: Array2<f64>,

    /// Speed of sound
    pub c_sound: Array1<f64>,

    /// Possible equation index for density
    jdensity: usize,

    /// Possible equation index for velocity in the xi direction
    jxivelocity: usize,

    /// Possible equation index for momentum in the xi direction
    jximomentum: usize,

    /// Possible equation index for velocity in the xi direction
    jetavelocity: usize,

    /// Possible equation index for momentum in the xi direction
    jetamomentum: usize,

    /// Possible equation index for momentum in the xi direction
    jenergy: usize,

    /// Possible equation index for momentum in the xi direction
    jpressure: usize,

    /// Adiabatic index
    pub adiabatic_index: f64,

    /// Whether this type of Physics / these variables are adiabatic
    pub is_adiabatic: bool,

    /// Whether this type of Physics / these variables are isothermal
    pub is_isothermal: bool,

    /// Helper struct for handling unit systems
    pub units: Units,
}

impl<const S: usize, const EQ: usize> Physics<S, EQ> {
    /// Updates the conservative variables in [self.cons] using the primitive variables [self.prim]
    pub fn update_cons(&mut self) {
        match self.mode {
            PhysicsMode::Euler1DIsot => {
                self.update_cons_euler1d_isot();
            },
            PhysicsMode::Euler1DAdiabatic => {
                self.update_cons_euler1d_adiabatic();
            },
            PhysicsMode::Euler2DIsot => {
                self.update_cons_euler2d_isot();
            },
        }
    }

    /// Updates the primitive variables in [self.cons] using the conservative variables [self.prim]
    pub fn update_prim(&mut self) {
        match self.mode {
            PhysicsMode::Euler1DIsot => {
                self.update_prim_euler1d_isot();
            },
            PhysicsMode::Euler1DAdiabatic => {
                self.update_prim_euler1d_adiabatic();
            },
            PhysicsMode::Euler2DIsot => {
                self.update_prim_euler2d_isot();
            },
        }
    }
}

/// Initialises a [Physics] object with mesh size `S` and `EQ` equations.
///
/// # Arguments
///
/// * `physicsconf` - Contains configuration for the [Physics] object
pub fn init_physics<const S: usize, const EQ: usize>(physicsconf: &PhysicsConfig) -> Physics<S, EQ> {
    let mode = physicsconf.mode;

    let prim = Array2::from_elem((EQ, S), 1.1);
    let cons = Array2::from_elem((EQ, S), 1.1);
    let c_sound = Array1::from_elem(S, 1.1);

    let jdensity = match mode {
        PhysicsMode::Euler1DAdiabatic => 0,
        PhysicsMode::Euler1DIsot => 0,
        PhysicsMode::Euler2DIsot => 0,
    };
    let jxivelocity = match mode {
        PhysicsMode::Euler1DAdiabatic => 1,
        PhysicsMode::Euler1DIsot => 1,
        PhysicsMode::Euler2DIsot => 1,
    };
    let jximomentum = match mode {
        PhysicsMode::Euler1DAdiabatic => 1,
        PhysicsMode::Euler1DIsot => 1,
        PhysicsMode::Euler2DIsot => 1,
    };
    let jetavelocity = match mode {
        PhysicsMode::Euler1DAdiabatic => usize::MAX,
        PhysicsMode::Euler1DIsot => usize::MAX,
        PhysicsMode::Euler2DIsot => 2,
    };
    let jetamomentum = match mode {
        PhysicsMode::Euler1DAdiabatic => usize::MAX,
        PhysicsMode::Euler1DIsot => usize::MAX,
        PhysicsMode::Euler2DIsot => 2,
    };
    let jenergy = match mode {
        PhysicsMode::Euler1DAdiabatic => 2,
        PhysicsMode::Euler1DIsot => usize::MAX,
        PhysicsMode::Euler2DIsot => usize::MAX,
    };
    let jpressure = match mode {
        PhysicsMode::Euler1DAdiabatic => 2,
        PhysicsMode::Euler1DIsot => usize::MAX,
        PhysicsMode::Euler2DIsot => usize::MAX,
    };

    let adiabatic_index = physicsconf.adiabatic_index;
    let is_adiabatic = match mode {
        PhysicsMode::Euler1DIsot | PhysicsMode::Euler2DIsot => false,
        PhysicsMode::Euler1DAdiabatic => true,
    };
    let is_isothermal = !is_adiabatic;
    let units = Units::new(physicsconf.units_mode);
    return Physics {
        mode,
        prim,
        cons,
        c_sound,
        jdensity,
        jxivelocity,
        jximomentum,
        jetavelocity,
        jetamomentum,
        jenergy,
        jpressure,
        adiabatic_index,
        is_adiabatic,
        is_isothermal,
        units,
    };
}
