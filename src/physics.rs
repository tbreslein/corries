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
    pub jdensity: usize,

    /// Possible equation index for velocity in the xi direction
    pub jxivelocity: usize,

    /// Possible equation index for momentum in the xi direction
    pub jximomentum: usize,

    /// Possible equation index for velocity in the xi direction
    pub jetavelocity: usize,

    /// Possible equation index for momentum in the xi direction
    pub jetamomentum: usize,

    /// Possible equation index for momentum in the xi direction
    pub jangmomentum: usize,

    /// Possible equation index for momentum in the xi direction
    pub jenergy: usize,

    /// Possible equation index for momentum in the xi direction
    pub jpressure: usize,

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

    // /// Updates conservative density and xi momentum according to the isothermal Euler equations
    // #[inline(always)]
    // fn update_cons_euler1d_isot(&mut self) {
    //     self.cons.row_mut(0).assign(&self.prim.row(0));
    //     self.calc_cons_linear_momentum_euler(1);
    // }
    //
    // /// Updates conservative density, xi momentum, and energy according to the adiabatic Euler equations
    // #[inline(always)]
    // fn update_cons_euler1d_adiabatic(&mut self) {
    //     self.update_cons_euler1d_isot();
    //     self.calc_energy_euler1d();
    // }
    //
    // /// Updates conservative density, xi momentum, and eta momentum according to the adiabatic Euler equations
    // #[inline(always)]
    // fn update_cons_euler2d_isot(&mut self) {
    //     self.update_cons_euler1d_isot();
    //     self.calc_cons_linear_momentum_euler(2);
    // }
    //
    // /// Updates primitive density and xi velocity according to the isothermal Euler equations
    // #[inline(always)]
    // fn update_prim_euler1d_isot(&mut self) {
    //     self.prim.row_mut(0).assign(&self.cons.row(0));
    //     self.calc_prim_linear_velocity_euler(1);
    // }
    //
    // /// Updates primitive density, xi velocity, and pressure according to the adiabatic Euler equations
    // #[inline(always)]
    // fn update_prim_euler1d_adiabatic(&mut self) {
    //     self.update_prim_euler1d_isot();
    //     self.calc_pressure_euler1d();
    // }
    //
    // /// Updates primitive density, xi velocity, and eta velocity according to the adiabatic Euler equations
    // #[inline(always)]
    // fn update_prim_euler2d_isot(&mut self) {
    //     self.update_prim_euler1d_isot();
    //     self.calc_prim_linear_velocity_euler(2);
    // }
    //
    // /// Calculates and updates the linear momentum in [self.cons] at row `j`.
    // #[inline(always)]
    // fn calc_cons_linear_momentum_euler(&mut self, j: usize) {
    //     self.cons.row_mut(j).assign(&(&self.prim.row(j) * &self.prim.row(0)));
    // }
    //
    // /// Calculates and updates the linear velocity in [self.prim] at row `j`.
    // #[inline(always)]
    // fn calc_prim_linear_velocity_euler(&mut self, j: usize) {
    //     self.prim.row_mut(j).assign(&(&self.cons.row(j) / &self.cons.row(0)));
    // }
    //
    // /// Calculates and updates the energy in [self.cons].
    // #[inline(always)]
    // fn calc_energy_euler1d(&mut self) {
    //     self.cons.row_mut(2).assign(
    //         &(&self.prim.row(2) * self.adiabatic_index.recip()
    //             + 0.5 * &self.prim.row(0) * &self.prim.row(1) * &self.prim.row(1)),
    //     );
    // }
    //
    // /// Calculates and updates the pressure in [self.prim].
    // #[inline(always)]
    // fn calc_pressure_euler1d(&mut self) {
    //     self.prim.row_mut(2).assign(
    //         &(self.adiabatic_index
    //             * (&self.cons.row(2) - 0.5 / &self.cons.row(0) * &self.cons.row(1) * &self.cons.row(1))),
    //     );
    // }
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
    let jangmomentum = match mode {
        PhysicsMode::Euler1DAdiabatic => usize::MAX,
        PhysicsMode::Euler1DIsot => usize::MAX,
        PhysicsMode::Euler2DIsot => usize::MAX,
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
        jangmomentum,
        jenergy,
        jpressure,
        adiabatic_index,
        is_adiabatic,
        is_isothermal,
        units,
    };
}
