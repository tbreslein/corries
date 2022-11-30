// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [CorriesConfig] structs and its nested structs for configuring Corries simulations.

use color_eyre::eyre::Context;
use color_eyre::Result;
use serde::Serialize;

use self::meshconfig::MeshConfig;
// use self::numericsconfig::NumericsConfig;
// use self::outputconfig::OutputConfig;
use self::physicsconfig::PhysicsConfig;
use crate::errorhandling::Validation;

pub mod meshconfig;
// pub mod numericsconfig;
// pub mod outputconfig;
pub mod physicsconfig;

/// Enumerates the different boundary conditions
#[derive(Debug, Serialize)]
pub enum BoundaryMode {
    /// Set of custom boundary conditions applied to each variable
    Custom(Vec<(PhysicsVariable, CustomBoundaryMode)>),
}

/// Enumerates the possible custom boundary conditions
#[derive(Debug, Serialize, Clone, Copy)]
pub enum CustomBoundaryMode {
    /// Extrapolate the values near the boundary into the ghost cells
    Extrapolate,

    /// Specialised version of Extrapolate for mass density in the Kepler case
    ExtrapolateDensityKepler,

    /// Specialised version of Extrapolate for eta velocity in the Kepler case
    ExtrapolateEtaVelocityKepler,

    /// Like NoGradients, but multiplies the value in the ghost cell with a very small number
    NearZero,

    /// No gradients condition; copies the value closest to the boundary into the ghost cells.
    NoGradients,

    /// Reflects outgoing flow, and applies NoGradients to incoming flow
    OutflowNoGradients,

    /// Reflects outgoing flow, and applies Extrapolate to incoming flow
    OutFlowExtrapolate,

    /// Like NoGradients, but switches the sign of the values in the ghost cells
    Reflecting,
}

/// Enumerates the different variables held by `Physics` objects
#[derive(Debug, Serialize, Clone, Copy)]
pub enum PhysicsVariable {
    /// Mass density
    Density,

    /// Velocity along the xi direction
    XiVelocity,

    /// Velocity along the eta direction
    EtaVelocity,

    /// Pressure, duh...
    Pressure,
}

/// Struct that carries the full configuration info for a simulation.
///
/// This struct is used in the beginning of a run to initialise all the runtime-objects that are
/// used throughout the simulation.
#[derive(Debug, Serialize)]
pub struct CorriesConfig {
    /// Whether to print the Corries banner to stdout
    pub print_banner: bool,

    /// Config for Mesh objects
    pub meshconfig: MeshConfig,

    /// Config for Physics objects
    pub physicsconfig: PhysicsConfig,

    /// boundary condition on the west border of the computational area
    pub boundary_condition_west: BoundaryMode,

    /// boundary condition on the east border of the computational area
    pub boundary_condition_east: BoundaryMode,

    // /// Config for everything related to numerics
    // pub numericsconfig: NumericsConfig,
    /// The number of outputs to write during the simulation, not counting output for the initial
    /// state
    pub output_counter_max: usize,
    // /// Config for Writer objects
    // pub writerconfig: Vec<OutputConfig>,
}

// impl Validation for CorriesConfig {
//     fn validate(&self) -> Result<()> {
//         self.meshconfig.validate().context("Validating config.meshconfig")?;
//         self.physicsconfig
//             .validate()
//             .context("Validating config.physicsconfig")?;
//         self.numericsconfig
//             .validate()
//             .context("Validating config.numericsconfig")?;
//         for outputconf in self.writerconfig.iter() {
//             outputconf.validate().context("Validating config.writerconf")?;
//         }
//         return Ok(());
//     }
// }
