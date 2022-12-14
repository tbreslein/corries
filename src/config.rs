// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [CorriesConfig] structs and its nested structs for configuring Corries simulations.

pub mod meshconfig;
pub mod numericsconfig;
pub mod outputconfig;
pub mod physicsconfig;

use color_eyre::{eyre::Context, Result};
use serde::Serialize;

pub use meshconfig::*;
pub use numericsconfig::*;
pub use outputconfig::*;
pub use physicsconfig::*;

use crate::errorhandling::Validation;

/// Enumerates the different boundary conditions
#[derive(Debug, Serialize)]
pub enum BoundaryMode {
    /// Set of custom boundary conditions applied to each variable
    Custom(Vec<(usize, CustomBoundaryMode)>),
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

/// Struct that carries the full configuration info for a simulation.
///
/// This struct is used in the beginning of a run to initialise all the runtime-objects that are
/// used throughout the simulation.
#[derive(Debug, Serialize)]
pub struct CorriesConfig {
    /// Whether to print the Corries banner to stdout
    pub print_banner: bool,

    /// Config for Mesh objects
    pub mesh_config: MeshConfig,

    /// Config for Physics objects
    pub physics_config: PhysicsConfig,

    /// boundary condition on the west border of the computational area
    pub boundary_condition_west: BoundaryMode,

    /// boundary condition on the east border of the computational area
    pub boundary_condition_east: BoundaryMode,
    /// Config for everything related to numerics
    pub numerics_config: NumericsConfig,

    /// The number of outputs to write during the simulation, not counting output for the initial
    /// state
    pub output_counter_max: usize,

    /// Config for Writer objects
    pub writer_config: Vec<OutputConfig>,
}

impl Validation for CorriesConfig {
    fn validate(&self) -> Result<()> {
        self.mesh_config.validate().context("Validating config.meshconfig")?;
        self.physics_config
            .validate()
            .context("Validating config.physicsconfig")?;
        self.numerics_config
            .validate()
            .context("Validating config.numericsconfig")?;
        for outputconf in self.writer_config.iter() {
            outputconf.validate().context("Validating config.writerconf")?;
        }
        return Ok(());
    }
}
