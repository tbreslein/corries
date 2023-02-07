// Copyright (c) 2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports [NumFluxConfig] for configuring the numerical flux schemes

use crate::errorhandling::Validation;
use color_eyre::{eyre::bail, Result};
use serde::Serialize;

/// Enumerates the different configurations for the different types of numerical flux schemes
///
/// Defaults to Kt with the VanLeer limiter function
#[derive(Debug, Serialize, Copy, Clone, PartialEq)]
pub enum NumFluxConfig {
    /// Hll configuration (i.e., no configuration)
    Hll,
    /// Kt configuration
    Kt {
        /// The type of limiter to be used in reconstruction in the Kt scheme
        limiter_mode: LimiterMode,
    },
}

impl Default for NumFluxConfig {
    fn default() -> Self {
        Self::Kt {
            limiter_mode: LimiterMode::VanLeer,
        }
    }
}

unsafe impl Send for NumFluxConfig {}
unsafe impl Sync for NumFluxConfig {}

impl Validation for NumFluxConfig {
    fn validate(&self) -> Result<()> {
        match self {
            Self::Hll => Ok(()),
            Self::Kt { limiter_mode } => limiter_mode.validate(),
        }
    }
}

/// Enumerates the different kinds of limiter functions used during reconstruction of cell boundary
/// values
///
/// Defaults to VanLeer
#[derive(Debug, Serialize, Copy, Clone, Default, PartialEq)]
pub enum LimiterMode {
    /// No limiter function, just average the differences between the cells
    NoLimiter,
    /// First order MinMod limiter
    MinMod,

    /// Superbee limiter function
    Superbee,

    /// Monocentric limiter (aka MinMod3); needs a parameter passed in that acts as the weight
    /// between differences in neighbouring cells (i.e. one index apart) and the difference between
    /// values in cells arching over one cell (i.e. two indeces apart)
    Monocent(f64),

    /// VanLeer limiter function
    #[default]
    VanLeer,
}

impl Validation for LimiterMode {
    fn validate(&self) -> Result<()> {
        match self {
            Self::NoLimiter | Self::MinMod | Self::Superbee | Self::VanLeer => Ok(()),
            Self::Monocent(p) => {
                if p > &1.0 {
                    Ok(())
                } else {
                    bail!("This must hold: Monocent parameter p > 1.0 ! Got {}", p)
                }
            },
        }
    }
}
