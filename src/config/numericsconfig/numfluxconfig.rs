// Copyright (c) 2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports [NumFluxConfig] for configuring the numerical flux schemes

use color_eyre::{eyre::bail, Result};
use serde::Serialize;

use crate::errorhandling::Validation;

/// todo
#[derive(Debug, Serialize)]
pub enum NumFluxConfig {
    /// todo
    Hll,
    /// todo
    Kt {
        /// todo
        limiter_mode: LimiterMode,
    },
}

impl Validation for NumFluxConfig {
    fn validate(&self) -> Result<()> {
        return match self {
            Self::Hll => Ok(()),
            Self::Kt { limiter_mode } => limiter_mode.validate(),
        };
    }
}

/// todo
#[derive(Debug, Serialize, Copy, Clone)]
pub enum LimiterMode {
    /// todo
    NoLimiter,
    /// todo
    MinMod,
    /// todo
    Superbee,
    /// todo
    Monocent(f64),
    /// todo
    VanLeer,
}

impl Validation for LimiterMode {
    fn validate(&self) -> Result<()> {
        return match self {
            Self::NoLimiter | Self::MinMod | Self::Superbee | Self::VanLeer => Ok(()),
            Self::Monocent(p) => {
                if p > &1.0 {
                    Ok(())
                } else {
                    bail!("This must hold: Monocent parameter p > 1.0 ! Got {}", p)
                }
            },
        };
    }
}
