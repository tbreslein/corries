// Copyright (c) 2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports [NumFluxConfig] for configuring [NumFlux](crate::rhs::numflux::NumFlux) objects.

use crate::errorhandling::Validation;
use color_eyre::{eyre::bail, Result};
use serde::Serialize;

/// Enumerates the different configurations for the different types of numerical flux schemes
///
/// Defaults to [Kt](NumFluxConfig::Kt) with the [VanLeer](LimiterMode::VanLeer) limiter function.
#[derive(Debug, Serialize, Copy, Clone, PartialEq)]
pub enum NumFluxConfig {
    /// Configuration for the Harten-Lax-van-Leer solver, i.e. the [Hll](crate::rhs::numflux::Hll)
    /// struct (no further configuration needed)
    Hll,
    /// Configuration for the Kurganov-Tadmor solver, i.e. the [Kt](crate::rhs::numflux::Kt)
    /// struct. Carries one field `limiter_mode` which controls the reconstruction scheme to use.
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
/// values.
///
/// Let
///
/// * `i`: mesh cell index
/// * `j`: equation index
/// * `dxi`: line differential along the `xi` coordinate
/// * `uc`: conservative variables at cell centres
/// * `a`: `uc[[j, i]] - uc[[j, i-1]]`
/// * `b`: `uc[[j, i+1]] - uc[[j, i]]`
/// * `c`: `uc[[j, i+1]] - uc[[j, i-1]]`
/// * `theta`: monocent parameter
///
/// Then the different limiter functions calculate the slop at each mesh cell index `i` and each
/// equation index `j` by:
///
/// * `NoLimiter`: `1 / dxi * 0.5 * (a + b)`
/// * `MinMod`: `if signum(a) * signum(b) > 0 { signum(a) * min(abs(a), abs(b)) } else { 0 }`
/// * `Superbee`:
/// ```text
/// if signum(a) * signum(b) > 0 {
///     signum(a) * min(abs(a), abs(b), 0.5 * max(abs(a), abs(b))
/// } else { 0 }`
/// ```
/// * `Monocent`:
/// ```text
/// if signum(a) * signum(b) && signum(b) * signum(c) > 0 {
///     signum(a) * min(abs(theta * a), abs(theta * b), abs(c))
/// } else { 0 }`
/// ```
/// * `VanLeer`: `(a * abs(b) + b * abs(a)) / (abs(a) + abs(b))`
///
/// Defaults to [VanLeer](LimiterMode::VanLeer)
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
