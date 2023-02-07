// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [Directions] enum

use serde::Serialize;

/// Enumerates the different directions inside a cell we can pull values from
///
/// Defaults to Cent
#[repr(u8)]
#[derive(Debug, Serialize, Copy, Clone, Default, PartialEq, Eq)]
pub enum Direction {
    /// west facing cell face
    West,
    /// cell centre (default)
    #[default]
    Cent,
    /// east facing cell face
    East,
}

unsafe impl Send for Direction {}
unsafe impl Sync for Direction {}
