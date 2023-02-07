// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [Directions] enum

/// Enumerates the different directions inside a cell we can pull values from
#[repr(u8)]
pub enum Direction {
    /// west facing cell face
    West,
    /// cell centre
    Cent,
    /// east facing cell face
    East,
}
