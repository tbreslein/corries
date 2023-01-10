// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [BoundaryCondition] trait and [init_boundary_condition] function.

use crate::{mesh::Mesh, BoundaryMode, CorriesConfig, Physics};

use self::custom::CustomBoundaryConditions;

mod custom;

/// Enumerates the directions a boundary condition can be applied to
pub enum Direction {
    /// West or inner border of the computational area
    West,

    /// East or outer border of the computational area
    East,
}

/// Identifies an object that can apply boundary condition to a `Physics` object
pub trait BoundaryCondition<P, const S: usize> {
    //<const S: usize, const EQ: usize> {
    /// Applies the condition
    fn apply(&mut self, u: &mut P, mesh: &Mesh<S>);
}

pub fn init_boundary_condition<P: Physics, const S: usize>(
    direction: Direction,
    config: &CorriesConfig,
) -> impl BoundaryCondition<P, S> {
    return match direction {
        Direction::West => match &config.boundary_condition_west {
            BoundaryMode::Custom(v) => CustomBoundaryConditions::new(direction, v),
        },
        Direction::East => match &config.boundary_condition_east {
            BoundaryMode::Custom(v) => CustomBoundaryConditions::new(direction, v),
        },
    };
}
