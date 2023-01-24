// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [BoundaryCondition] trait and [init_boundary_condition] function.

use crate::{mesh::Mesh, variables::Variables, BoundaryMode, CorriesConfig};

use self::custom::CustomBoundaryConditions;

mod custom;

/// Enumerates the directions a boundary condition can be applied to
pub enum Direction {
    /// West or inner border of the computational area
    West,

    /// East or outer border of the computational area
    East,
}

/// Identifies an object that can apply boundary condition to a [Variables] object
pub trait BoundaryCondition<const E: usize, const S: usize> {
    /// Applies the boundary condition to the primitive variables.
    fn apply(&mut self, vars: &mut Variables<E, S>, mesh: &Mesh<S>);
}

/// Initialises a [BoundaryCondition] object
pub fn init_boundary_condition<const E: usize, const S: usize>(
    direction: Direction,
    config: &CorriesConfig,
) -> impl BoundaryCondition<E, S> {
    return match direction {
        Direction::West => match &config.boundary_condition_west {
            BoundaryMode::Custom(v) => CustomBoundaryConditions::new(direction, v),
        },
        Direction::East => match &config.boundary_condition_east {
            BoundaryMode::Custom(v) => CustomBoundaryConditions::new(direction, v),
        },
    };
}
