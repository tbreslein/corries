// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [BoundaryCondition] trait and [init_boundary_condition] function.

use self::custom::CustomBoundaryConditions;
use crate::{mesh::Mesh, variables::Variables, BoundaryMode, CorriesConfig, CustomBoundaryMode, Direction};

mod custom;

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
    match direction {
        Direction::West => match &config.boundary_condition_west {
            BoundaryMode::Custom(v) => CustomBoundaryConditions::new(direction, v),
            BoundaryMode::NoGradients => {
                let v = (0..E)
                    .into_iter()
                    .map(|j| (j, CustomBoundaryMode::NoGradients))
                    .collect::<Vec<(usize, CustomBoundaryMode)>>();
                CustomBoundaryConditions::new(direction, &v)
            },
        },
        Direction::East => match &config.boundary_condition_east {
            BoundaryMode::Custom(v) => CustomBoundaryConditions::new(direction, v),
            BoundaryMode::NoGradients => {
                let v = (0..E)
                    .into_iter()
                    .map(|j| (j, CustomBoundaryMode::NoGradients))
                    .collect::<Vec<(usize, CustomBoundaryMode)>>();
                CustomBoundaryConditions::new(direction, &v)
            },
        },
        Direction::Cent => panic!("Cannot build BoundaryCondition for Direction::Cent!"),
    }
}
