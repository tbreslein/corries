// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [BoundaryConditionContainer] that applies boundary conditions to a `Physics`
//! object.

use crate::{
    config::{BoundaryMode, CorriesConfig},
    mesh::Mesh,
    physics::Physics,
};

use self::custom::CustomBoundaryConditions;

mod custom;

/// Enumerates the directions a boundary condition can be applied to
pub enum Direction {
    /// West or inner border of the computational area
    West,

    /// East or outer border of the computational area
    East,
}

/// Exposes the `apply` method for applying boundary condition to a `Physics` object.
pub struct BoundaryConditionContainer<const S: usize, const EQ: usize> {
    /// west-side boundary condition
    west: Box<dyn BoundaryCondition<S, EQ>>,

    /// east-side boundary condition
    east: Box<dyn BoundaryCondition<S, EQ>>,
}

/// Identifies an object that can apply boundary condition to a `Physics` object
trait BoundaryCondition<const S: usize, const EQ: usize> {
    /// Applies the condition
    fn apply(&mut self, u: &mut Physics<S, EQ>, mesh: &Mesh<S>);
}

impl<const S: usize, const EQ: usize> BoundaryConditionContainer<S, EQ> {
    pub fn new(config: &CorriesConfig, u: &Physics<S, EQ>) -> Self {
        return Self {
            west: match &config.boundary_condition_west {
                BoundaryMode::Custom(modes) => Box::new(CustomBoundaryConditions::new(Direction::West, modes, u)),
            },
            east: match &config.boundary_condition_east {
                BoundaryMode::Custom(modes) => Box::new(CustomBoundaryConditions::new(Direction::East, modes, u)),
            },
        };
    }
    pub fn apply(&mut self, u: &mut Physics<S, EQ>, mesh: &Mesh<S>) {
        self.west.as_mut().apply(u, mesh);
        self.east.as_mut().apply(u, mesh);
    }
}
