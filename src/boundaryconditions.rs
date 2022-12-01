// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [BoundaryConditionContainer] that applies boundary conditions to a `Physics`
//! object.

use crate::mesh::Mesh ;

mod custom;

/// Enumerates the directions a boundary condition can be applied to
pub enum Direction {
    /// West or inner border of the computational area
    West,

    /// East or outer border of the computational area
    East,
}

/// Identifies an object that can apply boundary condition to a `Physics` object
trait BoundaryCondition<P, const S: usize> {
    //<const S: usize, const EQ: usize> {
    /// Applies the condition
    fn apply(&mut self, u: &mut P, mesh: &Mesh<S>);
}
