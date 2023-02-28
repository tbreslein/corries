// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [BoundaryCondition] trait and [init_boundary_condition] function.

use self::custom::CustomBoundaryConditions;
use crate::{mesh::Mesh, variables::Variables, BoundaryMode, CorriesConfig, CustomBoundaryMode, Direction};

mod custom;

/// Identifies an object that can apply boundary condition to a [Variables] object
pub trait BoundaryCondition<const E: usize, const S: usize> {
    /// Applies the [BoundaryCondition] to the primitive variables.
    ///
    /// # Arguments
    ///
    /// * `vars` - The set of [Variables] to apply the condition to
    /// * `mesh` - The [Mesh] for the simulation
    fn apply(&mut self, vars: &mut Variables<E, S>, mesh: &Mesh<S>);
}

/// Initialises a [BoundaryCondition] object
///
/// # Arguments
///
/// * `direction` - Sets whether this set of custom conditions is used for the west or east end
/// of the computational area.
/// * `config` - The full configuration for the simulation, which also holds info abou which
/// boundary conditions to construct
///
/// # Panics
///
/// Panics if `direction` is set to [Direction::Cent].
///
/// # Examples
///
/// ```ignore
/// use corries::prelude::*;
///
/// // set up constants
/// set_Physics_and_E!(Euler1DAdiabatic);
/// const S: usize = 100;
/// type N = Hll<E,S>;
///
/// let t_end = 0.5;
/// let folder_name = "results";
/// let file_name = "noh";
///
/// // define the config instance
/// let config = CorriesConfig::default_riemann_test::<N, E, S>(t_end, folder_name, file_name);
///
/// let boundary_west = init_boundary_condition::<E,S>(Direction::West, &config);
/// let boundary_east = init_boundary_condition::<E,S>(Direction::East, &config);
/// ```
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
