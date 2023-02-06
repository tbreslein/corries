// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [State] struct

pub mod physics;
pub mod systems;
pub mod variables;

use std::marker::PhantomData;

use color_eyre::Result;

use crate::{boundaryconditions::BoundaryCondition, errorhandling::Validation, Mesh, PhysicsConfig};

pub use self::{physics::Physics, systems::*};

use self::variables::Variables;

/// Central struct for corries, holding the simulation state.
///
/// The first generic parameter determines how the state is handled internally, like are primitive
/// and conservative variables updated, how is the physical flux calculated, etc.
pub struct State<P: Physics<E, S>, const E: usize, const S: usize> {
    /// Variables at the centre of each cell in the mesh
    pub cent: Variables<E, S>,

    /// Variables at the west border of each cell in the mesh
    pub west: Variables<E, S>,

    /// Variables at the east border of each cell in the mesh
    pub east: Variables<E, S>,
    methods: PhantomData<P>,
}

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

impl<P: Physics<E, S>, const E: usize, const S: usize> State<P, E, S> {
    /// Construct a new [State] object
    pub fn new(physics_config: &PhysicsConfig) -> Self {
        Self {
            cent: Variables::new(physics_config),
            west: Variables::new(physics_config),
            east: Variables::new(physics_config),
            methods: PhantomData,
        }
    }

    /// Initialises self.west and self.east
    pub fn init_west_east(&mut self) {
        self.west.assign(&self.cent);
        self.east.assign(&self.cent);
    }

    /// Returns whether this state is handled adiabatically or isothermally
    pub const fn is_adiabatic(&self) -> bool {
        P::IS_ADIABATIC
    }

    const fn get_vars<const D: u8>(&self) -> &Variables<E, S> {
        if D == Direction::West as u8 {
            &self.west
        } else if D == Direction::East as u8 {
            &self.east
        } else {
            &self.cent
        }
    }

    // TODO: make this a const fn, once the feature is added to rustlang to enable mut refs in
    // const contexts.
    fn get_vars_mut<const D: u8>(&mut self) -> &mut Variables<E, S> {
        if D == Direction::West as u8 {
            &mut self.west
        } else if D == Direction::East as u8 {
            &mut self.east
        } else {
            &mut self.cent
        }
    }

    /// Update the primitive variables
    ///
    /// The template parameter denotes whether to update the cell centres, west cell faces, or east
    /// cell faces. Pass the parameter as the Direction enum cast to u8 to use this function.
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    /// let physics_config = PhysicsConfig { adiabatic_index = 1.66, units_mode = UnitsMode::SI };
    /// let mut u: State<Euler1DIsot, E, S> = State::new(physics_config);
    ///
    /// // Updates the primitive variables in the cell centres
    /// u.update_prim::<Direction::Cent as u8>();
    /// ```
    pub fn update_prim<const D: u8>(&mut self) {
        P::update_prim(self.get_vars_mut::<D>());
    }

    /// Update the conservative variables
    pub fn update_cons<const D: u8>(&mut self) {
        P::update_cons(self.get_vars_mut::<D>());
    }

    /// Updates variables that are not part of the primitive or conservative variables, like speed
    /// of sound, at the centre of the mesh's cells
    pub fn update_derived_variables<const D: u8>(&mut self) {
        P::update_derived_variables(self.get_vars_mut::<D>());
    }

    /// Updates the speed of the sound at the centre of the mesh's cells
    pub fn update_c_sound<const D: u8>(&mut self) {
        P::update_c_sound(self.get_vars_mut::<D>());
    }

    /// Updates the eigen values at the centre of the mesh's cells
    pub fn update_eigen_vals<const D: u8>(&mut self) {
        P::update_eigen_vals(self.get_vars_mut::<D>());
    }

    /// Updates the minimal eigen values at the centre of the mesh's cells
    pub fn update_eigen_vals_min<const D: u8>(&mut self) {
        P::update_eigen_vals_min(self.get_vars_mut::<D>());
    }

    /// Updates the maximal eigen values at the centre of the mesh's cells
    pub fn update_eigen_vals_max<const D: u8>(&mut self) {
        P::update_eigen_vals_max(self.get_vars_mut::<D>());
    }

    /// Updates the physical flux at the centre of the mesh's cells
    pub fn update_flux<const D: u8>(&mut self) {
        P::update_flux(self.get_vars_mut::<D>());
    }

    /// Update the full .cent field, assuming that .cent.prim is already up-to-date
    pub fn update_cent_from_prim(
        &mut self,
        boundary_west: &mut Box<dyn BoundaryCondition<E, S>>,
        boundary_east: &mut Box<dyn BoundaryCondition<E, S>>,
        mesh: &Mesh<S>,
    ) {
        P::update_vars_from_prim(&mut self.cent, boundary_west, boundary_east, mesh);
    }

    /// Update the full .west field, assuming that .west.prim is already up-to-date
    pub fn update_west_from_prim(
        &mut self,
        boundary_west: &mut Box<dyn BoundaryCondition<E, S>>,
        boundary_east: &mut Box<dyn BoundaryCondition<E, S>>,
        mesh: &Mesh<S>,
    ) {
        P::update_vars_from_prim(&mut self.west, boundary_west, boundary_east, mesh);
    }

    /// Update the full .east field, assuming that .east.prim is already up-to-date
    pub fn update_east_from_prim(
        &mut self,
        boundary_west: &mut Box<dyn BoundaryCondition<E, S>>,
        boundary_east: &mut Box<dyn BoundaryCondition<E, S>>,
        mesh: &Mesh<S>,
    ) {
        P::update_vars_from_prim(&mut self.east, boundary_west, boundary_east, mesh);
    }

    /// Update the full .cent field, assuming that .cent.cons is already up-to-date
    pub fn update_cent_from_cons(
        &mut self,
        boundary_west: &mut Box<dyn BoundaryCondition<E, S>>,
        boundary_east: &mut Box<dyn BoundaryCondition<E, S>>,
        mesh: &Mesh<S>,
    ) {
        P::update_vars_from_cons(&mut self.cent, boundary_west, boundary_east, mesh);
    }

    /// Update the full .west field, assuming that .west.cons is already up-to-date
    pub fn update_west_from_cons(
        &mut self,
        boundary_west: &mut Box<dyn BoundaryCondition<E, S>>,
        boundary_east: &mut Box<dyn BoundaryCondition<E, S>>,
        mesh: &Mesh<S>,
    ) {
        P::update_vars_from_cons(&mut self.west, boundary_west, boundary_east, mesh);
    }

    /// Update the full .east field, assuming that .east.cons is already up-to-date
    pub fn update_east_from_cons(
        &mut self,
        boundary_west: &mut Box<dyn BoundaryCondition<E, S>>,
        boundary_east: &mut Box<dyn BoundaryCondition<E, S>>,
        mesh: &Mesh<S>,
    ) {
        P::update_vars_from_cons(&mut self.east, boundary_west, boundary_east, mesh);
    }

    /// Assigns the fields of rhs to self
    pub fn assign(&mut self, rhs: &Self) {
        self.assign_cent(rhs);
        self.west.assign(&rhs.west);
        self.east.assign(&rhs.east);
    }

    /// Assigns rhs.cent to self.cent
    pub fn assign_cent(&mut self, rhs: &Self) {
        self.cent.assign(&rhs.cent);
    }

    /// Calculates the time step width according to the CFL criterium
    pub fn calc_dt_cfl(&self, c_cfl: f64, mesh: &Mesh<S>) -> Result<f64> {
        P::calc_dt_cfl(&self.cent.eigen_max(), c_cfl, mesh)
    }
}

impl<P: Physics<E, S>, const E: usize, const S: usize> Validation for State<P, E, S> {
    fn validate(&self) -> Result<()> {
        self.cent.validate()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::{PhysicsConfig, UnitsMode};

    use super::*;
    use proptest::prelude::*;
    const S: usize = 2;
    const PHYSICS_CONFIG: PhysicsConfig = PhysicsConfig {
        adiabatic_index: 1.4,
        units_mode: UnitsMode::SI,
    };

    mod euler1dadiabatic {
        use super::*;
        use approx::assert_relative_eq;
        set_Physics_and_E!(Euler1DAdiabatic);
        proptest! {
            #[test]
            fn conversion(p0 in 0.1f64..10.0, p1 in -10.0f64..10.0, p2 in 0.1f64..10.0) {
                // converting to cons and back to prim should be idempotent
                let mut u0 = State::<P, E, S>::new(&PHYSICS_CONFIG);
                u0.cent.prim.row_mut(0).fill(p0);
                u0.cent.prim.row_mut(1).fill(p1);
                u0.cent.prim.row_mut(2).fill(p2);
                let mut u = State::<P,E,S>::new(&PHYSICS_CONFIG);
                u.cent.prim.row_mut(0).fill(p0);
                u.cent.prim.row_mut(1).fill(p1);
                u.cent.prim.row_mut(2).fill(p2);
                u.update_cons_cent();
                u.update_prim_cent();
                assert_relative_eq!(u.cent.prim.row(0), u0.cent.prim.row(0), max_relative = 1.0e-12);
                assert_relative_eq!(u.cent.prim.row(1), u0.cent.prim.row(1), max_relative = 1.0e-12);
                assert_relative_eq!(u.cent.prim.row(2), u0.cent.prim.row(2), max_relative = 1.0e-8);
            }
        }
    }
    mod euler1disot {
        use super::*;
        use approx::assert_relative_eq;
        set_Physics_and_E!(Euler1DIsot);
        proptest! {
            #[test]
            fn conversion(p0 in 0.1f64..100_000.0, p1 in -100_000.0f64..100_000.0) {
                // converting to cons and back to prim should be idempotent
                let mut u0 = State::<P, E, S>::new(&PHYSICS_CONFIG);
                u0.cent.prim.row_mut(0).fill(p0);
                u0.cent.prim.row_mut(1).fill(p1);
                let mut u = State::<P, E, S>::new(&PHYSICS_CONFIG);
                u.cent.prim.row_mut(0).fill(p0);
                u.cent.prim.row_mut(1).fill(p1);
                u.update_cons_cent();
                u.update_prim_cent();
                assert_relative_eq!(u.cent.prim, u0.cent.prim, max_relative = 1.0e-12);
            }
        }
    }
}
