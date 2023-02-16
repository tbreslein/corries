// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [State] struct

use self::variables::Variables;
pub use self::{physics::Physics, systems::*};
use crate::{
    boundaryconditions::BoundaryCondition, directions::Direction, errorhandling::Validation, Mesh, PhysicsConfig,
};
use color_eyre::Result;
use std::marker::PhantomData;

pub mod physics;
pub mod systems;
pub mod variables;

/// Central struct for [corries](crate), holding the simulation state.
///
/// The first generic parameter determines which system of physics equations we are dealing with.
/// This parameter needs to be a type that implements [Physics].
///
/// These types determine things like how primitive and conservative variables are updated, how
/// accessors function, how physical flux is calculated for that system, etc.
#[derive(Debug, Clone, PartialEq)]
pub struct State<P: Physics<E, S>, const E: usize, const S: usize> {
    /// Variables at the centre of each cell in the mesh
    pub cent: Variables<E, S>,

    /// Variables at the west border of each cell in the mesh
    pub west: Variables<E, S>,

    /// Variables at the east border of each cell in the mesh
    pub east: Variables<E, S>,
    methods: PhantomData<P>,
}

unsafe impl<P: Physics<E, S>, const E: usize, const S: usize> Send for State<P, E, S> {}
unsafe impl<P: Physics<E, S>, const E: usize, const S: usize> Sync for State<P, E, S> {}

impl<P: Physics<E, S>, const E: usize, const S: usize> Default for State<P, E, S> {
    fn default() -> Self {
        Self::new(&PhysicsConfig::default())
    }
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

    /// Initialises self.west and self.east to make sure that values like the speed of sound are
    /// sound (pun intended).
    ///
    /// This is needed for numerical flux schemes that use linear reconstruction to calculate
    /// variable values on cell faces like the [Kurganov-Tadmor](crate::rhs::numflux::kt::Kt)
    /// scheme.
    ///
    /// If your simulation uses that scheme, you should call this method on the [State] object
    /// after applying your initial conditions to the mesh central variables.
    pub fn init_west_east(&mut self) {
        self.west.assign(&self.cent);
        self.east.assign(&self.cent);
    }

    /// Returns whether this [State] is handled adiabatically or isothermally
    pub const fn is_adiabatic(&self) -> bool {
        P::IS_ADIABATIC
    }

    const fn get_vars<const D: u8>(&self) -> &Variables<E, S> {
        match D {
            d if d == Direction::West as u8 => &self.west,
            d if d == Direction::Cent as u8 => &self.cent,
            d if d == Direction::East as u8 => &self.east,
            _ => panic!("called get vars with a template parameter that cannot be cast into from a Direction!"),
        }
    }

    // TODO: make this a const fn, once the feature is added to rustlang to enable mut refs in
    // const contexts.
    fn get_vars_mut<const D: u8>(&mut self) -> &mut Variables<E, S> {
        match D {
            d if d == Direction::West as u8 => &mut self.west,
            d if d == Direction::Cent as u8 => &mut self.cent,
            d if d == Direction::East as u8 => &mut self.east,
            _ => panic!("called get vars with a template parameter that cannot be cast into from a Direction!"),
        }
    }

    /// Update the primitive variables at a certain part of each mesh cell
    ///
    /// The template parameter denotes whether to update the cell centres, west cell faces, or east
    /// cell faces. Pass the parameter as the [Direction] enum cast to u8 to use this function.
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    /// set_Physics_and_E!(Euler1DIsot);
    /// const S: usize = 10;
    /// let mut u: State<P, E, S> = State::default();
    ///
    /// // Updates the primitive variables in the cell centres
    /// u.update_prim_d::<{Direction::Cent as u8}>();
    ///
    /// // Updates the primitive variables at west cell faces
    /// u.update_prim_d::<{Direction::West as u8}>();
    ///
    /// // Updates the primitive variables at east cell faces
    /// u.update_prim_d::<{Direction::East as u8}>();
    /// ```
    pub fn update_prim_d<const D: u8>(&mut self) {
        P::update_prim(self.get_vars_mut::<D>());
    }

    /// Update the primitive variables at the cell centres
    ///
    /// This is a specialisation of [State::update_prim_d()], where the direction is set to be
    /// [Direction::Cent]. The reasoning was simply that interacting with cell centre values is the
    /// most common case, so they should have the most straight forward interface.
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    /// set_Physics_and_E!(Euler1DIsot);
    /// const S: usize = 10;
    /// let mut u: State<P, E, S> = State::default();
    ///
    /// // Updates the primitive variables in the cell centres
    /// u.update_prim();
    /// ```
    pub fn update_prim(&mut self) {
        self.update_prim_d::<{ Direction::Cent as u8 }>();
    }

    /// Update the conservative variables at a certain part of each mesh cell
    ///
    /// The template parameter denotes whether to update the cell centres, west cell faces, or east
    /// cell faces. Pass the parameter as the [Direction] enum cast to u8 to use this function.
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    /// set_Physics_and_E!(Euler1DIsot);
    /// const S: usize = 10;
    /// let mut u: State<P, E, S> = State::default();
    ///
    /// // Updates the conservative variables in the cell centres
    /// u.update_cons_d::<{Direction::Cent as u8}>();
    ///
    /// // Updates the conservative variables at west cell faces
    /// u.update_cons_d::<{Direction::West as u8}>();
    ///
    /// // Updates the conservative variables at east cell faces
    /// u.update_cons_d::<{Direction::East as u8}>();
    /// ```
    pub fn update_cons_d<const D: u8>(&mut self) {
        P::update_cons(self.get_vars_mut::<D>());
    }

    /// Update the conservative variables at the cell centres
    ///
    /// This is a specialisation of [State::update_cons_d()], where the direction is set to be
    /// [Direction::Cent]. The reasoning was simply that interacting with cell centre values is the
    /// most common case, so they should have the most straight forward interface.
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    /// set_Physics_and_E!(Euler1DIsot);
    /// const S: usize = 10;
    /// let mut u: State<P, E, S> = State::default();
    ///
    /// // Updates the conservative variables in the cell centres
    /// u.update_cons();
    /// ```
    pub fn update_cons(&mut self) {
        self.update_cons_d::<{ Direction::Cent as u8 }>();
    }

    /// Updates variables that are not part of the primitive or conservative variables, like speed
    /// of sound, at a certain part of each mesh cell
    ///
    /// The template parameter denotes whether to update the cell centres, west cell faces, or east
    /// cell faces. Pass the parameter as the [Direction] enum cast to u8 to use this function.
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    /// set_Physics_and_E!(Euler1DIsot);
    /// const S: usize = 10;
    /// let mut u: State<P, E, S> = State::default();
    ///
    /// // Updates the derived variables in the cell centres
    /// u.update_derived_variables_d::<{Direction::Cent as u8}>();
    ///
    /// // Updates the derived variables at west cell faces
    /// u.update_derived_variables_d::<{Direction::West as u8}>();
    ///
    /// // Updates the derived variables at east cell faces
    /// u.update_derived_variables_d::<{Direction::East as u8}>();
    /// ```
    pub fn update_derived_variables_d<const D: u8>(&mut self) {
        P::update_derived_variables(self.get_vars_mut::<D>());
    }

    /// Updates variables that are not part of the primitive or conservative variables, like speed
    /// of sound, at cell centres
    ///
    /// This is a specialisation of [State::update_derived_variables_d()], where the direction is
    /// set to [Direction::Cent]. The reasoning was simply that interacting with cell centre values
    /// is the most common case, so they should have the most straight forward interface.
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    /// set_Physics_and_E!(Euler1DIsot);
    /// const S: usize = 10;
    /// let mut u: State<P, E, S> = State::default();
    ///
    /// // Updates the derived variables in the cell centres
    /// u.update_derived_variables();
    /// ```
    pub fn update_derived_variables(&mut self) {
        self.update_derived_variables_d::<{ Direction::Cent as u8 }>();
    }

    /// Updates the speed of the sound at certain parts of each mesh cell
    ///
    /// The template parameter denotes whether to update the cell centres, west cell faces, or east
    /// cell faces. Pass the parameter as the [Direction] enum cast to u8 to use this function.
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    /// set_Physics_and_E!(Euler1DIsot);
    /// const S: usize = 10;
    /// let mut u: State<P, E, S> = State::default();
    ///
    /// // Updates the speed of sound in the cell centres
    /// u.update_c_sound_d::<{Direction::Cent as u8}>();
    ///
    /// // Updates the speed of sound at west cell faces
    /// u.update_c_sound_d::<{Direction::West as u8}>();
    ///
    /// // Updates the speed of sound at east cell faces
    /// u.update_c_sound_d::<{Direction::East as u8}>();
    /// ```
    pub fn update_c_sound_d<const D: u8>(&mut self) {
        P::update_c_sound(self.get_vars_mut::<D>());
    }

    /// Updates the speed of the sound at the cell centres
    ///
    /// This is a specialisation of [State::update_c_sound_d()], where the direction is set to be
    /// [Direction::Cent]. The reasoning was simply that interacting with cell centre values is the
    /// most common case, so they should have the most straight forward interface.
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    /// set_Physics_and_E!(Euler1DIsot);
    /// const S: usize = 10;
    /// let mut u: State<P, E, S> = State::default();
    ///
    /// // Updates the speed of sound in the cell centres
    /// u.update_c_sound();
    /// ```
    pub fn update_c_sound(&mut self) {
        self.update_c_sound_d::<{ Direction::Cent as u8 }>();
    }

    /// Updates the eigen values at certain parts of each mesh cell
    ///
    /// The template parameter denotes whether to update the cell centres, west cell faces, or east
    /// cell faces. Pass the parameter as the [Direction] enum cast to u8 to use this function.
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    /// set_Physics_and_E!(Euler1DIsot);
    /// const S: usize = 10;
    /// let mut u: State<P, E, S> = State::default();
    ///
    /// // Updates the eigen values in the cell centres
    /// u.update_eigen_vals_d::<{Direction::Cent as u8}>();
    ///
    /// // Updates the eigen values at west cell faces
    /// u.update_eigen_vals_d::<{Direction::West as u8}>();
    ///
    /// // Updates the eigen values at east cell faces
    /// u.update_eigen_vals_d::<{Direction::East as u8}>();
    /// ```
    pub fn update_eigen_vals_d<const D: u8>(&mut self) {
        P::update_eigen_vals(self.get_vars_mut::<D>());
    }

    /// Updates the eigen values at cell centres
    ///
    /// This is a specialisation of [State::update_eigen_vals_d()], where the direction is set to
    /// be [Direction::Cent]. The reasoning was simply that interacting with cell centre values is
    /// the most common case, so they should have the most straight forward interface.
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    /// set_Physics_and_E!(Euler1DIsot);
    /// const S: usize = 10;
    /// let mut u: State<P, E, S> = State::default();
    ///
    /// // Updates the eigen values in the cell centres
    /// u.update_eigen_vals();
    /// ```
    pub fn update_eigen_vals(&mut self) {
        self.update_eigen_vals_d::<{ Direction::Cent as u8 }>();
    }

    /// Updates the minimal eigen values at certain parts of each mesh cell
    ///
    /// The template parameter denotes whether to update the cell centres, west cell faces, or east
    /// cell faces. Pass the parameter as the [Direction] enum cast to u8 to use this function.
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    /// set_Physics_and_E!(Euler1DIsot);
    /// const S: usize = 10;
    /// let mut u: State<P, E, S> = State::default();
    ///
    /// // Updates the minimal eigen values in the cell centres
    /// u.update_eigen_vals_min_d::<{Direction::Cent as u8}>();
    ///
    /// // Updates the minimal eigen values at west cell faces
    /// u.update_eigen_vals_min_d::<{Direction::West as u8}>();
    ///
    /// // Updates the minimal eigen values at east cell faces
    /// u.update_eigen_vals_min_d::<{Direction::East as u8}>();
    /// ```
    pub fn update_eigen_vals_min_d<const D: u8>(&mut self) {
        P::update_eigen_vals_min(self.get_vars_mut::<D>());
    }

    /// Updates the minimal eigen values at cell centres
    ///
    /// This is a specialisation of [State::update_eigen_vals_min_d()], where the direction is set
    /// to be [Direction::Cent]. The reasoning was simply that interacting with cell centre values
    /// is the most common case, so they should have the most straight forward interface.
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    /// set_Physics_and_E!(Euler1DIsot);
    /// const S: usize = 10;
    /// let mut u: State<P, E, S> = State::default();
    ///
    /// // Updates the minimal eigen values in the cell centres
    /// u.update_eigen_vals_min();
    /// ```
    pub fn update_eigen_vals_min(&mut self) {
        self.update_eigen_vals_min_d::<{ Direction::Cent as u8 }>();
    }

    /// Updates the maximal eigen values at certain parts of each mesh cell
    ///
    /// The template parameter denotes whether to update the cell centres, west cell faces, or east
    /// cell faces. Pass the parameter as the [Direction] enum cast to u8 to use this function.
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    /// set_Physics_and_E!(Euler1DIsot);
    /// const S: usize = 10;
    /// let mut u: State<P, E, S> = State::default();
    ///
    /// // Updates the maximal eigen values in the cell centres
    /// u.update_eigen_vals_max_d::<{Direction::Cent as u8}>();
    ///
    /// // Updates the maximal eigen values at west cell faces
    /// u.update_eigen_vals_max_d::<{Direction::West as u8}>();
    ///
    /// // Updates the maximal eigen values at east cell faces
    /// u.update_eigen_vals_max_d::<{Direction::East as u8}>();
    /// ```
    pub fn update_eigen_vals_max_d<const D: u8>(&mut self) {
        P::update_eigen_vals_max(self.get_vars_mut::<D>());
    }

    /// Updates the maximal eigen values at cell centres
    ///
    /// This is a specialisation of [State::update_eigen_vals_max_d()], where the direction is set
    /// to be [Direction::Cent]. The reasoning was simply that interacting with cell centre values
    /// is the most common case, so they should have the most straight forward interface.
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    /// set_Physics_and_E!(Euler1DIsot);
    /// const S: usize = 10;
    /// let mut u: State<P, E, S> = State::default();
    ///
    /// // Updates the maximal eigen values in the cell centres
    /// u.update_eigen_vals_max();
    /// ```
    pub fn update_eigen_vals_max(&mut self) {
        self.update_eigen_vals_max_d::<{ Direction::Cent as u8 }>();
    }

    /// Updates the physical flux at certain parts of each mesh cell
    ///
    /// The template parameter denotes whether to update the cell centres, west cell faces, or east
    /// cell faces. Pass the parameter as the [Direction] enum cast to u8 to use this function.
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    /// set_Physics_and_E!(Euler1DIsot);
    /// const S: usize = 10;
    /// let mut u: State<P, E, S> = State::default();
    ///
    /// // Updates the physical flux in the cell centres
    /// u.update_flux_d::<{Direction::Cent as u8}>();
    ///
    /// // Updates the physical flux at west cell faces
    /// u.update_flux_d::<{Direction::West as u8}>();
    ///
    /// // Updates the physical flux at east cell faces
    /// u.update_flux_d::<{Direction::East as u8}>();
    /// ```
    pub fn update_flux_d<const D: u8>(&mut self) {
        P::update_flux(self.get_vars_mut::<D>());
    }

    /// Updates the physical flux at cell centres
    ///
    /// This is a specialisation of [State::update_flux_d()], where the direction is set to be
    /// [Direction::Cent]. The reasoning was simply that interacting with cell centre values is the
    /// most common case, so they should have the most straight forward interface.
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    /// set_Physics_and_E!(Euler1DIsot);
    /// const S: usize = 10;
    /// let mut u: State<P, E, S> = State::default();
    ///
    /// // Updates the physical flux in the cell centres
    /// u.update_flux();
    /// ```
    pub fn update_flux(&mut self) {
        self.update_flux_d::<{ Direction::Cent as u8 }>();
    }

    /// Updates all variables at certain parts of each mesh cell, assuming that primitive variables
    /// are up to date
    ///
    /// The template parameter denotes whether to update the cell centres, west cell faces, or east
    /// cell faces. Pass the parameter as the [Direction] enum cast to u8 to use this function.
    pub fn update_vars_from_prim_d<const D: u8>(
        &mut self,
        boundary_west: &mut Box<dyn BoundaryCondition<E, S>>,
        boundary_east: &mut Box<dyn BoundaryCondition<E, S>>,
        mesh: &Mesh<S>,
    ) {
        P::update_vars_from_prim(self.get_vars_mut::<D>(), boundary_west, boundary_east, mesh);
    }

    /// Updates all variables at cell centres, assuming that primitive variables are up to date
    ///
    /// This is a specialisation of [State::update_vars_from_prim_d()], where the direction is set
    /// to be [Direction::Cent]. The reasoning was simply that interacting with cell centre values
    /// is the most common case, so they should have the most straight forward interface.
    pub fn update_vars_from_prim(
        &mut self,
        boundary_west: &mut Box<dyn BoundaryCondition<E, S>>,
        boundary_east: &mut Box<dyn BoundaryCondition<E, S>>,
        mesh: &Mesh<S>,
    ) {
        self.update_vars_from_prim_d::<{ Direction::Cent as u8 }>(boundary_west, boundary_east, mesh);
    }

    /// Updates all variables at certain parts of each mesh cell, assuming that conservative variables
    /// are up to date
    ///
    /// The template parameter denotes whether to update the cell centres, west cell faces, or east
    /// cell faces. Pass the parameter as the [Direction] enum cast to u8 to use this function.
    pub fn update_vars_from_cons_d<const D: u8>(
        &mut self,
        boundary_west: &mut Box<dyn BoundaryCondition<E, S>>,
        boundary_east: &mut Box<dyn BoundaryCondition<E, S>>,
        mesh: &Mesh<S>,
    ) {
        P::update_vars_from_cons(self.get_vars_mut::<D>(), boundary_west, boundary_east, mesh);
    }

    /// Updates all variables at cell centres, assuming that conservative variables are up to date
    ///
    /// This is a specialisation of [State::update_vars_from_cons_d], where the direction is set
    /// to be [Direction::Cent]. The reasoning was simply that interacting with cell centre values
    /// is the most common case, so they should have the most straight forward interface.
    pub fn update_vars_from_cons(
        &mut self,
        boundary_west: &mut Box<dyn BoundaryCondition<E, S>>,
        boundary_east: &mut Box<dyn BoundaryCondition<E, S>>,
        mesh: &Mesh<S>,
    ) {
        self.update_vars_from_cons_d::<{ Direction::Cent as u8 }>(boundary_west, boundary_east, mesh);
    }

    /// Assigns the fields of the argument to `self`.
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    /// set_Physics_and_E!(Euler1DIsot);
    /// const S: usize = 10;
    /// let mut u1: State<P, E, S> = State::default();
    /// let mut u2: State<P, E, S> = State::default();
    ///
    /// u1.cent.prim[[0,0]] = 0.0;
    /// u2.cent.prim[[0,0]] = 1.0;
    ///
    /// u1.west.prim[[0,0]] = 0.0;
    /// u2.west.prim[[0,0]] = 2.0;
    ///
    /// u1.east.prim[[0,0]] = 0.0;
    /// u2.east.prim[[0,0]] = 3.0;
    ///
    /// // Updates all fields of u2 to u1
    /// u1.assign(&u2);
    ///
    /// assert_eq!(u1.cent.prim[[0,0]], 1.0);
    /// assert_eq!(u1.west.prim[[0,0]], 2.0);
    /// assert_eq!(u1.east.prim[[0,0]], 3.0);
    /// ```
    pub fn assign(&mut self, rhs: &Self) {
        self.assign_vars_d::<{ Direction::Cent as u8 }>(rhs);
        self.assign_vars_d::<{ Direction::West as u8 }>(rhs);
        self.assign_vars_d::<{ Direction::East as u8 }>(rhs);
    }

    /// Assigns the argument `rhs` to `self` at certain parts of each mesh cell.
    ///
    /// The template parameter denotes whether to update the cell centres, west cell faces, or east
    /// cell faces. Pass the parameter as the Direction enum cast to u8 to use this function.
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    /// set_Physics_and_E!(Euler1DIsot);
    /// const S: usize = 10;
    /// let physics_config = PhysicsConfig { adiabatic_index: 1.66, units_mode: UnitsMode::SI };
    /// let mut u1: State<P, E, S> = State::default();
    /// let mut u2: State<P, E, S> = State::default();
    ///
    /// u1.cent.prim[[0,0]] = 0.0;
    /// u2.cent.prim[[0,0]] = 1.0;
    ///
    /// u1.west.prim[[0,0]] = 0.0;
    /// u2.west.prim[[0,0]] = 2.0;
    ///
    /// u1.east.prim[[0,0]] = 0.0;
    /// u2.east.prim[[0,0]] = 3.0;
    ///
    /// // Assigns vars from u2 to u1 in the cell centres
    /// u1.assign_vars_d::<{Direction::Cent as u8}>(&u2);
    ///
    /// // Assigns vars from u2 to u1 at west cell faces
    /// u1.assign_vars_d::<{Direction::West as u8}>(&u2);
    ///
    /// // Assigns vars from u2 to u1 at east cell faces
    /// u1.assign_vars_d::<{Direction::East as u8}>(&u2);
    ///
    /// assert_eq!(u1.cent.prim[[0,0]], 1.0);
    /// assert_eq!(u1.west.prim[[0,0]], 2.0);
    /// assert_eq!(u1.east.prim[[0,0]], 3.0);
    /// ```
    pub fn assign_vars_d<const D: u8>(&mut self, rhs: &Self) {
        self.get_vars_mut::<D>().assign(rhs.get_vars::<D>());
    }

    /// Assigns the fields of the argument `rhs` to `self` for cell centric values.
    ///
    /// This is a specialisation of [State::assign_vars_d()], where the direction is set to be
    /// [Direction::Cent]. The reasoning was simply that interacting with cell centre values is the
    /// most common case, so they should have the most straight forward interface.
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    /// set_Physics_and_E!(Euler1DIsot);
    /// const S: usize = 10;
    /// let physics_config = PhysicsConfig { adiabatic_index: 1.66, units_mode: UnitsMode::SI };
    /// let mut u1: State<P, E, S> = State::default();
    /// let mut u2: State<P, E, S> = State::default();
    ///
    /// u1.cent.prim[[0,0]] = 0.0;
    /// u2.cent.prim[[0,0]] = 1.0;
    ///
    /// // Updates all fields of u2 to u1
    /// u1.assign_vars(&u2);
    ///
    /// assert_eq!(u1.cent.prim[[0,0]], 1.0);
    /// ```
    pub fn assign_vars(&mut self, rhs: &Self) {
        self.assign_vars_d::<{ Direction::Cent as u8 }>(rhs);
    }

    /// Calculates the time step width according to the CFL criterium.
    ///
    /// This function may error, for example when the time step width would hit infinity.
    ///
    /// # Arguments
    ///
    /// * `c_cfl`: CFL safety parameter
    /// * `mesh`: The [Mesh] of the simulation
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
                u.update_cons();
                u.update_prim();
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
                u.update_cons();
                u.update_prim();
                assert_relative_eq!(u.cent.prim, u0.cent.prim, max_relative = 1.0e-12);
            }
        }
    }
}
