// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

pub mod physics;
pub mod systems;
pub mod variables;

use std::marker::PhantomData;

use color_eyre::Result;

use crate::{boundaryconditions::BoundaryCondition, errorhandling::Validation, Mesh, PhysicsConfig};

pub use self::{physics::Physics, systems::*};

use self::variables::Variables;

pub struct State<P: Physics<E, S>, const E: usize, const S: usize> {
    pub cent: Variables<E, S>,
    pub west: Variables<E, S>,
    pub east: Variables<E, S>,
    methods: PhantomData<P>,
}

impl<P: Physics<E, S>, const E: usize, const S: usize> State<P, E, S> {
    pub fn new(physics_config: &PhysicsConfig) -> Self {
        return Self {
            cent: Variables::new(physics_config),
            west: Variables::new(physics_config),
            east: Variables::new(physics_config),
            methods: PhantomData,
        };
    }

    pub const fn is_adiabatic(&self) -> bool {
        P::IS_ADIABATIC
    }

    pub fn update_prim_cent(&mut self) {
        P::update_prim(&mut self.cent);
    }

    pub fn update_prim_west(&mut self) {
        P::update_prim(&mut self.west);
    }

    pub fn update_prim_east(&mut self) {
        P::update_prim(&mut self.east);
    }

    pub fn update_cons_cent(&mut self) {
        P::update_cons(&mut self.cent);
    }

    pub fn update_cons_west(&mut self) {
        P::update_cons(&mut self.west);
    }

    pub fn update_cons_east(&mut self) {
        P::update_cons(&mut self.east);
    }

    pub fn update_derived_variables_cent(&mut self) {
        P::update_derived_variables(&mut self.cent);
    }

    pub fn update_derived_variables_west(&mut self) {
        P::update_derived_variables(&mut self.west);
    }

    pub fn update_derived_variables_east(&mut self) {
        P::update_derived_variables(&mut self.east);
    }

    pub fn update_c_sound_cent(&mut self) {
        P::update_c_sound(&mut self.cent);
    }

    pub fn update_c_sound_west(&mut self) {
        P::update_c_sound(&mut self.west);
    }

    pub fn update_c_sound_east(&mut self) {
        P::update_c_sound(&mut self.east);
    }

    pub fn update_eigen_vals_cent(&mut self) {
        P::update_eigen_vals(&mut self.cent);
    }

    pub fn update_eigen_vals_west(&mut self) {
        P::update_eigen_vals(&mut self.west);
    }

    pub fn update_eigen_vals_east(&mut self) {
        P::update_eigen_vals(&mut self.east);
    }

    pub fn update_eigen_vals_min_cent(&mut self) {
        P::update_eigen_vals_min(&mut self.cent);
    }

    pub fn update_eigen_vals_min_west(&mut self) {
        P::update_eigen_vals_min(&mut self.west);
    }

    pub fn update_eigen_vals_min_east(&mut self) {
        P::update_eigen_vals_min(&mut self.east);
    }

    pub fn update_eigen_vals_max_cent(&mut self) {
        P::update_eigen_vals_max(&mut self.cent);
    }

    pub fn update_eigen_vals_max_west(&mut self) {
        P::update_eigen_vals_max(&mut self.west);
    }

    pub fn update_eigen_vals_max_east(&mut self) {
        P::update_eigen_vals_max(&mut self.east);
    }

    pub fn update_flux_cent(&mut self) {
        P::update_flux(&mut self.cent);
    }

    pub fn update_flux_west(&mut self) {
        P::update_flux(&mut self.west);
    }

    pub fn update_flux_east(&mut self) {
        P::update_flux(&mut self.east);
    }

    pub fn update_cent_from_prim(
        &mut self,
        boundary_west: &mut Box<dyn BoundaryCondition<E, S>>,
        boundary_east: &mut Box<dyn BoundaryCondition<E, S>>,
        mesh: &Mesh<S>,
    ) {
        P::update_vars_from_prim(&mut self.cent, boundary_west, boundary_east, mesh);
    }

    pub fn update_west_from_prim(
        &mut self,
        boundary_west: &mut Box<dyn BoundaryCondition<E, S>>,
        boundary_east: &mut Box<dyn BoundaryCondition<E, S>>,
        mesh: &Mesh<S>,
    ) {
        P::update_vars_from_prim(&mut self.west, boundary_west, boundary_east, mesh);
    }

    pub fn update_east_from_prim(
        &mut self,
        boundary_west: &mut Box<dyn BoundaryCondition<E, S>>,
        boundary_east: &mut Box<dyn BoundaryCondition<E, S>>,
        mesh: &Mesh<S>,
    ) {
        P::update_vars_from_prim(&mut self.east, boundary_west, boundary_east, mesh);
    }

    pub fn update_cent_from_cons(
        &mut self,
        boundary_west: &mut Box<dyn BoundaryCondition<E, S>>,
        boundary_east: &mut Box<dyn BoundaryCondition<E, S>>,
        mesh: &Mesh<S>,
    ) {
        P::update_vars_from_cons(&mut self.cent, boundary_west, boundary_east, mesh);
    }

    pub fn update_west_from_cons(
        &mut self,
        boundary_west: &mut Box<dyn BoundaryCondition<E, S>>,
        boundary_east: &mut Box<dyn BoundaryCondition<E, S>>,
        mesh: &Mesh<S>,
    ) {
        P::update_vars_from_cons(&mut self.west, boundary_west, boundary_east, mesh);
    }

    pub fn update_east_from_cons(
        &mut self,
        boundary_west: &mut Box<dyn BoundaryCondition<E, S>>,
        boundary_east: &mut Box<dyn BoundaryCondition<E, S>>,
        mesh: &Mesh<S>,
    ) {
        P::update_vars_from_cons(&mut self.east, boundary_west, boundary_east, mesh);
    }

    pub fn assign(&mut self, rhs: &Self) {
        self.assign_cent(rhs);
        self.west.assign(&rhs.west);
        self.east.assign(&rhs.east);
    }

    pub fn assign_cent(&mut self, rhs: &Self) {
        self.cent.assign(&rhs.cent);
    }

    pub fn calc_dt_cfl(&self, c_cfl: f64, mesh: &Mesh<S>) -> Result<f64> {
        return P::calc_dt_cfl(&self.cent.eigen_max(), c_cfl, mesh);
    }
}

impl<P: Physics<E, S>, const E: usize, const S: usize> Validation for State<P, E, S> {
    fn validate(&self) -> Result<()> {
        self.cent.validate()?;
        return Ok(());
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
