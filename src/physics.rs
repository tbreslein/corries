// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [Physics] struct that handles the variables and physical state of the simulation.

use color_eyre::{
    eyre::{bail, ensure, Context},
    Result,
};
use ndarray::{Array1, Array2, Zip};

use crate::{
    boundaryconditions::BoundaryConditionContainer,
    config::{
        physicsconfig::{PhysicsConfig, PhysicsMode},
        PhysicsVariable,
    },
    errorhandling::Validation,
    mesh::Mesh,
    units::Units,
    writer::{
        data::{Data, DataName, StructAssociation},
        Collectable,
    },
};

mod systems;

/// Struct that governs the variables and state for a system of differential equations
#[derive(Debug)]
pub struct Physics<const S: usize, const EQ: usize> {
    /// The type of physics equations we are solving
    pub mode: PhysicsMode,

    /// Primitive variables
    pub prim: Array2<f64>,

    /// Conservative variables
    pub cons: Array2<f64>,

    /// Speed of sound
    pub c_sound: Array1<f64>,
    c_sound_helper: Array1<f64>,

    /// Eigen values
    pub eigen_vals: Array2<f64>,

    /// Possible equation index for density
    pub jdensity: usize,

    /// Possible equation index for velocity in the xi direction
    pub jxivelocity: usize,

    /// Possible equation index for momentum in the xi direction
    pub jximomentum: usize,

    /// Possible equation index for velocity in the xi direction
    pub jetavelocity: usize,

    /// Possible equation index for momentum in the xi direction
    pub jetamomentum: usize,

    /// Possible equation index for momentum in the xi direction
    pub jenergy: usize,

    /// Possible equation index for momentum in the xi direction
    pub jpressure: usize,

    /// Adiabatic index
    pub adiabatic_index: f64,

    /// Whether this type of Physics / these variables are adiabatic
    pub is_adiabatic: bool,

    /// Whether this type of Physics / these variables are isothermal
    pub is_isothermal: bool,

    /// Helper struct for handling unit systems
    pub units: Units,
}

impl<const S: usize, const EQ: usize> Physics<S, EQ> {
    /// Initialises a [Physics] object with mesh size `S` and `EQ` equations.
    ///
    /// # Arguments
    ///
    /// * `physicsconf` - Contains configuration for the [Physics] object
    pub fn new(physicsconf: &PhysicsConfig) -> Self {
        let mode = physicsconf.mode;

        let prim = Array2::zeros((EQ, S));
        let cons = Array2::zeros((EQ, S));
        let c_sound = Array1::zeros(S);
        let c_sound_helper = Array1::zeros(S);
        let eigen_vals = Array2::zeros((EQ, S));

        let jdensity = match mode {
            PhysicsMode::Euler1DAdiabatic => 0,
            PhysicsMode::Euler1DIsot => 0,
            PhysicsMode::Euler2DIsot => 0,
        };
        let jxivelocity = match mode {
            PhysicsMode::Euler1DAdiabatic => 1,
            PhysicsMode::Euler1DIsot => 1,
            PhysicsMode::Euler2DIsot => 1,
        };
        let jximomentum = match mode {
            PhysicsMode::Euler1DAdiabatic => 1,
            PhysicsMode::Euler1DIsot => 1,
            PhysicsMode::Euler2DIsot => 1,
        };
        let jetavelocity = match mode {
            PhysicsMode::Euler1DAdiabatic => usize::MAX,
            PhysicsMode::Euler1DIsot => usize::MAX,
            PhysicsMode::Euler2DIsot => 2,
        };
        let jetamomentum = match mode {
            PhysicsMode::Euler1DAdiabatic => usize::MAX,
            PhysicsMode::Euler1DIsot => usize::MAX,
            PhysicsMode::Euler2DIsot => 2,
        };
        let jenergy = match mode {
            PhysicsMode::Euler1DAdiabatic => 2,
            PhysicsMode::Euler1DIsot => usize::MAX,
            PhysicsMode::Euler2DIsot => usize::MAX,
        };
        let jpressure = match mode {
            PhysicsMode::Euler1DAdiabatic => 2,
            PhysicsMode::Euler1DIsot => usize::MAX,
            PhysicsMode::Euler2DIsot => usize::MAX,
        };

        let adiabatic_index = physicsconf.adiabatic_index;
        let is_adiabatic = match mode {
            PhysicsMode::Euler1DIsot | PhysicsMode::Euler2DIsot => false,
            PhysicsMode::Euler1DAdiabatic => true,
        };
        let is_isothermal = !is_adiabatic;
        let units = Units::new(physicsconf.units_mode);

        return Self {
            mode,
            prim,
            cons,
            c_sound,
            c_sound_helper,
            eigen_vals,
            jdensity,
            jxivelocity,
            jximomentum,
            jetavelocity,
            jetamomentum,
            jenergy,
            jpressure,
            adiabatic_index,
            is_adiabatic,
            is_isothermal,
            units,
        };
    }

    /// Assigns a [Physics] to `self`.
    ///
    /// This is basically a copy function, but I did not want to implement `Copy` because that
    /// could lead to unwanted implicit copies.
    ///
    /// # Arguments
    ///
    /// * `rhs` - The object whose fields are assigned to `self`
    pub fn assign(&mut self, rhs: &Physics<S, EQ>) {
        self.prim.assign(&rhs.prim);
        self.cons.assign(&rhs.cons);
        self.c_sound.assign(&rhs.c_sound);
    }

    /// Updates the conservative variables in [self.cons] using the primitive variables [self.prim]
    pub fn update_cons(&mut self) {
        match self.mode {
            PhysicsMode::Euler1DIsot => {
                self.update_cons_euler1d_isot();
            },
            PhysicsMode::Euler1DAdiabatic => {
                self.update_cons_euler1d_adiabatic();
            },
            PhysicsMode::Euler2DIsot => {
                self.update_cons_euler2d_isot();
            },
        }
    }

    /// Updates the primitive variables in [self.cons] using the conservative variables [self.prim]
    pub fn update_prim(&mut self) {
        match self.mode {
            PhysicsMode::Euler1DIsot => {
                self.update_prim_euler1d_isot();
            },
            PhysicsMode::Euler1DAdiabatic => {
                self.update_prim_euler1d_adiabatic();
            },
            PhysicsMode::Euler2DIsot => {
                self.update_prim_euler2d_isot();
            },
        }
    }

    /// Calculates the physical flux and saves it in `flux`.
    pub fn calc_physical_flux(&self, flux: &mut Array2<f64>) {
        match self.mode {
            PhysicsMode::Euler1DIsot => {
                self.calc_physical_flux_euler1d_isot(flux);
            },
            PhysicsMode::Euler1DAdiabatic => {
                self.calc_physical_flux_euler1d_adiabatic(flux);
            },
            PhysicsMode::Euler2DIsot => {
                self.calc_physical_flux_euler2d_isot(flux);
            },
        };
    }

    /// Updates values in [Physics] outside of `self.prim` and `self.cons`
    pub fn update_derived_variables(&mut self) {
        self.update_c_sound();
        self.update_eigen_vals();
    }

    /// Updates the speed of sound vector
    fn update_c_sound(&mut self) {
        match self.mode {
            PhysicsMode::Euler1DIsot | PhysicsMode::Euler2DIsot => {
                self.update_c_sound_isot();
            },
            PhysicsMode::Euler1DAdiabatic => {
                self.update_c_sound_adiabatic();
            },
        }
    }

    /// Updates the eigen values
    fn update_eigen_vals(&mut self) {
        self.eigen_vals
            .row_mut(0)
            .assign(&(&self.prim.row(self.jxivelocity) - &self.c_sound));
        self.eigen_vals
            .row_mut(EQ - 1)
            .assign(&(&self.prim.row(self.jxivelocity) + &self.c_sound));
        for j in 1..(EQ - 1) {
            self.eigen_vals
                .row_mut(j)
                .assign(&(&self.prim.row(self.jxivelocity) - &self.c_sound));
        }
    }

    /// Calculates the CFL time step width
    pub fn calc_dt_cfl(&mut self, dt_cfl_param: f64, mesh: &Mesh<S>) -> Result<f64> {
        self.update_eigen_vals();
        let dt = dt_cfl_param
            / Zip::from(self.eigen_vals.row(EQ - 1))
                .and(&mesh.cell_width_inv)
                .fold(0.0f64, |acc, eigenval, cw| acc.max(f64::abs(eigenval * cw)));
        if !dt.is_finite() {
            bail!("dt_cfl turned non-finite! Got dt_cfl = {}", dt);
        }
        return Ok(dt);
    }

    /// Converts a [PhysicsVariable] `var` to the corresponding index for `self.mode`.
    ///
    /// Returns an `Err` if the `var` is not part of this mode.
    pub fn convert_physics_variable_to_index(&self, var: PhysicsVariable) -> Result<usize> {
        return match var {
            PhysicsVariable::Density => Ok(self.jdensity),
            PhysicsVariable::XiVelocity => Ok(self.jxivelocity),
            PhysicsVariable::EtaVelocity => match self.mode {
                PhysicsMode::Euler1DIsot | PhysicsMode::Euler1DAdiabatic => bail!(
                    "var {:?} cannot be converted to usize, because it is not part of the PhysicsMode {:?}!",
                    var,
                    self.mode
                ),
                PhysicsMode::Euler2DIsot => Ok(self.jetavelocity),
            },
            PhysicsVariable::Pressure => match self.mode {
                PhysicsMode::Euler1DIsot | PhysicsMode::Euler2DIsot => bail!(
                    "var {:?} cannot be converted to usize, because it is not part of the PhysicsMode {:?}!",
                    var,
                    self.mode
                ),
                PhysicsMode::Euler1DAdiabatic => Ok(self.jpressure),
            },
        };
    }

    /// Assume that `self.prim` is up to date, update derived variables, boundaries and
    /// conservative variables.
    pub fn update_everything_from_prim(
        &mut self,
        boundary_conditions: &mut BoundaryConditionContainer<S, EQ>,
        mesh: &Mesh<S>,
    ) -> Result<()> {
        self.update_derived_variables();
        boundary_conditions.apply(self, mesh);
        self.update_cons();
        self.validate()
            .context("Calling u.validate in u.update_everything_from_prim")?;
        return Ok(());
    }

    /// Assume that `self.cons` is up to date, update derived variables, boundaries and
    /// primitive variables.
    pub fn update_everything_from_cons(
        &mut self,
        boundary_conditions: &mut BoundaryConditionContainer<S, EQ>,
        mesh: &Mesh<S>,
    ) -> Result<()> {
        self.update_prim();
        self.update_everything_from_prim(boundary_conditions, mesh)
            .context("Calling u.update_everything_from_prim inside u.update_everything_from_cons")?;
        return Ok(());
    }
}

impl<const S: usize, const EQ: usize> Collectable for Physics<S, EQ> {
    fn collect_data(&self, data: &mut Data, mesh_offset: usize) -> Result<()> {
        match (data.association, data.name) {
            (StructAssociation::Physics, DataName::Prim(j)) => self.write_vector(&self.prim.row(j), data, mesh_offset),
            (StructAssociation::Physics, DataName::Cons(j)) => self.write_vector(&self.cons.row(j), data, mesh_offset),
            (StructAssociation::Physics, DataName::CSound) => {
                self.write_vector(&self.c_sound.view(), data, mesh_offset)
            },
            (StructAssociation::Physics, x) => bail!("Tried associating {:?} with Physics!", x),
            (StructAssociation::Mesh, x) | (StructAssociation::TimeStep, x) => {
                bail!("name.association() for {:?} returned {:?}", x, data.association)
            },
        }?;
        return Ok(());
    }
}

impl<const S: usize, const EQ: usize> Validation for Physics<S, EQ> {
    fn validate(&self) -> Result<()> {
        ensure!(
            self.prim.fold(true, |acc, x| acc && x.is_finite()),
            "Physics::prim must be finite! Got: {}",
            self.prim
        );
        ensure!(
            self.cons.fold(true, |acc, x| acc && x.is_finite()),
            "Physics::cons must be finite! Got: {}",
            self.cons
        );
        match self.mode {
            PhysicsMode::Euler1DAdiabatic => {
                ensure!(
                    self.prim.row(self.jdensity).fold(true, |acc, x| acc && x > &0.0),
                    "Mass density must be positive! Got: {}",
                    self.prim.row(self.jdensity)
                );
                ensure!(
                    self.prim.row(self.jpressure).fold(true, |acc, x| acc && x > &0.0),
                    "Pressure must be positive! Got: {}",
                    self.prim.row(self.jpressure)
                );
                ensure!(
                    self.c_sound.fold(true, |acc, x| acc && x.is_finite()),
                    "Physics::c_sound must be finite! Got: {}",
                    self.c_sound
                );
            },
            PhysicsMode::Euler1DIsot => {
                ensure!(
                    self.prim.row(self.jdensity).fold(true, |acc, x| acc && x > &0.0),
                    "Mass density must be positive! Got: {}",
                    self.prim.row(self.jdensity)
                );
            },
            PhysicsMode::Euler2DIsot => {
                ensure!(
                    self.prim.row(self.jdensity).fold(true, |acc, x| acc && x > &0.0),
                    "Mass density must be positive! Got: {}",
                    self.prim.row(self.jdensity)
                );
            },
        }
        return Ok(());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;
    const S: usize = 1004;

    mod euler1dadiabatic {
        use crate::get_n_equations;

        use super::*;
        use approx::assert_relative_eq;
        const PHYSICS_MODE: PhysicsMode = PhysicsMode::Euler1DAdiabatic;
        const EQ: usize = get_n_equations(PHYSICS_MODE);
        proptest! {
            #[test]
            fn conversion(p0 in 0.1f64..100.0, p1 in -100.0f64..100.0, p2 in 0.1f64..100.0, gamma in 0.1f64..0.9) {
                let physicsconf = PhysicsConfig {
                    adiabatic_index: gamma,
                    mode: PHYSICS_MODE,
                    units_mode: crate::units::UnitsMode::SI,
                };
                let mut u0: Physics<S, EQ> = Physics::new(&physicsconf);
                u0.prim.row_mut(0).fill(p0);
                u0.prim.row_mut(1).fill(p1);
                u0.prim.row_mut(2).fill(p2);
                let mut u: Physics<S, EQ> = Physics::new(&physicsconf);
                u.prim.row_mut(0).fill(p0);
                u.prim.row_mut(1).fill(p1);
                u.prim.row_mut(2).fill(p2);

                // converting to cons and back to prim should be idempotent
                u.update_cons();
                u.update_prim();
                assert_relative_eq!(u.prim, u0.prim, max_relative = 1.0e-8);
            }
        }
    }

    mod euler1disot {
        use crate::get_n_equations;

        use super::*;
        use approx::assert_relative_eq;
        const PHYSICS_MODE: PhysicsMode = PhysicsMode::Euler1DIsot;
        const EQ: usize = get_n_equations(PHYSICS_MODE);
        proptest! {
            #[test]
            fn conversion(p0 in 0.1f64..100_000.0, p1 in -100_000.0f64..100_000.0) {
                let physicsconf = PhysicsConfig {
                    adiabatic_index: 0.5,
                    mode: PHYSICS_MODE,
                    units_mode: crate::units::UnitsMode::SI,
                };
                let mut u0: Physics<S, EQ> = Physics::new(&physicsconf);
                u0.prim.row_mut(0).fill(p0);
                u0.prim.row_mut(1).fill(p1);
                let mut u: Physics<S, EQ> = Physics::new(&physicsconf);
                u.prim.row_mut(0).fill(p0);
                u.prim.row_mut(1).fill(p1);

                // converting to cons and back to prim should be idempotent
                u.update_cons();
                u.update_prim();
                assert_relative_eq!(u.prim, u0.prim, max_relative = 1.0e-12);
            }
        }
    }

    mod euler2disot {
        use crate::get_n_equations;

        use super::*;
        use approx::assert_relative_eq;
        const PHYSICS_MODE: PhysicsMode = PhysicsMode::Euler2DIsot;
        const EQ: usize = get_n_equations(PHYSICS_MODE);
        proptest! {
            #[test]
            fn conversion(p0 in 0.1f64..100_000.0, p1 in -100_000.0f64..100_000.0, p2 in -100_000.0f64..100_000.0) {
                let physicsconf = PhysicsConfig {
                    adiabatic_index: 0.5,
                    mode: PHYSICS_MODE,
                    units_mode: crate::units::UnitsMode::SI,
                };
                let mut u0: Physics<S, EQ> = Physics::new(&physicsconf);
                u0.prim.row_mut(0).fill(p0);
                u0.prim.row_mut(1).fill(p1);
                u0.prim.row_mut(2).fill(p2);
                let mut u: Physics<S, EQ> = Physics::new(&physicsconf);
                u.prim.row_mut(0).fill(p0);
                u.prim.row_mut(1).fill(p1);
                u.prim.row_mut(2).fill(p2);

                // converting to cons and back to prim should be idempotent
                u.update_cons();
                u.update_prim();
                assert_relative_eq!(u.prim, u0.prim, max_relative = 1.0e-12);
            }
        }
    }
}
