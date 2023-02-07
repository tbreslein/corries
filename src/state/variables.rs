// Copyright (c) 2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [Variables] struct

use crate::{errorhandling::Validation, Collectable, Data, DataName, PhysicsConfig, StructAssociation};
use color_eyre::{
    eyre::{bail, ensure},
    Result,
};
use ndarray::{Array1, Array2, ArrayView1};

/// Contains the variables of a physics system
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Variables<const E: usize, const S: usize> {
    /// Primitive variables
    pub prim: Array2<f64>,

    /// Conservative variables
    pub cons: Array2<f64>,

    /// Speed of sound
    pub c_sound: Array1<f64>,

    /// Eigen values
    pub eigen_vals: Array2<f64>,

    /// Physical flux
    pub flux: Array2<f64>,

    /// Adiabatic index
    pub gamma: f64,
}

unsafe impl<const E: usize, const S: usize> Send for Variables<E, S> {}
unsafe impl<const E: usize, const S: usize> Sync for Variables<E, S> {}

impl<const E: usize, const S: usize> Variables<E, S> {
    /// Constructs a new [Variables] object
    pub fn new(physics_config: &PhysicsConfig) -> Self {
        Variables {
            prim: Array2::zeros((E, S)),
            cons: Array2::zeros((E, S)),
            c_sound: Array1::zeros(S),
            eigen_vals: Array2::zeros((E, S)),
            flux: Array2::zeros((E, S)),
            gamma: physics_config.adiabatic_index,
        }
    }

    /// Assign the fields of rhs to self.
    pub fn assign(&mut self, rhs: &Self) {
        self.prim.assign(&rhs.prim);
        self.cons.assign(&rhs.cons);
        self.c_sound.assign(&rhs.c_sound);
    }

    /// Return a view to the vector of minimal eigen values.
    pub fn eigen_min(&self) -> ArrayView1<f64> {
        self.eigen_vals.row(0)
    }

    /// Return a view to the vector of maximal eigen values.
    pub fn eigen_max(&self) -> ArrayView1<f64> {
        self.eigen_vals.row(E - 1)
    }
}

impl<const E: usize, const S: usize> Collectable for Variables<E, S> {
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
        Ok(())
    }
}

impl<const E: usize, const S: usize> Validation for Variables<E, S> {
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
        Ok(())
    }
}
