// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [Data] struct and associated enums

use super::DataValue;
use ndarray::Array1;
use serde::Serialize;
use std::fmt;

/// Represents a piece of data that can be written into an output stream
#[derive(Debug, Clone, PartialEq)]
pub struct Data {
    /// The identifier of the piece of data
    pub name: DataName,

    /// The payload to be written into the stream
    pub payload: DataValue,

    /// The data type this piece of data is associated with
    pub data_type: DataType,

    /// The struct that the piece of data is extracted from
    pub association: StructAssociation,
}

unsafe impl Send for Data {}
unsafe impl Sync for Data {}

impl Data {
    /// Builds a new `Data` object.
    ///
    /// # Arguments
    ///
    /// * `data` - `DataName` that identifies what this struct will be holding
    /// * `mesh_offset` - The possible offset for vector-like values
    pub fn new<const S: usize>(name: &DataName, mesh_offset: usize) -> Self {
        Self {
            name: *name,
            payload: match name.data_type() {
                DataType::Usize => DataValue::Usize(0),
                DataType::Float => DataValue::Float(0.0),
                DataType::String => DataValue::String("".to_string()),
                DataType::VectorFloat => DataValue::VectorFloat(Array1::zeros([S - 2 * mesh_offset])),
            },
            data_type: name.data_type(),
            association: name.association(),
        }
    }
}

impl fmt::Display for Data {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self.name {
            DataName::XiCent => write!(f, "xi_cent"),
            DataName::XiWest => write!(f, "xi_west"),
            DataName::XiEast => write!(f, "xi_east"),
            DataName::Prim(j) => write!(f, "prim_{}", j),
            DataName::Cons(j) => write!(f, "cons_{}", j),
            DataName::CSound => write!(f, "c_sound"),
            DataName::T => write!(f, "t"),
            DataName::Iter => write!(f, "iter"),
            DataName::Dt => write!(f, "dt"),
            DataName::DtKind => write!(f, "dt_kind"),
        }
    }
}

/// Enumerates the different pieces of data that can be written to output
#[derive(Debug, Serialize, Clone, Copy, PartialEq, Eq)]
pub enum DataName {
    // ====
    // Mesh
    // ====
    /// xi coordinates at the cell centre (Vector)
    XiCent,

    /// xi coordinates at a cell's west border (Vector)
    XiWest,

    /// xi coordinates at a cell's east border (Vector)
    XiEast,
    // =======
    // Physics
    // =======
    /// Primitive variables at index given through the payload
    Prim(usize), // usize -> matrix index

    /// Conservative variables at index given through the payload
    Cons(usize),

    /// Speed of sound
    CSound,

    // ===============
    // TimeIntegration
    // ===============
    /// Time coordinate
    T,

    /// Current iteration count
    Iter,

    /// Last time step width
    Dt,

    /// What the limiting factor of the last time step width was
    DtKind,
}

unsafe impl Send for DataName {}
unsafe impl Sync for DataName {}

/// Enumerates the different types of data that be written to output.
///
/// This has overlap with corries::writer::DataValue, but I needed an enum without a payload to
/// implement DataName::datatype().
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DataType {
    /// Value is a `usize`
    Usize,

    /// Value is an `f64`
    Float,

    /// Value is a vector-like container of `f64`
    VectorFloat,

    /// Value is a `String`
    String,
}

unsafe impl Send for DataType {}
unsafe impl Sync for DataType {}

/// Enumerates the structs / traits a `DataName` may be associated to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StructAssociation {
    /// Value is associated with `Mesh` objects
    Mesh,

    /// Value is associated with `Physics` objects
    Physics,

    /// Valus is associated with `TimeStep` objects
    TimeStep,
}

unsafe impl Send for StructAssociation {}
unsafe impl Sync for StructAssociation {}

impl DataName {
    /// Returns the `DataType` associated with the `DataName`.
    pub fn data_type(&self) -> DataType {
        match &self {
            Self::XiCent => DataType::VectorFloat,
            Self::XiWest => DataType::VectorFloat,
            Self::XiEast => DataType::VectorFloat,
            Self::Prim(_) => DataType::VectorFloat,
            Self::Cons(_) => DataType::VectorFloat,
            Self::CSound => DataType::VectorFloat,
            Self::Iter => DataType::Usize,
            Self::T => DataType::Float,
            Self::Dt => DataType::Float,
            Self::DtKind => DataType::String,
        }
    }

    /// Returns the `StructAssociation` associated with the `DataName`
    pub fn association(&self) -> StructAssociation {
        match &self {
            Self::XiCent => StructAssociation::Mesh,
            Self::XiWest => StructAssociation::Mesh,
            Self::XiEast => StructAssociation::Mesh,
            Self::Prim(_) => StructAssociation::Physics,
            Self::Cons(_) => StructAssociation::Physics,
            Self::CSound => StructAssociation::Physics,
            Self::Iter => StructAssociation::TimeStep,
            Self::T => StructAssociation::TimeStep,
            Self::Dt => StructAssociation::TimeStep,
            Self::DtKind => StructAssociation::TimeStep,
        }
    }
}
