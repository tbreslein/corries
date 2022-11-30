// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports [OutputConfig] for configuring the objects that write output.

use std::fmt;

use color_eyre::eyre::ensure;

use crate::errorhandling::Validation;

/// Carries information about how an `Output` is built
#[derive(Debug)]
pub struct OutputConfig {
    /// Type of stream the `Output` writes to
    pub stream_mode: StreamMode,

    /// How the `Output` formats its data
    pub formatting_mode: FormattingMode,

    /// Whether the `Output` converts its incoming data to a single `String` or a `Vec<String>`
    pub string_conversion_mode: ToStringConversionMode,

    /// Which folder to write file output to
    pub folder_name: String,

    /// Whether the folder contents should be cleared out at the beginning of the run
    pub should_clear_out_folder: bool,

    /// File name for the file(s) to write file output to
    pub file_name: String,

    /// Precision for printing floating point numbers
    pub precision: usize,

    /// Whether to include ghost cells in vector data
    pub should_print_ghostcells: bool,

    /// Whether to print metadata into the stream on simulation startup
    pub should_print_metadata: bool,

    /// Identifiers for the data being written by the `Output`
    pub data_names: Vec<DataName>,
}

impl Validation for OutputConfig {
    fn validate(&self) -> color_eyre::Result<()> {
        if self.stream_mode == StreamMode::File {
            ensure!(
                !self.folder_name.is_empty(),
                "When writing to a file, folder_name may not be empty!"
            );
            ensure!(
                !self.file_name.is_empty(),
                "When writing to a file, file_name may not be empty!"
            );
        }
        if self.string_conversion_mode == ToStringConversionMode::Scalar {
            for name in self.data_names.iter() {
                ensure!(
                    name.datatype() != DataType::VectorFloat,
                    "You cannot write vectors into an Output in Scalar mode! Got data_name = {:?}",
                    name
                );
            }
        }
        return Ok(());
    }
}

/// Enumerates the different streams an `Output` may write to
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StreamMode {
    /// Streams to stdout
    Stdout,

    /// Streams to a file
    File,
}

/// Enumerates whether an `Output` writes scalar or vector values
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ToStringConversionMode {
    /// Writes scalars
    Scalar,

    /// Writes vectors
    Vector,
}

/// Enumerates how an `Output` formats its output
#[derive(Debug, Clone, Copy)]
pub enum FormattingMode {
    /// Comma seperated output
    CSV,

    /// Tab seperated output
    TSV,
}

/// Enumerates the different pieces of data that can be written to output
#[derive(Debug, Clone, Copy)]
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

impl fmt::Display for DataName {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        return match self {
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
        };
    }
}

/// Enumerates the different types of data that be written to output.
///
/// This has overlap with corries::writer::DataValue, but I needed an enum without a payload to
/// implement DataName::datatype().
#[derive(Debug, PartialEq, Eq)]
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

/// Enumerates the structs / traits a `DataName` may be associated to.
#[derive(Debug)]
pub enum StructAssociation {
    /// Value is associated with `Mesh` objects
    Mesh,

    /// Value is associated with `Physics` objects
    Physics,

    /// Valus is associated with `TimeStep` objects
    TimeStep,
}

impl DataName {
    /// Returns the `DataType` associated with the `DataName`.
    pub fn datatype(&self) -> DataType {
        return match &self {
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
        };
    }

    /// Returns the `StructAssociation` associated with the `DataName`
    pub fn association(&self) -> StructAssociation {
        return match &self {
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
        };
    }
}
