// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::eyre::ensure;

use crate::errorhandling::Validation;

/// Carries information about how an `Output` is built
#[derive(Debug)]
pub struct OutputConfig {
    /// Type of stream the `Output` writes to
    pub stream_mode: StreamMode,

    /// How the `Output` formats its data
    pub formatter_mode: FormatterMode,

    /// Whether the `Output` converts its incoming data to a single `String` or a `Vec<String>`
    pub string_conversion_mode: ToStringConversionMode,

    /// Which folder to write file output to
    pub folder_name: String,

    /// Which folder to write file output to
    pub file_name: String,

    /// Precision for printing floating point numbers
    pub precision: usize,

    /// Whether to include ghost cells in vector data
    pub should_print_ghostcells: bool,

    /// Whether to print metadata into the stream on simulation startup
    pub should_print_metadata: bool,

    /// Identifiers for the data being written by the `Output`
    pub data: Vec<DataName>,
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
        return Ok(());
    }
}

/// Enumerates the different streams an `Output` may write to
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum StreamMode {
    Stdout,
    File,
}

/// Enumerates whether an `Output` writes scalar or vector values
#[derive(Debug, Clone, Copy)]
pub enum ToStringConversionMode {
    Scalar,
    Vector,
}

/// Enumerates how an `Output` formats its output
#[derive(Debug, Clone, Copy)]
pub enum FormatterMode {
    /// Comma seperated output
    CSV,

    /// Tab seperated output
    TSV,
}

#[derive(Debug, Clone)]
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
}
