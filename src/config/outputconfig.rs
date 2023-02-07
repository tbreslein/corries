// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports [OutputConfig] for configuring the objects that write output.

use color_eyre::eyre::ensure;
use serde::Serialize;
use crate::{
    errorhandling::Validation,
    writer::data::{DataName, DataType},
};

/// Carries information about how an `Output` is built
#[derive(Debug, Serialize, Clone, Default)]
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

unsafe impl Send for OutputConfig {}
unsafe impl Sync for OutputConfig {}

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
                    name.data_type() != DataType::VectorFloat,
                    "You cannot write vectors into an Output in Scalar mode! Got data_name = {:?}",
                    name
                );
            }
        }
        Ok(())
    }
}

/// Enumerates the different streams an `Output` may write to
///
/// Defaults to Stdout.
#[derive(Debug, Serialize, Clone, Copy, Default, PartialEq, Eq)]
pub enum StreamMode {
    /// Streams to stdout
    #[default]
    Stdout,

    /// Streams to a file
    File,
}

unsafe impl Send for StreamMode {}
unsafe impl Sync for StreamMode {}

/// Enumerates whether an `Output` writes scalar or vector values
///
/// Defaults to Scalar.
#[derive(Debug, Serialize, Clone, Copy, Default, PartialEq, Eq)]
pub enum ToStringConversionMode {
    /// Writes scalars
    #[default]
    Scalar,

    /// Writes vectors
    Vector,
}

unsafe impl Send for ToStringConversionMode {}
unsafe impl Sync for ToStringConversionMode {}

/// Enumerates how an `Output` formats its output
///
/// Defaults to TSV.
#[derive(Debug, Serialize, Clone, Copy, Default, PartialEq, Eq)]
pub enum FormattingMode {
    /// Comma seperated output
    CSV,

    /// Tab seperated output
    #[default]
    TSV,
}

unsafe impl Send for FormattingMode {}
unsafe impl Sync for FormattingMode {}
