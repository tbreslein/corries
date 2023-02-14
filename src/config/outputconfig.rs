// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports [OutputConfig], which is used to configure different output streams
//! for [Writer](crate::writer::Writer) objects.

use crate::{
    errorhandling::Validation,
    writer::data::{DataName, DataType},
};
use color_eyre::eyre::ensure;
use serde::Serialize;

/// Carries information about how a single `Output` is built
#[derive(Debug, Serialize, Clone)]
pub struct OutputConfig {
    /// Type of stream the `Output` writes to
    pub stream_mode: StreamMode,

    /// How the `Output` formats its data
    pub formatting_mode: FormattingMode,

    /// Whether the `Output` converts its incoming data to a single [String] or a [Vec<String>]
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

impl OutputConfig {
    /// Constructs a [OutputConfig] default configuration that writes to stdout.
    ///
    /// This will set up the output such that the `data_names` field in [OutputConfig] consists of:
    ///
    /// `DataName::Iter, DataName::T, DataName::Dt, DataName::DtKind`
    ///
    /// Check out the asserts in the example to see what the other defaults are.
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    ///
    /// let stdout_config = OutputConfig::default_stdout();
    ///
    /// assert_eq!(stdout_config.stream_mode, StreamMode::Stdout);
    /// assert_eq!(stdout_config.formatting_mode, FormattingMode::TSV);
    /// assert_eq!(stdout_config.string_conversion_mode, ToStringConversionMode::Scalar);
    /// assert_eq!(stdout_config.folder_name, "".to_string());
    /// assert_eq!(stdout_config.should_clear_out_folder, false);
    /// assert_eq!(stdout_config.file_name, "".to_string());
    /// assert_eq!(stdout_config.precision, 3);
    /// assert_eq!(stdout_config.should_print_ghostcells, false);
    /// assert_eq!(stdout_config.should_print_metadata, false);
    /// assert_eq!(stdout_config.data_names[0], DataName::Iter);
    /// assert_eq!(stdout_config.data_names[1], DataName::T);
    /// assert_eq!(stdout_config.data_names[2], DataName::Dt);
    /// assert_eq!(stdout_config.data_names[3], DataName::DtKind);
    /// ```
    pub fn default_stdout() -> Self {
        Self::default_stdout_with_names(vec![DataName::Iter, DataName::T, DataName::Dt, DataName::DtKind])
    }

    /// Constructs an [OutputConfig] default configuration that writes to stdout for a given set of
    /// output values.
    ///
    /// Check out the asserts in the example to see what those defaults are.
    ///
    /// # Arguments
    ///
    /// * `data_names` - a vector of [DataName] objects that will be printed to stdout
    ///
    /// NOTE: `data_names` cannot contain identifiers for vector-like values, like
    /// `DataName::XiCent`, because this configuration has `string_conversion_mode` set to
    /// `ToStringConversionMode::Scalar`!
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    ///
    /// let data_names = vec![DataName::Iter, DataName::T];
    /// let stdout_config = OutputConfig::default_stdout_with_names(data_names);
    ///
    /// assert_eq!(stdout_config.stream_mode, StreamMode::Stdout);
    /// assert_eq!(stdout_config.formatting_mode, FormattingMode::TSV);
    /// assert_eq!(stdout_config.string_conversion_mode, ToStringConversionMode::Scalar);
    /// assert_eq!(stdout_config.folder_name, "".to_string());
    /// assert_eq!(stdout_config.should_clear_out_folder, false);
    /// assert_eq!(stdout_config.file_name, "".to_string());
    /// assert_eq!(stdout_config.precision, 3);
    /// assert_eq!(stdout_config.should_print_ghostcells, false);
    /// assert_eq!(stdout_config.should_print_metadata, false);
    /// assert_eq!(stdout_config.data_names[0], DataName::Iter);
    /// assert_eq!(stdout_config.data_names[1], DataName::T);
    /// ```
    pub fn default_stdout_with_names(data_names: Vec<DataName>) -> Self {
        Self {
            stream_mode: StreamMode::Stdout,
            formatting_mode: FormattingMode::TSV,
            string_conversion_mode: ToStringConversionMode::Scalar,
            folder_name: "".to_string(),
            should_clear_out_folder: false,
            file_name: "".to_string(),
            precision: 3,
            should_print_ghostcells: false,
            should_print_metadata: false,
            data_names,
        }
    }

    /// Constructs a [OutputConfig] default configuration that writes to files for a given set of
    /// output values, folder path and file name.
    ///
    /// The files will be written to the path described by `folder_name`, and will be named
    /// `<file_name>_XXX`, where `XXX` is the output count for that specific write. For example,
    /// let's say `file_name` is set to `noh`, then the files will be named `noh_00.csv`,
    /// `noh_01.csv`, `noh_02.csv`, and so on.
    ///
    /// Check out the asserts in the example to see what the defaults are.
    ///
    /// # Arguments
    ///
    /// * `folder_name` - A string slice describing the relative path to write the results to
    /// * `file_name` - A string slice giving the base file name for the output files
    /// * `num_equations` - Tells the function how many equations you are working with in order to
    /// print all primitive variables
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    ///
    /// let folder_name = "results";
    /// let file_name = "noh";
    /// let stdout_config = OutputConfig::default_file(folder_name, file_name, 3);
    ///
    /// assert_eq!(stdout_config.stream_mode, StreamMode::File);
    /// assert_eq!(stdout_config.formatting_mode, FormattingMode::CSV);
    /// assert_eq!(stdout_config.string_conversion_mode, ToStringConversionMode::Vector);
    /// assert_eq!(stdout_config.folder_name, folder_name.to_string());
    /// assert_eq!(stdout_config.should_clear_out_folder, true);
    /// assert_eq!(stdout_config.file_name, file_name.to_string());
    /// assert_eq!(stdout_config.precision, 7);
    /// assert_eq!(stdout_config.should_print_ghostcells, true);
    /// assert_eq!(stdout_config.should_print_metadata, true);
    /// assert_eq!(stdout_config.data_names[0], DataName::T);
    /// assert_eq!(stdout_config.data_names[1], DataName::XiCent);
    /// assert_eq!(stdout_config.data_names[2], DataName::Prim(0));
    /// assert_eq!(stdout_config.data_names[3], DataName::Prim(1));
    /// assert_eq!(stdout_config.data_names[4], DataName::Prim(2));
    /// ```
    pub fn default_file(folder_name: &str, file_name: &str, num_equations: usize) -> Self {
        let mut data_names = vec![DataName::T, DataName::XiCent];
        for j in 0..num_equations {
            data_names.push(DataName::Prim(j));
        }
        Self::default_file_with_names(folder_name, file_name, data_names)
    }

    /// Constructs a [OutputConfig] default configuration that writes to files for a given set of
    /// output values, folder path and file name.
    ///
    /// The files will be written to the path described by `folder_name`, and will be named
    /// `<file_name>_XXX`, where `XXX` is the output count for that specific write. For example,
    /// let's say `file_name` is set to `noh`, then the files will be named `noh_00.csv`,
    /// `noh_01.csv`, `noh_02.csv`, and so on.
    ///
    /// Check out the asserts in the example to see what the defaults are.
    ///
    /// # Arguments
    ///
    /// * `folder_name` - A string slice describing the relative path to write the results to
    /// * `file_name` - A string slice giving the base file name for the output files
    /// * `data_names` - a vector of `DataName` objects that will be printed to stdout
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    ///
    /// let data_names = vec![DataName::Iter, DataName::T, DataName::Prim(0), DataName::Prim(1)];
    /// let folder_name = "results";
    /// let file_name = "noh";
    /// let stdout_config = OutputConfig::default_file_with_names(folder_name, file_name, data_names);
    ///
    /// assert_eq!(stdout_config.stream_mode, StreamMode::File);
    /// assert_eq!(stdout_config.formatting_mode, FormattingMode::CSV);
    /// assert_eq!(stdout_config.string_conversion_mode, ToStringConversionMode::Vector);
    /// assert_eq!(stdout_config.folder_name, folder_name.to_string());
    /// assert_eq!(stdout_config.should_clear_out_folder, true);
    /// assert_eq!(stdout_config.file_name, file_name.to_string());
    /// assert_eq!(stdout_config.precision, 7);
    /// assert_eq!(stdout_config.should_print_ghostcells, true);
    /// assert_eq!(stdout_config.should_print_metadata, true);
    /// assert_eq!(stdout_config.data_names[0], DataName::Iter);
    /// assert_eq!(stdout_config.data_names[1], DataName::T);
    /// assert_eq!(stdout_config.data_names[2], DataName::Prim(0));
    /// assert_eq!(stdout_config.data_names[3], DataName::Prim(1));
    /// ```
    pub fn default_file_with_names(folder_name: &str, file_name: &str, data_names: Vec<DataName>) -> Self {
        Self {
            stream_mode: StreamMode::File,
            formatting_mode: FormattingMode::CSV,
            string_conversion_mode: ToStringConversionMode::Vector,
            folder_name: folder_name.to_string(),
            should_clear_out_folder: true,
            file_name: file_name.to_string(),
            precision: 7,
            should_print_ghostcells: true,
            should_print_metadata: true,
            data_names,
        }
    }
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

/// Enumerates the different streams an `Output` may write to.
///
/// Defaults to [Stdout](StreamMode::Stdout).
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
/// Defaults to [Scalar](ToStringConversionMode::Scalar).
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
/// Defaults to [TSV](FormattingMode::TSV).
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
