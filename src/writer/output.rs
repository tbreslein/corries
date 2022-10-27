// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use crate::{
    config::outputconfig::{FormatterMode, OutputConfig, StreamMode, ToStringConversionMode, DataName},
    mesh::Mesh,
};

use super::DataValue;

pub struct Output<'a> {
    /// Index offset for printing spatially resolved data, useful for skipping ghost cells
    mesh_offset: usize,

    /// Whether this is this `Output`'s first output in this run
    first_output: bool,

    /// Whether this `Output` should write metadata into its `Stream`
    should_print_metadata: bool,

    /// Precision of floating point numbers
    precision: usize,

    /// Matrix of raw, unformatted data, still retaining their original data types
    data_matrix: Vec<&'a DataValue>,

    /// Intermediate string representation of the elements of the `data_matrix`
    string_matrix: Vec<Vec<String>>,

    /// Vector of strings ready to be written to a stream
    stream_strings: Vec<&'a str>,

    /// Writes arbitrary data from `corries` objects to `string_matrix`
    string_conversion_mode: ToStringConversionMode,

    /// Formats the `string_matrix` into `stream_strings` which can be written to a stream
    formatter_mode: FormatterMode,

    /// Handles writing into a stream
    stream_mode: StreamMode,

    /// Identifiers for the data being written to the stream
    data: Vec<DataName>,
}

impl Output<'_> {
    pub fn new(outputconfig: &OutputConfig, mesh: &Mesh) -> Self {
        let rows = match outputconfig.string_conversion_mode {
            ToStringConversionMode::Scalar => 1,
            ToStringConversionMode::Vector => if outputconfig.should_print_ghostcells {mesh.n_all} else {mesh.n_comp},
        };
        let columns = outputconfig.data.len();
        return Output {
            mesh_offset: if outputconfig.should_print_ghostcells {
                mesh.imin
            } else {
                mesh.ixi_in
            },
            first_output: true,
            should_print_metadata: outputconfig.should_print_metadata,
            precision: outputconfig.precision,
            data_matrix: Vec::with_capacity(rows),
            string_matrix: vec![Vec::with_capacity(columns); rows],
            stream_strings: vec![""],
            string_conversion_mode: outputconfig.string_conversion_mode,
            formatter_mode: outputconfig.formatter_mode,
            stream_mode: outputconfig.stream_mode,
            data: outputconfig.data.clone(),
        };
    }
}
