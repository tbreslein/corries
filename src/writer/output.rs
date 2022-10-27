// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use std::path::Path;

use color_eyre::{eyre::bail, Result};
use ndarray::{ArrayD, IxDyn};

use crate::{
    config::{
        outputconfig::{
            DataName, DataType, FormatterMode, OutputConfig, StreamMode, StructAssociation,
            ToStringConversionMode,
        },
        CorriesConfig,
    },
    mesh::Mesh,
};

use super::{DataValue, Write};

pub struct Output {
    /// Index offset for printing spatially resolved data, useful for skipping ghost cells
    mesh_offset: usize,

    /// Whether this is this `Output`'s first output in this run
    first_output: bool,

    /// Current output counter
    output_counter: usize,

    /// Whether this `Output` should write metadata into its `Stream`
    pub should_print_metadata: bool,

    /// Precision of floating point numbers
    precision: usize,

    /// Precision of floating point numbers
    width: usize,

    /// Matrix of raw, unformatted data, still retaining their original data types
    pub data_matrix: Vec<DataValue>,

    /// Intermediate string representation of the elements of the `data_matrix`
    delimiter: char,

    /// Vector of strings ready to be written to a stream
    stream_strings: Vec<String>,

    /// Writes arbitrary data from `corries` objects to `string_matrix`
    string_conversion_mode: ToStringConversionMode,

    /// Handles writing into a stream
    stream_mode: StreamMode,

    /// File name for the file(s) to write file output to
    file_name: String,

    /// How many digits the file counter for file output needs
    file_counter_width: usize,

    /// Identifiers for the data being written to the stream
    data_names: Vec<DataName>,
}

impl Output {
    pub fn new(outputconfig: &OutputConfig, mesh: &Mesh, output_count_max: usize) -> Self {
        let rows = match outputconfig.string_conversion_mode {
            ToStringConversionMode::Scalar => 1,
            ToStringConversionMode::Vector => {
                if outputconfig.should_print_ghostcells {
                    mesh.n_all
                } else {
                    mesh.n_comp
                }
            }
        };
        let mut data_matrix = vec![];
        for name in outputconfig.data_names.iter() {
            match name.datatype() {
                DataType::Int => data_matrix.push(DataValue::Int(0)),
                DataType::Usize => data_matrix.push(DataValue::Usize(0)),
                DataType::Float => data_matrix.push(DataValue::Float(0.0)),
                DataType::String => data_matrix.push(DataValue::String("".to_string())),
                DataType::VectorFloat => data_matrix.push(DataValue::VectorFloat(
                    ArrayD::from_elem(IxDyn(&[rows]), 0.0),
                )),
            };
        }
        return Output {
            mesh_offset: if outputconfig.should_print_ghostcells {
                mesh.imin
            } else {
                mesh.ixi_in
            },
            first_output: true,
            output_counter: 0,
            should_print_metadata: outputconfig.should_print_metadata,
            precision: outputconfig.precision,
            width: match outputconfig.formatter_mode {
                FormatterMode::TSV => 12usize.min(outputconfig.precision + 7),
                FormatterMode::CSV => 0,
            },
            data_matrix,
            delimiter: match outputconfig.formatter_mode {
                FormatterMode::TSV => ' ',
                FormatterMode::CSV => ',',
            },
            stream_strings: vec!["".to_string(); rows],
            string_conversion_mode: outputconfig.string_conversion_mode,
            stream_mode: outputconfig.stream_mode,
            file_name: if outputconfig.folder_name.ends_with("/") {
                outputconfig.folder_name.clone() + &outputconfig.file_name.clone() + "_"
            } else {
                outputconfig.folder_name.clone() + "/" + &outputconfig.file_name.clone() + "_"
            },
            file_counter_width: f64::log10(output_count_max as f64) as usize + 1,
            data_names: outputconfig.data_names.clone(),
        };
    }

    pub fn update_data_matrix(&mut self, mesh: &Mesh) -> Result<()> {
        for (i, name) in self.data_names.iter().enumerate() {
            match name.association() {
                StructAssociation::Mesh => {
                    mesh.collect_data(name, &mut self.data_matrix[i], self.mesh_offset)
                }
            }?;
        }
        return Ok(());
    }

    pub fn write_output(&mut self) -> Result<()> {
        self.data_matrix_to_stream_strings()?;
        self.stream()?;
        self.first_output = false;
        return Ok(());
    }

    pub fn write_metadata(&self, config: &CorriesConfig) -> Result<()> {
        match self.stream_mode {
            StreamMode::Stdout => {
                println!("{}", config.metadata_dump());
            }
            StreamMode::File => {}
        }
        return Ok(());
    }

    fn data_matrix_to_stream_strings(&mut self) -> Result<()> {
        self.stream_strings.iter_mut().for_each(|line| line.clear());
        match self.string_conversion_mode {
            ToStringConversionMode::Scalar => {
                for value in self.data_matrix.iter() {
                    match value {
                        DataValue::Int(x) => {
                            self.stream_strings[0] += &format!("{:>width$}", x, width = self.width)
                        }
                        DataValue::Usize(x) => {
                            self.stream_strings[0] += &format!("{:>width$}", x, width = self.width)
                        }
                        DataValue::Float(x) => {
                            self.stream_strings[0] +=
                                &format!("{:>width$.*}", self.precision, x, width = self.width)
                        }
                        DataValue::String(x) => {
                            self.stream_strings[0] += &format!("{:>width$}", x, width = self.width)
                        }
                        DataValue::VectorFloat(_) => {
                            bail!("Tried writing a vector into a scalar output!")
                        }
                    };
                    self.stream_strings[0].push(self.delimiter);
                }
            }
            ToStringConversionMode::Vector => {
                for value in self.data_matrix.iter() {
                    for i in 0..self.stream_strings.len() {
                        match value {
                            DataValue::Int(x) => {
                                self.stream_strings[i] +=
                                    &format!("{:>width$}", x, width = self.width)
                            }
                            DataValue::Usize(x) => {
                                self.stream_strings[i] +=
                                    &format!("{:>width$}", x, width = self.width)
                            }
                            DataValue::Float(x) => {
                                self.stream_strings[i] +=
                                    &format!("{:>width$.*}", self.precision, x, width = self.width)
                            }
                            DataValue::String(x) => {
                                self.stream_strings[i] +=
                                    &format!("{:>width$}", x, width = self.width)
                            }
                            DataValue::VectorFloat(v) => {
                                self.stream_strings[i] += &format!(
                                    "{:>width$.*}",
                                    self.precision,
                                    v[i],
                                    width = self.width
                                )
                            }
                        };
                        self.stream_strings[i].push(self.delimiter);
                    }
                }
            }
        }
        self.stream_strings.iter_mut().for_each(|line| {
            line.pop();
        });
        return Ok(());
    }

    fn stream(&mut self) -> Result<()> {
        match self.stream_mode {
            StreamMode::Stdout => {
                if self.first_output {
                    let mut s = "#".to_string();
                    for name in self.data_names.iter() {
                        let name_string = format!("{:?}", name);
                        s += &format!("{:>width$}", name_string, width = self.width);
                        s.push(self.delimiter);
                    }
                    s.pop();
                    println!("{}", s);
                }
                for line in self.stream_strings.iter() {
                    println!("{}", line);
                }
            },
            // TODO: implement this!
            StreamMode::File => {
                todo!();
                // let path = Path::new(&(self.file_name.clone() + &format!("{:0w$}", self.output_counter, w=self.file_counter_width)));
                // path.display();
            }
        }
        return Ok(());
    }
}
