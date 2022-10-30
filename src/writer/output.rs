// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use std::fs::{create_dir_all, remove_dir_all, File};
use std::io::Write;
use std::path::Path;

use color_eyre::{
    eyre::{bail, WrapErr},
    Result,
};
use ndarray::Array1;

use crate::{
    config::{
        outputconfig::{
            DataName, DataType, FormatterMode, OutputConfig, StreamMode, StructAssociation, ToStringConversionMode,
        },
        CorriesConfig,
    },
    mesh::Mesh,
};

use super::{CorriesWrite, DataValue};

/// Struct that handles writing to output into a file or stdout
pub struct Output {
    /// Index offset for printing spatially resolved data, useful for skipping ghost cells
    mesh_offset: usize,

    /// Whether this is this `Output`'s first output in this run
    first_output: bool,

    /// Current output counter
    output_counter: usize,

    /// How many digits the file counter for file output needs
    output_counter_width: usize,

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

    /// Leading character for comments
    leading_comment_symbol: String,

    /// Vector of strings ready to be written to a stream
    stream_strings: Vec<String>,

    /// Writes arbitrary data from `corries` objects to `string_matrix`
    string_conversion_mode: ToStringConversionMode,

    /// Handles writing into a stream
    pub stream_mode: StreamMode,

    /// File name for the file(s) to write file output to
    file_name: String,

    /// Ending of the file name (like .csv)
    file_name_ending: String,

    /// Identifiers for the data being written to the stream
    data_names: Vec<DataName>,
}

impl Output {
    /// Builds a new `color_eyre::Result<Output>` object.
    ///
    /// # Arguments
    ///
    /// * `outputconf` - `OutputConfig` containing configuration to build `Output` objects
    /// * `mesh` - provides information about the mesh in this simulation
    /// * `output_count_max` - How many distinct outputs should be written during this run
    pub fn new<const S: usize>(outputconfig: &OutputConfig, mesh: &Mesh<S>, output_count_max: usize) -> Result<Self> {
        let rows = match outputconfig.string_conversion_mode {
            ToStringConversionMode::Scalar => 1,
            ToStringConversionMode::Vector => {
                if outputconfig.should_print_ghostcells {
                    mesh.n_all
                } else {
                    mesh.n_comp
                }
            },
        };
        let mut data_matrix = vec![];
        for name in outputconfig.data_names.iter() {
            match name.datatype() {
                DataType::Int => data_matrix.push(DataValue::Int(0)),
                DataType::Usize => data_matrix.push(DataValue::Usize(0)),
                DataType::Float => data_matrix.push(DataValue::Float(0.0)),
                DataType::String => data_matrix.push(DataValue::String("".to_string())),
                DataType::VectorFloat => data_matrix.push(DataValue::VectorFloat(Array1::zeros(S))),
            };
        }

        let file_name = if outputconfig.folder_name.ends_with('/') {
            outputconfig.folder_name.clone() + &outputconfig.file_name.clone()
        } else {
            // make sure there is a / between folder and file names
            outputconfig.folder_name.clone() + "/" + &outputconfig.file_name.clone()
        } + {
            if let ToStringConversionMode::Vector = outputconfig.string_conversion_mode {
                // vector output files are seperated into multiple files, and have to enumerated.
                // for example, the format of one of those file names could be:
                // /path/to/my/output_0420.csv
                // where the _ seperates the filename root and the output counter, and output
                // counter is an integer with leading zeros as padding to make sure you can
                // correctly sort the files.
                "_"
            } else {
                ""
            }
        };

        if let StreamMode::File = outputconfig.stream_mode {
            if outputconfig.should_clear_out_folder && Path::new(&outputconfig.folder_name).is_dir() {
                remove_dir_all(&outputconfig.folder_name)
                    .wrap_err_with(|| format!("Failed to remove directory {}", outputconfig.folder_name))?;
            }
            if !Path::new(&outputconfig.folder_name).is_dir() {
                create_dir_all(&outputconfig.folder_name)
                    .wrap_err_with(|| format!("Failed to create directory {}", outputconfig.folder_name))?;
            }
        }

        return Ok(Output {
            mesh_offset: if outputconfig.should_print_ghostcells {
                mesh.imin
            } else {
                mesh.ixi_in
            },
            first_output: true,
            output_counter: 0,
            output_counter_width: f64::log10(output_count_max as f64) as usize + 1,
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
            leading_comment_symbol: match outputconfig.formatter_mode {
                FormatterMode::TSV => "#".to_string(),
                FormatterMode::CSV => "".to_string(),
            },
            stream_strings: vec!["".to_string(); rows],
            string_conversion_mode: outputconfig.string_conversion_mode,
            stream_mode: outputconfig.stream_mode,
            // folder_name: outputconfig.folder_name.clone(),
            file_name,
            file_name_ending: match outputconfig.formatter_mode {
                FormatterMode::TSV => ".tsv".to_string(),
                FormatterMode::CSV => ".csv".to_string(),
            },
            data_names: outputconfig.data_names.clone(),
        });
    }

    /// Goes through `self.datanames` and pulls data from the other arguments writing that data
    /// into `self.data_matrix`.
    ///
    /// # Arguments
    ///
    /// * `mesh` - `Mesh` object to pull data from
    pub fn update_data_matrix<const S: usize>(&mut self, mesh: &Mesh<S>) -> Result<()> {
        for (i, name) in self.data_names.iter().enumerate() {
            match name.association() {
                StructAssociation::Mesh => mesh.collect_data(name, &mut self.data_matrix[i], self.mesh_offset),
            }?;
        }
        return Ok(());
    }

    /// Makes this `Output` object write data in `self.data_matrix` into a stream.
    pub fn write_output(&mut self) -> Result<()> {
        self.data_matrix_to_stream_strings()?;
        self.stream()?;
        self.first_output = false;
        return Ok(());
    }

    /// Makes this `Output` object write meta data in `self.data_matrix` into a stream.
    pub fn write_metadata<const S: usize>(&self, config: &CorriesConfig) -> Result<()> {
        match self.stream_mode {
            StreamMode::Stdout => {
                println!("{}", config.metadata_dump::<S>());
            },
            StreamMode::File => {
                let full_path_string = self.file_name.clone() + "_metadata.dat";
                let path = Path::new(&full_path_string);
                let mut file = File::options()
                    .write(true)
                    .append(false)
                    .create(true)
                    .open(path)
                    .wrap_err_with(|| format!("Failed to open file: {}!", path.display()))?;
                file.write_all(config.metadata_dump::<S>().as_bytes())
                    .wrap_err_with(|| format!("Failed to write to file: {}!", path.display()))?;
            },
        }
        return Ok(());
    }

    /// Transforms `self.data_matrix` into `self.stream_strings`.
    fn data_matrix_to_stream_strings(&mut self) -> Result<()> {
        self.stream_strings.iter_mut().for_each(|line| line.clear());
        match self.string_conversion_mode {
            ToStringConversionMode::Scalar => {
                for value in self.data_matrix.iter() {
                    match value {
                        DataValue::Int(x) => self.stream_strings[0] += &format!("{:>width$}", x, width = self.width),
                        DataValue::Usize(x) => self.stream_strings[0] += &format!("{:>width$}", x, width = self.width),
                        DataValue::Float(x) => {
                            self.stream_strings[0] += &format!("{:>width$.*}", self.precision, x, width = self.width)
                        },
                        DataValue::String(x) => self.stream_strings[0] += &format!("{:>width$}", x, width = self.width),
                        DataValue::VectorFloat(_) => {
                            bail!("Tried writing a vector into a scalar output!")
                        },
                    };
                    self.stream_strings[0].push(self.delimiter);
                }
            },
            ToStringConversionMode::Vector => {
                for value in self.data_matrix.iter() {
                    for i in 0..self.stream_strings.len() {
                        match value {
                            DataValue::Int(x) => {
                                self.stream_strings[i] += &format!("{:>width$}", x, width = self.width)
                            },
                            DataValue::Usize(x) => {
                                self.stream_strings[i] += &format!("{:>width$}", x, width = self.width)
                            },
                            DataValue::Float(x) => {
                                self.stream_strings[i] +=
                                    &format!("{:>width$.*}", self.precision, x, width = self.width)
                            },
                            DataValue::String(x) => {
                                self.stream_strings[i] += &format!("{:>width$}", x, width = self.width)
                            },
                            DataValue::VectorFloat(v) => {
                                self.stream_strings[i] +=
                                    &format!("{:>width$.*}", self.precision, v[i], width = self.width)
                            },
                        };
                        self.stream_strings[i].push(self.delimiter);
                    }
                }
            },
        }
        self.stream_strings.iter_mut().for_each(|line| {
            line.pop();
            line.push('\n');
        });
        return Ok(());
    }

    /// Writes `self.stream_strings` into a stream depending on `self.stream_mode`.
    fn stream(&mut self) -> Result<()> {
        match self.stream_mode {
            StreamMode::Stdout => {
                if self.first_output {
                    print!("{}", self.get_header());
                }
                for line in self.stream_strings.iter() {
                    print!("{}", line);
                }
            },
            StreamMode::File => match self.string_conversion_mode {
                ToStringConversionMode::Scalar => {
                    let full_path_string = self.file_name.clone() + &self.file_name_ending;
                    let path = Path::new(&full_path_string);
                    let mut file = File::options()
                        .write(true)
                        .append(!self.first_output)
                        .create(true)
                        .open(path)
                        .wrap_err_with(|| format!("Failed to open file: {}!", path.display()))?;

                    if self.first_output {
                        file.write_all(self.get_header().as_bytes())
                            .wrap_err_with(|| format!("Failed to write to file: {}!", path.display()))?;
                    }
                    file.write_all(self.stream_strings[0].as_bytes())
                        .wrap_err_with(|| format!("Failed to write to file: {}!", path.display()))?;
                },
                ToStringConversionMode::Vector => {
                    let full_path_string =
                        self.file_name.clone() + &self.get_output_count_string() + &self.file_name_ending;
                    let path = Path::new(&full_path_string);
                    let mut file = File::options()
                        .write(true)
                        .append(false)
                        .create(true)
                        .open(path)
                        .wrap_err_with(|| format!("Failed to open file: {}!", path.display()))?;

                    file.write_all(self.get_header().as_bytes())
                        .wrap_err_with(|| format!("Failed to write to file: {}!", path.display()))?;
                    for line in self.stream_strings.iter() {
                        file.write_all(line.as_bytes())
                            .wrap_err_with(|| format!("Failed to write to file: {}!", path.display()))?;
                    }
                },
            },
        }
        return Ok(());
    }

    /// Returns the header for writing output.
    fn get_header(&self) -> String {
        let mut s = self.leading_comment_symbol.clone();
        for name in self.data_names.iter() {
            let name_string = format!("{:?}", name);
            s += &format!("{:>width$}", name_string, width = self.width);
            s.push(self.delimiter);
        }
        s.pop();
        s.push('\n');
        return s;
    }

    /// Returns the `String` containing the output counter padded with leading zeros.
    fn get_output_count_string(&self) -> String {
        return format!("{:0w$}", self.output_counter, w = self.output_counter_width);
    }
}
