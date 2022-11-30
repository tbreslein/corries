// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use std::fs::{create_dir_all, remove_dir_all, File};
use std::io::Write;
use std::path::Path;

use color_eyre::{eyre::WrapErr, Result};

use crate::physics::Physics;
use crate::timeintegration::TimeIntegration;
use crate::{
    config::outputconfig::{FormattingMode, OutputConfig, StreamMode, ToStringConversionMode},
    mesh::Mesh,
};

use super::data::{Data, StructAssociation};
use super::{Collectable, DataValue};

/// Struct that handles writing to output into a file or stdout
pub struct Output {
    /// Matrix of raw, unformatted data, still retaining their original data types
    pub data: Vec<Data>,

    /// Whether this is this `Output`'s first output in this run
    first_output: bool,

    /// Number of rows to print per output step
    rows: usize,

    /// Writes arbitrary data from `corries` objects to `string_matrix`
    string_conversion_mode: ToStringConversionMode,

    /// How the data is generally formatted per line, i.e. CSV vs TSV
    formatting_mode: FormattingMode,

    /// Handles writing into a stream
    stream_mode: StreamMode,

    /// File name for the file(s) to write file output to
    file_name: String,

    /// Ending of the file name (like .csv)
    file_name_ending: String,

    /// Index offset for printing spatially resolved data, useful for skipping ghost cells
    mesh_offset: usize,

    /// How many digits the file counter for file output needs
    output_counter_width: usize,

    /// Whether this `Output` should write metadata into its `Stream`
    should_print_metadata: bool,

    /// Precision of floating point numbers
    precision: usize,

    /// Width of the written strings
    width: usize,

    /// Intermediate string representation of the elements of the `data_matrix`
    delimiter: char,

    /// Leading character for comments
    leading_comment_symbol: String,
}

impl Output {
    /// Builds a new `color_eyre::Result<Output>` object.
    ///
    /// # Arguments
    ///
    /// * `outputconf` - `OutputConfig` containing configuration to build `Output` objects
    /// * `mesh` - provides information about the mesh in this simulation
    /// * `output_counter_max` - How many distinct outputs should be written during this run
    pub fn new<const S: usize>(outputconfig: &OutputConfig, mesh: &Mesh<S>, output_counter_max: usize) -> Result<Self> {
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
        let mesh_offset = if outputconfig.should_print_ghostcells {
            // TODO: turn these into global constants!
            mesh.imin
        } else {
            mesh.ixi_in
        };
        return Ok(Output {
            data: outputconfig
                .data_names
                .iter()
                .map(|name| Data::new::<S>(name, mesh_offset))
                .collect(),
            first_output: true,
            rows: match outputconfig.string_conversion_mode {
                ToStringConversionMode::Scalar => 1,
                ToStringConversionMode::Vector => {
                    if outputconfig.should_print_ghostcells {
                        mesh.n_all
                    } else {
                        mesh.n_comp
                    }
                },
            },
            string_conversion_mode: outputconfig.string_conversion_mode,
            formatting_mode: outputconfig.formatting_mode,
            stream_mode: outputconfig.stream_mode,
            file_name: if outputconfig.folder_name.ends_with('/') {
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
            },
            file_name_ending: match outputconfig.formatting_mode {
                FormattingMode::TSV => ".tsv".to_string(),
                FormattingMode::CSV => ".csv".to_string(),
            },
            mesh_offset,
            output_counter_width: f64::log10(output_counter_max as f64) as usize + 1,
            should_print_metadata: outputconfig.should_print_metadata,
            precision: outputconfig.precision,
            width: match outputconfig.formatting_mode {
                FormattingMode::TSV => 12usize.min(outputconfig.precision + 7),
                FormattingMode::CSV => 0,
            },
            delimiter: match outputconfig.formatting_mode {
                FormattingMode::TSV => ' ',
                FormattingMode::CSV => ',',
            },
            leading_comment_symbol: match outputconfig.formatting_mode {
                FormattingMode::TSV => "#".to_string(),
                FormattingMode::CSV => "".to_string(),
            },
        });
    }

    /// Goes through `self.datanames` and pulls data from the other arguments writing that data
    /// into `self.data_matrix`.
    ///
    /// # Arguments
    ///
    /// * `u` - Provides data for the state of the simulation
    /// * `time` - Provides data on the time coordinates
    /// * `mesh` - Provides mesh data
    pub fn update_data<const S: usize, const EQ: usize>(
        &mut self,
        u: &Physics<S, EQ>,
        time: &TimeIntegration<S, EQ>,
        mesh: &Mesh<S>,
    ) -> Result<()> {
        self.data.iter_mut().try_for_each(|data| match data.association {
            StructAssociation::Mesh => mesh.collect_data(data, self.mesh_offset),
            StructAssociation::Physics => u.collect_data(data, self.mesh_offset),
            StructAssociation::TimeStep => time.timestep.collect_data(data, self.mesh_offset),
        })?;
        return Ok(());
    }

    /// Makes this `Output` object write data in `self.data_matrix` into a stream.
    pub fn write_output(&mut self, output_counter: usize) -> Result<()> {
        match self.stream_mode {
            StreamMode::File => {
                match self.string_conversion_mode {
                    ToStringConversionMode::Scalar => {
                        let full_path_string = self.file_name.clone() + &self.file_name_ending;
                        let path = Path::new(&full_path_string);
                        let mut file = File::options()
                            .write(true)
                            .append(!self.first_output)
                            .create(true)
                            .open(path)
                            .wrap_err_with(|| format!("Failed to open file: {}!", path.display()))?;

                        self.write_into_stream(&mut file)?;
                    },
                    ToStringConversionMode::Vector => {
                        let full_path_string = self.file_name.clone()
                            + &format!("{:0w$}", output_counter, w = self.output_counter_width)
                            + &self.file_name_ending;
                        let path = Path::new(&full_path_string);
                        let mut file = File::options()
                            .write(true)
                            .append(false)
                            .create(true)
                            .open(path)
                            .wrap_err_with(|| format!("Failed to open file: {}!", path.display()))?;

                        self.write_into_stream(&mut file)?;
                    },
                };
            },
            StreamMode::Stdout => {
                self.write_into_stream(&mut std::io::stdout())?;
            },
        };
        self.first_output = false;
        return Ok(());
    }

    /// Take a buffer, and write this object's data into it
    fn write_into_stream<T: std::io::Write>(&self, buffer: &mut T) -> Result<()> {
        if (self.string_conversion_mode == ToStringConversionMode::Scalar && self.first_output)
            || self.string_conversion_mode == ToStringConversionMode::Vector
        {
            write!(buffer, "{}", self.get_header())?;
        }
        for j in 0..self.rows {
            match self.formatting_mode {
                // TSV output needs on space of padding to align with the # at the beginning of the header
                FormattingMode::TSV => {
                    write!(buffer, " ")?;
                },
                FormattingMode::CSV => {},
            };

            // Write the first piece of data by itself, so that afterwards we can always print the delimiter as
            // preceding to the next piece of data.
            self.format_data_to_buffer(buffer, &self.data[0], j)?;
            for data in self.data.iter().skip(1) {
                write!(buffer, "{}", self.delimiter)?;
                self.format_data_to_buffer(buffer, data, j)?;
            }
            writeln!(buffer)?;
        }
        return Ok(());
    }

    /// Format a piece of data (at index j, if it is indexable) and write it into the buffer
    fn format_data_to_buffer<T: std::io::Write>(&self, buffer: &mut T, data: &Data, j: usize) -> Result<()> {
        return match &data.payload {
            DataValue::String(x) => {
                write!(buffer, "{:>width$}", x, width = self.width).wrap_err("Unable to write to buffer!")
            },
            DataValue::Usize(x) => {
                write!(buffer, "{:>width$}", x, width = self.width).wrap_err("Unable to write to buffer!")
            },
            DataValue::Float(x) => write!(buffer, "{:>width$.*e}", self.precision, x, width = self.width)
                .wrap_err("Unable to write to buffer!"),
            DataValue::VectorFloat(x) => write!(buffer, "{:>width$.*e}", self.precision, x[j], width = self.width)
                .wrap_err("Unable to write to buffer!"),
        };
    }

    /// Makes this `Output` object write meta_data into a stream.
    pub fn write_metadata<const S: usize>(&self, meta_data: &str) -> Result<()> {
        if self.should_print_metadata {
            match self.stream_mode {
                StreamMode::Stdout => {
                    println!("{}", meta_data);
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
                    file.write_all(meta_data.as_bytes())
                        .wrap_err_with(|| format!("Failed to write to file: {}!", path.display()))?;
                },
            }
        }
        return Ok(());
    }

    /// Returns the header for writing output.
    fn get_header(&self) -> String {
        let mut s = self
            .data
            .iter()
            .fold(self.leading_comment_symbol.clone(), |mut acc, name| {
                acc += &format!("{:>width$}", format!("{}", name), width = self.width);
                acc.push(self.delimiter);
                return acc;
            });
        s.pop();
        s.push('\n');
        return s;
    }
}
