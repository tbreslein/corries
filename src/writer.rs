// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::{eyre::bail, Result};
use ndarray::{s, Array1, ArrayView1};

use self::output::Output;
use crate::{
    config::{outputconfig::DataName, CorriesConfig},
    mesh::Mesh,
    physics::Physics,
};

mod output;

/// Wrapper for `Vec<Output>` that facilitates writing output into streams.
pub struct Writer {
    /// Contains the objects that handle writing into single streams.
    pub outputs: Vec<Output>,

    /// Current output counter
    output_counter: usize,

    /// Whether we should perform an output at this point
    pub should_perform_output: bool,
}

/// Wrapper enum around the different datatypes that can be written into an output stream.
#[derive(Debug)]
pub enum DataValue {
    Int(i32),
    Usize(usize),
    Float(f64),
    VectorFloat(Array1<f64>),
    String(String),
}

impl Writer {
    /// Builds a new `Result<Writer>` object.
    ///
    /// # Arguments
    ///
    /// * `config` - Configuration for a corries simulation
    /// * `mesh` - The mesh for this simulation
    /// * `outputer_count_max` - How many output steps this simulation should be going through
    pub fn new<const S: usize>(config: &CorriesConfig, mesh: &Mesh<S>) -> Result<Self> {
        let mut outputs = vec![];
        for outputconf in config.writerconfig.iter() {
            let output = Output::new(outputconf, mesh, config.output_counter_max)?;
            outputs.push(output);
        }
        return Ok(Writer {
            outputs,
            output_counter: 0,
            should_perform_output: true,
        });
    }

    /// Loops through `self.outputs` and calls their `update_data_matrix` methods.
    ///
    /// # Arguments
    ///
    /// * `mesh` - Provides mesh data
    pub fn update_data_matrices<const S: usize, const EQ: usize>(
        &mut self,
        mesh: &Mesh<S>,
        u: &Physics<S, EQ>,
    ) -> Result<()> {
        // TODO: can I safely thread this loop?
        for output in self.outputs.iter_mut() {
            output.update_data_matrix(mesh, u)?;
        }
        return Ok(());
    }

    /// Loops through `self.outputs` and calls their `write_output` methods.
    pub fn write_output(&mut self) -> Result<()> {
        // TODO: can I safely thread this loop?
        for output in self.outputs.iter_mut() {
            output.write_output(self.output_counter)?;
        }
        self.output_counter += 1;
        return Ok(());
    }

    /// Loops through `self.outputs` and calls their `write_metadata` methods if they are
    /// configured to do so.
    pub fn write_metadata<const S: usize>(&mut self, config: &CorriesConfig) -> Result<()> {
        for output in self.outputs.iter_mut() {
            if output.should_print_metadata {
                output.write_metadata::<S>(config)?;
            };
        }
        return Ok(());
    }
}

/// Describes objects can write to an `Output`
pub trait CorriesWrite {
    /// Collects the `value` corresponding to the `name`.
    ///
    /// # Arguments
    ///
    /// * `name` - The identifier of this piece of data
    /// * `value` - Where to write the data to
    /// * `mesh_offset` - In case of writing vector data, how many cells should be skipped at the
    /// beginning and the end of the vector
    fn collect_data(&self, name: &DataName, value: &mut DataValue, mesh_offset: usize) -> Result<()>;

    /// Writes a vector-like `field` to a `value`.
    ///
    /// # Arguments
    ///
    /// * `field` - The vector-like piece of data
    /// * `value` - Where to write the data to
    /// * `mesh_offset` - In case of writing vector data, how many cells should be skipped at the
    /// beginning and the end of the vector
    fn write_vector(&self, field: &ArrayView1<f64>, value: &mut DataValue, mesh_offset: usize) -> Result<()> {
        return match value {
            DataValue::VectorFloat(v) => {
                if mesh_offset == 0 {
                    v.assign(field);
                    Ok(())
                } else {
                    v.assign(&field.slice(s![mesh_offset..field.len() - mesh_offset]));
                    Ok(())
                }
            },
            _ => bail!(
                "Tried assigning vector field to non-vector DataValue!\nfield: {}\nvalue: {:?}",
                field,
                value
            ),
        };
    }
}
