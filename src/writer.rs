// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::{eyre::bail, Result};
use ndarray::{s, ArrayD};

use self::output::Output;
use crate::{
    config::{outputconfig::DataName, CorriesConfig},
    mesh::Mesh,
};

mod output;

pub struct Writer {
    pub outputs: Vec<Output>,
}

#[derive(Debug)]
pub enum DataValue {
    Int(i32),
    Usize(usize),
    Float(f64),
    VectorFloat(ArrayD<f64>),
    String(String),
}

impl Writer {
    pub fn new(config: &CorriesConfig, mesh: &Mesh) -> Self {
        let mut outputs = vec![];
        for outputconf in config.writerconf.iter() {
            outputs.push(Output::new(outputconf, mesh));
        }
        return Writer { outputs };
    }

    pub fn update_data_matrices(&mut self, mesh: &Mesh) -> Result<()> {
        // TODO: can I safely thread this loop?
        for output in self.outputs.iter_mut() {
            output.update_data_matrix(mesh)?;
        }
        return Ok(());
    }

    pub fn write_output(&mut self) -> Result<()> {
        // TODO: can I safely thread this loop?
        for output in self.outputs.iter_mut() {
            output.write_output()?;
        }
        return Ok(());
    }

    pub fn write_metadata(&mut self, config: &CorriesConfig) -> Result<()> {
        for output in self.outputs.iter_mut() {
            if output.should_print_metadata {
                output.write_metadata(config)?;
            };
        }
        return Ok(());
    }
}

/// Describes objects can write to an `Output`
pub trait Write {
    fn collect_data(
        &self,
        name: &DataName,
        value: &mut DataValue,
        mesh_offset: usize,
    ) -> Result<()>;

    fn write_vector(
        &self,
        field: &ArrayD<f64>,
        value: &mut DataValue,
        mesh_offset: usize,
    ) -> Result<()> {
        return match value {
            DataValue::VectorFloat(v) => {
                if mesh_offset == 0 {
                    v.assign(field);
                    Ok(())
                } else {
                    v.assign(&field.slice(s![mesh_offset..field.len() - mesh_offset]));
                    Ok(())
                }
            }
            _ => bail!(
                "Tried assigning vector field to non-vector DataValue!\nfield: {}\nvalue: {:?}",
                field,
                value
            ),
        };
    }
}
