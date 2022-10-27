// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use crate::{config::CorriesConfig, mesh::Mesh};
use self::output::Output;
use std::sync::{Arc, Mutex};

mod output;

pub struct Writer<'a> {
    outputs: Vec<Arc<Mutex<Output<'a>>>>,
}

#[derive(Debug)]
pub enum DataValue {
    Int(i32),
    Usize(usize),
    Float(f64),
    VectorFloat(ndarray::ArrayD<f64>),
    String(String),
}

impl Writer<'_> {
    pub fn new(config: &CorriesConfig, mesh: &Mesh) -> Self {
        let mut outputs = vec![];
        for outputconf in config.writerconf.iter() {
            outputs.push(Arc::new(Mutex::new(Output::new(&outputconf, &mesh))));
        }
        return Writer {outputs};
    }
}
