// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use self::output::Output;

mod output;

pub struct Writer<'a> {
    outputs: Vec<&'a mut Output<'a>>,
}

enum DataValue {
    Int(i32),
    Usize(usize),
    Double(f64),
    Vector(ndarray::ArrayD<f64>),
}
