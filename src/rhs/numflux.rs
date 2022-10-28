// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

pub trait Numflux {}

pub fn init_numflux() -> Box<dyn Numflux> {
    return Box::new(Hll {});
}

pub struct Hll {}

impl Numflux for Hll {}
