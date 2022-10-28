// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::Result;

use self::numflux::init_numflux;

mod numflux;

pub struct Rhs {
    _numflux: Box<dyn numflux::Numflux>,
}

impl Rhs {
    pub fn new() -> Result<Rhs> {
        return Ok(Rhs {
            _numflux: init_numflux(),
        });
    }
}
