// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use self::numflux::init_numflux;

mod numflux;

pub struct Rhs {
    _numflux: Box<dyn numflux::Numflux>,
}

impl Rhs {
    pub fn new() -> Self {
        return Rhs {
            _numflux: init_numflux(),
        };
    }
}
