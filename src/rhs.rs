// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use color_eyre::Result;

use self::numflux::init_numflux;

mod numflux;

pub struct RHS {
    numflux: Box<dyn numflux::Numflux>,
}

impl RHS {
    pub fn new() -> Result<RHS> {
        return Ok(RHS { numflux: init_numflux() })
    }
}

