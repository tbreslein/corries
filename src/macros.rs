// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

/// TODO
#[macro_export]
macro_rules! set_Physics_and_E {
    (Euler1DAdiabatic) => {
        type P = Euler1DAdiabatic<S>;
        const E: usize = 3;
    };
    (Euler1DIsot) => {
        type P = Euler1DIsot<S>;
        const E: usize = 2;
    };
    (Euler2DAdiabatic) => {
        type P = Euler2DAdiabatic<S>;
        const E: usize = 4;
    };
    (Euler2DIsot) => {
        type P = Euler2DIsot<S>;
        const E: usize = 3;
    };
}
