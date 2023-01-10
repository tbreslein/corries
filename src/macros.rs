// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports useful macros for public use

/// Expands to a type alias `P` for the type of physics you are using, and a constant `E` that
/// represents the number of equations in that system.
/// This macro is used to make sure that E is always set correctly for your type of physics, while
/// also giving you a useful type alias to make your code look more generic than it actually is.
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
