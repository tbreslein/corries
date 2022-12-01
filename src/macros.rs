// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

/// TODO
#[macro_export]
macro_rules! set_Physics_and_E {
    (Euler1DIsot) => {
        type P = corries::physics::systems::euler1disot::Euler1DIsot<S>;
        const E: usize = 2;
    };
}
