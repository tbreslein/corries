// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports useful macros for public use

/// Expands to a type alias `P` for the type of physics you are using, and a constant `E` that
/// represents the number of equations in that system.
/// This macro is used to make sure that E is always set correctly for your type of physics, while
/// also giving you a useful type alias to make your code look more generic than it actually is.
///
/// NOTE: You need to have set the constant S, which sets how many cells the mesh will have, before
/// calling this macro, otherwise it will not compile!
///
/// # Arguments
///
/// Only accepts either of:
///
/// * `Euler1DAdiabatic`
/// * `Euler1DIsot`
///
/// # Examples
///
/// ```
/// use corries::prelude::*;
/// use std::any::TypeId;
///
/// // Set up adiabatic 1d Euler physics
/// const S: usize = 100;
/// set_Physics_and_E!(Euler1DAdiabatic);
/// assert_eq!(TypeId::of::<P>(), TypeId::of::<Euler1DAdiabatic<S>>());
/// assert_eq!(E, 3);
/// ```
///
/// ```
/// use corries::prelude::*;
/// use std::any::TypeId;
///
/// // Set up isothermal 1d Euler physics
/// const S: usize = 100;
/// set_Physics_and_E!(Euler1DIsot);
/// assert_eq!(TypeId::of::<P>(), TypeId::of::<Euler1DIsot<S>>());
/// assert_eq!(E, 2);
/// ```
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
}
