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
        const E: usize = P::NUM_EQ;
    };
    (Euler1DIsot) => {
        type P = Euler1DIsot<S>;
        const E: usize = P::NUM_EQ;
    };
}

/// Macro to check that each double in a list is finite.
///
/// # Examples
///
/// ```
/// use color_eyre::{Result, eyre::ensure};
/// use corries::check_finite_double;
///
/// fn get_finite() -> Result<()> {
///     check_finite_double!(1.0_f64);
///     Ok(())
/// }
/// fn get_inf() -> Result<()> {
///     check_finite_double!(2.2_f64, f64::INFINITY, 1.0_f64);
///     Ok(())
/// }
///
/// fn main() {
///     let succeeds = get_finite();
///     let errors = get_inf();
///     assert!(succeeds.is_ok());
///     assert!(errors.is_err());
/// }
/// ```
#[macro_export]
macro_rules! check_finite_double {
    ($($s:ident.$x:ident),*) => {
       $(ensure!(
            $s.$x.is_finite(),
            "{0} turned non-finite! Got: {0} = {1}",
            stringify!($x),
            $s.$x
        );)*
    };
    ($($x:expr),*) => {
        $(
            ensure!(
                $x.is_finite(),
                "{0} turned non-finite! Got: {0} = {1}",
                stringify!($x),
                $x);
         )*
    };
}

/// Macro to check that each `Array` in a list is finite.
///
/// # Examples
///
/// ```
/// use color_eyre::{Result, eyre::ensure};
/// use corries::check_finite_arrayd;
/// use ndarray::Array1;
///
/// fn get_finites() -> Result<()> {
///     check_finite_arrayd!(Array1::<f64>::zeros(5));
///     Ok(())
/// }
/// fn get_infs() -> Result<()> {
///     check_finite_arrayd!(Array1::<f64>::ones(3), Array1::<f64>::from_elem(5, f64::INFINITY));
///     Ok(())
/// }
///
/// fn main() {
///     let succeeds = get_finites();
///     let errors = get_infs();
///     assert!(succeeds.is_ok());
///     assert!(errors.is_err());
/// }
/// ```
#[macro_export]
macro_rules! check_finite_arrayd {
    ($($s:ident.$x:ident),*) => {
       $(ensure!(
            $s.$x.fold(true, |acc, y| acc && y.is_finite()),
            "{0} turned non-finite! Got: {0} = {1}",
            stringify!($x),
            $s.$x
        );)*
    };
    ($($x:expr),*) => {
       $(ensure!(
            $x.fold(true, |acc, y| acc && y.is_finite()),
            "{0} turned non-finite! Got: {0} = {1}",
            stringify!($x),
            $x
        );)*
    };
}

/// Macro to check that each double in a list is finite.
///
/// # Examples
///
/// ```
/// use color_eyre::{Result, eyre::ensure};
/// use corries::check_positive_double;
///
/// fn get_positive() -> Result<()> {
///     check_positive_double!(1.0_f64);
///     Ok(())
/// }
/// fn get_zero() -> Result<()> {
///     check_positive_double!(2.2_f64, 0.0_f64, 1.0_f64);
///     Ok(())
/// }
///
/// fn main() {
///     let succeeds = get_positive();
///     let errors = get_zero();
///     assert!(succeeds.is_ok());
///     assert!(errors.is_err());
/// }
/// ```
#[macro_export]
macro_rules! check_positive_double {
    ($($s:ident.$x:ident),*) => {
       $(ensure!(
            $s.$x > 0.0,
            "{0} is non-positive! Got: {0} = {1}",
            stringify!($x),
            $s.$x
        );)*
    };
    ($($x:expr),*) => {
        $(
            ensure!(
                $x > 0.0,
                "{0} is non-positive! Got: {0} = {1}",
                stringify!($x),
                $x);
         )*
    };
}

/// Macro to check that each `Array` has only positive non-zero elements.
///
/// # Examples
///
/// ```
/// use color_eyre::{Result, eyre::ensure};
/// use corries::check_positive_arrayd;
/// use ndarray::Array1;
///
/// fn get_ones() -> Result<()> {
///     check_positive_arrayd!(Array1::<f64>::ones(5));
///     Ok(())
/// }
/// fn get_zeros() -> Result<()> {
///     check_positive_arrayd!(Array1::<f64>::ones(2), Array1::<f64>::zeros(2));
///     Ok(())
/// }
///
/// fn main() {
///     let succeeds = get_ones();
///     let errors = get_zeros();
///     assert!(succeeds.is_ok());
///     assert!(errors.is_err());
/// }
/// ```
#[macro_export]
macro_rules! check_positive_arrayd {
    ($($s:ident.$x:ident),*) => {
       $(ensure!(
            $s.$x.fold(true, |acc, y| acc && y > &0.0),
            "{0} turned non-positive! Got: {0} = {1}",
            stringify!($x),
            $s.$x
        );)*
    };
    ($($x:expr),*) => {
       $(ensure!(
            $x.fold(true, |acc, y| acc && y > &0.0),
            "{0} turned non-positive! Got: {0} = {1}",
            stringify!($x),
            $x
        );)*
    };
}
