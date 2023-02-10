// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

/// Macro to check that a double is finite.
///
/// # Examples
///
/// ```ignore
/// use color_eyre::{Result, eyre::ensure};
///
/// fn get_finite() -> Result<f64> {
///     check_finite_double!(1.0);
///     Ok(())
/// }
/// fn get_inf() -> Result<f64> {
///     check_finite_double!(f64::INFINITY);
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
macro_rules! check_finite_double {
    ($s:ident.$x:ident) => {
        ensure!(
            $s.$x.is_finite(),
            "{0} turned non-finite! Got: {0} = {1}",
            stringify!($x),
            $s.$x
        )
    };
}

/// Macro to check that a list of doubles is finite.
///
/// # Examples
///
/// ```ignore
/// use color_eyre::{Result, eyre::ensure};
///
/// fn get_finites() -> Result<f64> {
///     check_finite_multiple_doubles!(1.0, 2.0, 3.0);
///     Ok(())
/// }
/// fn get_non_finites() -> Result<f64> {
///     check_finite_multiple_doubles!(f64::INFINITY, f64::NAN);
///     Ok(())
/// }
///
/// fn main() {
///     let succeeds = get_finites();
///     let errors = get_non_finites();
///     assert!(succeeds.is_ok());
///     assert!(errors.is_err());
/// }
/// ```
macro_rules! check_finite_multiple_doubles {
    ($($s:ident.$x:ident),*) => {
       $(check_finite_double!($s.$x);)*
    };
}

/// Macro to check that the elements of an `Array` are finite.
///
/// # Examples
///
/// ```ignore
/// use color_eyre::{Result, eyre::ensure};
/// use ndarray::Array1;
///
/// fn get_finites() -> Result<f64> {
///     check_finite_arrayd!(Array1::zeros(5));
///     Ok(())
/// }
/// fn get_infs() -> Result<f64> {
///     check_finite_arrayd!(Array1::from_elem(5, f64::INFINITY));
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
macro_rules! check_finite_arrayd {
    ($s:ident.$x:ident) => {
        ensure!(
            $s.$x.fold(true, |acc, y| acc & y.is_finite()),
            "{0} turned non-finite! Got: {0} = {1}",
            stringify!($x),
            $s.$x
        )
    };
}

/// Macro to check that an `Array` is non-empty.
///
/// # Examples
///
/// ```ignore
/// use color_eyre::{Result, eyre::ensure};
/// use ndarray::Array1;
///
/// fn get_non_empty() -> Result<f64> {
///     check_nonempty_arrayd!(Array1::zeros(5));
///     Ok(())
/// }
/// fn get_empty() -> Result<f64> {
///     check_nonempty_arrayd!(Array1::zeros(0));
///     Ok(())
/// }
///
/// fn main() {
///     let succeeds = get_non_empty();
///     let errors = get_empty();
///     assert!(succeeds.is_ok());
///     assert!(errors.is_err());
/// }
/// ```
macro_rules! check_nonempty_arrayd {
    ($s:ident.$x:ident) => {
        ensure!(!$s.$x.is_empty(), "{0} is empty!", stringify!($x))
    };
}

/// Macro to check that an `Array` is non-empty and its elements are finite.
///
/// # Examples
///
/// ```ignore
/// use color_eyre::{Result, eyre::ensure};
/// use ndarray::Array1;
///
/// fn get_non_empty() -> Result<f64> {
///     check_nonempty_finite_arrayd!(Array1::zeros(5));
///     Ok(())
/// }
/// fn get_empty() -> Result<f64> {
///     check_nonempty_finite_arrayd!(Array1::zeros(0));
///     Ok(())
/// }
/// fn get_infs() -> Result<f64> {
///     check_nonempty_finite_arrayd!(Array1::from_elem(5, f64::INFINITY));
///     Ok(())
/// }
///
/// fn main() {
///     let succeeds = get_non_empty();
///     let errors = get_empty();
///     let errors_too = get_infs();
///     assert!(succeeds.is_ok());
///     assert!(errors.is_err());
///     assert!(errors_too.is_err());
/// }
/// ```
macro_rules! check_nonempty_finite_arrayd {
    ($s:ident.$x:ident) => {
        check_nonempty_arrayd!($s.$x);
        check_finite_arrayd!($s.$x);
    };
}

/// Macro to check that, for each `Array` in a list of `Array`, that it is non-empty and their elements are finite.
///
/// # Examples
///
/// ```ignore
/// use color_eyre::{Result, eyre::ensure};
/// use ndarray::Array1;
///
/// fn get_finite_non_empty() -> Result<f64> {
///     check_nonempty_finite_multiple_arrayd!(Array1::ones(7), Array1::zeros(5));
///     Ok(())
/// }
/// fn get_empty() -> Result<f64> {
///     check_nonempty_finite_multiple_arrayd!(Array1::zeros(3), Array1::zeros(0));
///     Ok(())
/// }
/// fn get_infs() -> Result<f64> {
///     check_nonempty_finite_multiple_arrayd!(Array1::ones(5), Array1::from_elem(5, f64::INFINITY));
///     Ok(())
/// }
///
/// fn main() {
///     let succeeds = get_finite_non_empty();
///     let errors = get_empty();
///     let errors_too = get_infs();
///     assert!(succeeds.is_ok());
///     assert!(errors.is_err());
///     assert!(errors_too.is_err());
/// }
/// ```
macro_rules! check_nonempty_finite_multiple_arrayd {
    ($($s:ident.$x:ident),*) => {
        $(check_nonempty_finite_arrayd!($s.$x);)*
    };
}
