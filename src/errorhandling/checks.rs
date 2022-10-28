// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

/// Macro to check that a double is finite.
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
macro_rules! check_finite_multiple_doubles {
    ($($s:ident.$x:ident),*) => {
       $(check_finite_double!($s.$x);)*
    };
}

/// Macro to check that the elements of an `ArrayD` is finite.
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

/// Macro to check that an `ArrayD` is non-empty.
macro_rules! check_nonempty_arrayd {
    ($s:ident.$x:ident) => {
        ensure!(!$s.$x.is_empty(), "{0} is empty!", stringify!($x))
    };
}

/// Macro to check that an `ArrayD` is non-empty and its elements are finite.
macro_rules! check_nonempty_finite_arrayd {
    ($s:ident.$x:ident) => {
        check_nonempty_arrayd!($s.$x);
        check_finite_arrayd!($s.$x);
    };
}

/// Macro to check that, for each `ArrayD` in a list of `ArrayD`, that it is non-empty and their elements are finite.
macro_rules! check_nonempty_finite_multiple_arrayd {
    ($($s:ident.$x:ident),*) => {
        $(check_nonempty_finite_arrayd!($s.$x);)*
    };
}
