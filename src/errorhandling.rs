// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Contains everything regarding error handling, and exports the [Validation] trait as well as the
//! [checks] module.

use color_eyre::Result;

#[macro_use]
pub mod checks;

/// Trait for all structs that can validate themselves.
///
/// These structs need to implement the `validate(&self) -> color_eyre::Result<()>` method. This
/// method's purpose is to make sure the fields of a given struct are coherent and adhere to
/// rules specific to this struct.
///
/// For example, a `Mesh`'s coordinate arrays should be finite and non-empty, so its `validate`
/// method should check its arrays for these properties.
pub trait Validation {
    /// Make sure that the fields of `&self` are coherent and adhere to struct specific internal rules.
    fn validate(&self) -> Result<()>;
}
