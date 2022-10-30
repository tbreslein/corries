// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

/// Length scale in CGS units
const LENGTH_CGS: f64 = 0.01;

/// Mass scale in CGS units
const MASS_CGS: f64 = 0.001;

/// Enum for the different kinds of unit systems available
#[derive(Debug, Copy, Clone)]
pub enum UnitsMode {
    SI,
    CGS,
}

/// Struct that serves to know the type units the simulation uses, and exposes methods to convert
/// between different units.
pub struct Units {
    mode: UnitsMode,
}

impl Units {
    /// Builds a new `Units` object.
    ///
    /// # Arguments
    ///
    /// * `mode` - The type of units used when setting up the simulation
    pub fn new(mode: UnitsMode) -> Self {
        return Units { mode };
    }

    /// Takes the argument representing a length, and converts it to SI units
    pub fn convert_length_to_si(&self, input: f64) -> f64 {
        return match self.mode {
            UnitsMode::SI => input,
            UnitsMode::CGS => input * LENGTH_CGS,
        };
    }

    /// Takes the argument representing a time, and converts it to SI units
    pub fn convert_time_to_si(&self, input: f64) -> f64 {
        return match self.mode {
            UnitsMode::SI => input,
            UnitsMode::CGS => input,
        };
    }

    /// Takes the argument representing a mass, and converts it to SI units
    pub fn convert_mass_to_si(&self, input: f64) -> f64 {
        return match self.mode {
            UnitsMode::SI => input,
            UnitsMode::CGS => input * MASS_CGS,
        };
    }

    /// Takes the argument representing a length, and converts it to CGS units
    pub fn convert_length_to_cgs(&self, input: f64) -> f64 {
        let si_input = self.convert_length_to_si(input);
        return match self.mode {
            UnitsMode::SI => si_input,
            UnitsMode::CGS => si_input / LENGTH_CGS,
        };
    }

    /// Takes the argument representing a time, and converts it to CGS units
    pub fn convert_time_to_cgs(&self, input: f64) -> f64 {
        let si_input = self.convert_time_to_si(input);
        return match self.mode {
            UnitsMode::SI => si_input,
            UnitsMode::CGS => si_input,
        };
    }

    /// Takes the argument representing a mass, and converts it to CGS units
    pub fn convert_mass_to_cgs(&self, input: f64) -> f64 {
        let si_input = self.convert_mass_to_si(input);
        return match self.mode {
            UnitsMode::SI => si_input,
            UnitsMode::CGS => si_input / MASS_CGS,
        };
    }
}
