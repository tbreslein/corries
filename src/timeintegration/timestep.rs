// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

use crate::config::numericsconfig::NumericsConfig;

use super::DtKind;

pub struct TimeStep {
    pub iter: usize,
    pub iter_max: usize,
    pub t: f64,
    t0: f64,
    pub t_end: f64,
    t_old: f64,
    dt: f64,
    dt_min: f64,
    dt_max: f64,
    dt_cfl_param: f64,
    dt_kind: DtKind,

    /// Time step width between outputs
    pub dt_output: f64,

    /// At which time to perform the next output
    pub t_next_output: f64,
}

impl TimeStep {
    pub fn new(numericsconfig: &NumericsConfig, output_counter_max: usize) -> Self {
        return Self {
            iter: 0,
            iter_max: numericsconfig.iter_max,
            t: numericsconfig.t0,
            t0: numericsconfig.t0,
            t_end: numericsconfig.t_end,
            t_old: numericsconfig.t0,
            dt: 0.0,
            dt_min: numericsconfig.dt_min,
            dt_max: numericsconfig.dt_max,
            dt_cfl_param: numericsconfig.dt_cfl_param,
            dt_kind: DtKind::Init,
            dt_output: (numericsconfig.t_end - numericsconfig.t0) / output_counter_max as f64,
            t_next_output: numericsconfig.t0,
        };
    }
}
