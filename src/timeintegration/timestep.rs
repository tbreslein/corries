// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

use color_eyre::{eyre::bail, Result};

use crate::{config::numericsconfig::NumericsConfig, mesh::Mesh, physics::Physics};

use super::DtKind;

pub struct TimeStep {
    pub iter: usize,
    pub iter_max: usize,
    pub t: f64,
    pub t_end: f64,
    pub dt: f64,
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
            t_end: numericsconfig.t_end,
            dt: 0.0,
            dt_min: numericsconfig.dt_min,
            dt_max: numericsconfig.dt_max,
            dt_cfl_param: numericsconfig.dt_cfl_param,
            dt_kind: DtKind::Init,
            dt_output: (numericsconfig.t_end - numericsconfig.t0) / output_counter_max as f64,
            t_next_output: numericsconfig.t0,
        };
    }

    pub fn calc_dt_expl<const S: usize, const EQ: usize>(
        &mut self,
        u: &mut Physics<S, EQ>,
        // rhs: &Rhs<S, EQ>,
        mesh: &Mesh<S>,
    ) -> Result<()> {
        // let dt_cfl_pair = (u.calc_dt_cfl(self.dt_cfl_param, &mesh), DtKind::CFL);
        //TODO: source term time steps
        self.dt = u.calc_dt_cfl(self.dt_cfl_param, mesh)?;
        if self.dt < self.dt_min {
            bail!(
                "Time step width dt dipped below dt_min! Got dt = {}, dt_min = {}",
                self.dt,
                self.dt_min
            );
        }
        self.dt_kind = DtKind::Cfl;
        return Ok(());
    }

    pub fn cap_dt(&mut self) {
        self.dt = self.dt.min(self.dt_max);
        if (self.dt + self.t) / self.t_next_output > 1.0 {
            self.dt = self.t_next_output - self.t;
        }
    }
}
