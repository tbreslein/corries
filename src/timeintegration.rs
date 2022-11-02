// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

use color_eyre::{Result, eyre::bail};

use crate::{config::numericsconfig::NumericsConfig, mesh::Mesh, physics::Physics, rhs::Rhs};

use self::{rkf::RungeKuttaFehlberg, timestep::TimeStep};

mod rkf;
mod timestep;

pub enum DtKind {
    Init,
    CFL,
}

trait TimeSolver<const S: usize, const EQ: usize> {
    fn next_solution(
        &mut self,
        time: &mut TimeStep,
        u: &mut Physics<S, EQ>,
        rhs: &mut Rhs<S, EQ>,
        mesh: &Mesh<S>,
    ) -> Result<()>;
}

pub struct TimeIntegration<const S: usize, const EQ: usize> {
    pub time: TimeStep,
    solver: Box<dyn TimeSolver<S, EQ>>,
}

impl<const S: usize, const EQ: usize> TimeIntegration<S, EQ> {
    pub fn new(numericsconfig: &NumericsConfig, output_counter_max: usize) -> Self {
        return Self {
            time: TimeStep::new(numericsconfig, output_counter_max),
            solver: Box::new(RungeKuttaFehlberg::new()),
        };
    }
    pub fn next_solution(&mut self, u: &mut Physics<S, EQ>, rhs: &mut Rhs<S, EQ>, mesh: &Mesh<S>) -> Result<()> {
        self.solver.next_solution(&mut self.time, u, rhs, mesh)?;
        self.time.iter += 1;
        if self.time.iter >= self.time.iter_max {
            bail!("time.iter reached time.iter_max! time.iter = {}, time.iter_max = {}", self.time.iter, self.time.iter_max);
        }
        return Ok(());
    }
}
