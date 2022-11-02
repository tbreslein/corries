// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

use crate::{physics::Physics, rhs::Rhs, mesh::Mesh};

mod rkf;

pub struct TimeIntegration<const S: usize, const EQ: usize> {
    pub time: TimeStep,
    solver: Box<dyn TimeSolver<S, EQ>>,
}

pub struct TimeStep {
    iter: usize,
    iter_max: usize,
    t: f64,
    t0: f64,
    t_end: f64,
    t_old: f64,
    dt: f64,
    dt_min: f64,
    dt_max: f64,
    dt_output: f64,
    t_next_output: f64,
    should_perform_output: f64,
    dt_cfl_param: f64,
    dt_kind: DtKind,
}

pub enum DtKind {
    CFL,
}

trait TimeSolver<const S: usize, const EQ: usize> {
    fn next_solution(&mut self, u: &mut Physics<S, EQ>, rhs: &mut Rhs<S, EQ>, mesh: &Mesh<S>);
}

struct RungeKuttaFehlberg<const S: usize, const EQ: usize> {

}

impl<const S: usize, const EQ: usize> TimeSolver<S, EQ> for RungeKuttaFehlberg<S, EQ> {
    fn next_solution(&mut self, u: &mut Physics<S, EQ>, rhs: &mut Rhs<S, EQ>, mesh: &Mesh<S>) {
        todo!();
    }
}

