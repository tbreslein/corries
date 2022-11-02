// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

use color_eyre::Result;

use crate::{mesh::Mesh, physics::Physics, rhs::Rhs};

use super::{timestep::TimeStep, TimeSolver};

pub struct RungeKuttaFehlberg<const S: usize, const EQ: usize> {}

impl<const S: usize, const EQ: usize> TimeSolver<S, EQ> for RungeKuttaFehlberg<S, EQ> {
    fn next_solution(
        &mut self,
        time: &mut TimeStep,
        _u: &mut Physics<S, EQ>,
        _rhs: &mut Rhs<S, EQ>,
        _mesh: &Mesh<S>,
    ) -> Result<()> {
        todo!();
    }
}

impl<const S: usize, const EQ: usize> RungeKuttaFehlberg<S, EQ> {
    pub fn new() -> Self {
        return Self {};
    }
}
