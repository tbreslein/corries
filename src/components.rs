// Copyright (c) 2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

use color_eyre::{eyre::Context, Result};

use crate::{DtKind, Mesh, NumFlux, Physics, Solver, State, TimeSolver, Writer};

/// todo
pub trait Runner {
    /// todo
    fn run_corries(&mut self) -> Result<()>;
}

/// todo
pub type CorriesComponents<P, N, T, const E: usize, const S: usize> =
    (State<P, E, S>, Solver<P, N, T, E, S>, Mesh<S>, Writer);

impl<P, N, T, const E: usize, const S: usize> Runner for CorriesComponents<P, N, T, E, S>
where
    P: Physics<E, S>,
    N: NumFlux<E, S>,
    T: TimeSolver<P, E, S>,
{
    fn run_corries(&mut self) -> Result<()> {
        // source: https://stackoverflow.com/questions/41207306/mutable-reference-to-a-tuple-as-input-parameter
        // This is where I learned what the ref keyword is...
        let (ref mut u, ref mut solver, ref mesh, ref mut writer) = *self;
        loop {
            if solver.timestep.t >= solver.timestep.t_next_output - solver.timestep.dt_min {
                solver.timestep.t_next_output += solver.timestep.dt_output;
                writer
                    .update_data(&u.cent, &solver.timestep, mesh)
                    .context("Calling writer.update_data_matrices in run_corries")?;
                writer
                    .write_output()
                    .context("Calling wirter.write_output in run_corries")?;
            }

            if solver.timestep.t + solver.timestep.dt_min * 0.01 >= solver.timestep.t_end {
                break;
            }

            if let err @ Err(_) = solver.next_solution(u, &mesh) {
                solver.timestep.dt_kind = DtKind::ErrorDump;
                writer
                    .update_data(&u.cent, &solver.timestep, mesh)
                    .context("Calling writer.update_data_matrices during the error dump in run_corries")?;
                writer
                    .write_output()
                    .context("Calling writer.write_output during the error dump in run_corries")?;
                return err;
            }
        }
        Ok(())
    }
}
