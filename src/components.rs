// Copyright (c) 2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [CorriesComponents] type alias, as well as [Runner] trait.

use color_eyre::{eyre::Context, Result};
use crate::{DtKind, Mesh, NumFlux, Physics, Solver, State, TimeSolver, Writer};

/// This trait is written solely for the purpose of implementing this method on CorriesComponents,
/// which itself is just a type alias.
///
/// Refer to [CorriesComponents] for a detailed usage example.
pub trait Runner {
    /// Runs a [corries](crate) simulation.
    ///
    /// Refer to [CorriesComponents] for a detailed usage example.
    fn run_corries(&mut self) -> Result<()>;
}

/// A type alias around a tuple of the four components of a [corries](crate) simulation:
///
/// * [State]
/// * [Solver]
/// * [Mesh]
/// * [Writer]
///
/// The main purpose of this tuple is to be a return type for
/// [CorriesConfig::init_corries()](crate::CorriesConfig::init_corries()). With that, we can chain
/// that method into the [Runner::run_corries()] method to run the simulation.
pub type CorriesComponents<P, N, T, const E: usize, const S: usize> =
    (State<P, E, S>, Solver<P, N, T, E, S>, Mesh<S>, Writer);

impl<P, N, T, const E: usize, const S: usize> Runner for CorriesComponents<P, N, T, E, S>
where
    P: Physics<E, S>,
    N: NumFlux<E, S>,
    T: TimeSolver<P, E, S>,
{
    /// Runs a [corries](crate) simulation.
    ///
    /// This assumes that you already set your initial conditions in `u` and that it is fully
    /// up-to-date.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use corries::prelude::*;
    ///
    /// // Set up a simulation using isothermal Euler physics on a 100 cell mesh, using the Hll and
    /// // Runge-Kutta-Fehlberg schemes for solving the equations.
    /// const S: usize = 100;
    /// set_Physics_and_E!(Euler1DIsot);
    /// type N = Hll<E, S>;
    /// type T = RungeKuttaFehlberg<P, E, S>;
    ///
    /// // use the default config for Riemann tests
    /// let t_end = 0.5;
    /// let folder_name = "results";
    /// let file_name = "noh";
    ///
    /// // define the config instance
    /// CorriesConfig::default_riemann_test::<N, E, S>(t_end, folder_name, file_name)
    ///     .init_corries::<P, N, T, E,S>(|_, _, _| Ok(()))
    ///     .unwrap()
    ///     .run_corries()
    ///     .unwrap()
    /// ```
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
