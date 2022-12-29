// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [RungeKuttaFehlberg] struct

use color_eyre::{
    eyre::{ensure, Context},
    Result,
};
use ndarray::{Array2, Array3, Axis};

use crate::{
    config::CorriesConfig,
    errorhandling::Validation,
    // errorhandling::Validation,
    mesh::Mesh,
    physics::Physics,
    rhs::Rhs,
    update_everything_from_cons,
    NumFlux,
    TimeIntegrationConfig,
};

use self::butchertableau::ButcherTableau;

use super::{timestep::TimeStep, TimeSolver};

mod butchertableau;

/// Struct for solving the time integration step using Runge-Kutta-Fehlberg methods.
pub struct RungeKuttaFehlberg<P: Physics + Validation> {
    /// Butcher Tableau for the Runge-Kutta method
    bt: ButcherTableau,

    /// Stores the different intermediate solutions to construct a new solution from
    k_bundle: Array3<f64>,

    /// Newly calculated error in this asc step
    err_new: f64,

    /// Error from the previous asc step
    err_old: f64,

    /// The number of asc steps during this iteration
    n_asc: usize,

    /// Relative tolerance for automated step control
    asc_relative_tolerance: f64,

    /// Absolute tolerance for automated step control
    asc_absolute_tolerance: f64,

    /// Timestep "friction" for automated step control
    asc_timestep_friction: f64,

    /// Whether the solution in the current iteration has been accepted
    solution_accepted: bool,

    /// Conservative variables for the low-order solution, used by automated step control
    u_cons_low: Array2<f64>,

    /// Primitive variables from the last time step; automated step control needs this to
    /// be saved at the beginning of each iteration to roll the state back from time to time.
    u_prim_old: Array2<f64>,

    /// Stores the full intermediate solution
    utilde: P,
}

impl<P: Physics + Validation + 'static> TimeSolver<P> for RungeKuttaFehlberg<P> {
    /// Constructs a new [RungeKuttaFehlberg] object
    ///
    /// # Arguments
    ///
    /// * `rkfconfig` - Configuration specifically for [RungeKuttaFehlberg] objects
    /// * `physicsconfig` - Configuration for [Physics] objects, needed because `utilde`
    fn new<const E: usize, const S: usize>(config: &CorriesConfig, u: &P) -> Result<Self> {
        let (bt, asc_relative_tolerance, asc_absolute_tolerance, asc_timestep_friction) =
            match &config.numerics_config.time_integration_config {
                TimeIntegrationConfig::Rkf(rkf_config) => (
                    ButcherTableau::new(rkf_config),
                    rkf_config.asc_relative_tolerance,
                    rkf_config.asc_absolute_tolerance,
                    rkf_config.asc_timestep_friction,
                ),
            };
        let order = bt.order;
        let mut utilde = P::new(&config.physics_config);
        utilde.assign(u);
        return Ok(Self {
            bt,
            k_bundle: Array3::zeros([order, E, S]),
            err_new: 0.0,
            err_old: 0.0,
            n_asc: 0,
            asc_relative_tolerance,
            asc_absolute_tolerance,
            asc_timestep_friction,
            solution_accepted: true,
            u_cons_low: Array2::zeros((E, S)),
            u_prim_old: Array2::zeros((E, S)),
            utilde,
        });
    }

    fn next_solution<N: NumFlux, const E: usize, const S: usize>(
        &mut self,
        time: &mut TimeStep,
        u: &mut P,
        rhs: &mut Rhs<P, N, S>,
        mesh: &Mesh<S>,
    ) -> Result<()> {
        time.iter += 1;
        time.calc_dt_expl(u, mesh)
            .context("time.calc_dt_expl at the beginning of RungeKuttaFehlberg::next_solution")?;
        time.cap_dt();

        if self.bt.asc {
            self.u_prim_old.assign(&u.prim());
            self.solution_accepted = false;
            self.n_asc = 0;

            while !self.solution_accepted {
                time.dt = self.calc_rkf_solution::<N, E, S>(time.dt, u, rhs, mesh).context(
                    "RungeKuttaFehlberg::calc_rkf_solution in the while loop in RungeKuttaFehlberg::next_solution",
                )?;
                if self.err_new < 1.0 {
                    time.cap_dt();
                    self.solution_accepted = true;
                }
                self.n_asc += 1;
            }
        }

        // NOTE: for embedded methods, this call calculates a solution with the finalised this->dt since the solution is
        // already considered accepted. This is necessary, since this->dt might have been capped again.
        let _ = self
            .calc_rkf_solution::<N, E, S>(time.dt, u, rhs, mesh)
            .context("RungeKuttaFehlberg::calc_rkf_solution at the end of RungeKuttaFehlberg::next_solution")?;
        u.assign_cons(&self.utilde.cons());
        update_everything_from_cons(u, &mut rhs.boundary_west, &mut rhs.boundary_east, mesh);
        u.validate()
            .context("Calling u.update_everything_from_prim at the end of RungeKuttaFehlberg::calc_rkf_solution")?;
        time.t += time.dt;
        return Ok(());
    }
}

impl<P: Physics + Validation + 'static> RungeKuttaFehlberg<P> {
    /// Calculates a single solution with an RKF method.
    ///
    /// # Arguments
    ///
    /// * `dt_in` - Input time step width
    /// * `u` - Input [Physics] state
    /// * `rhs` - Solves the right-hand side
    /// * `mesh` - Information about spatial properties
    fn calc_rkf_solution<N: NumFlux, const E: usize, const S: usize>(
        &mut self,
        dt_in: f64,
        u: &mut P,
        rhs: &mut Rhs<P, N, S>,
        mesh: &Mesh<S>,
    ) -> Result<f64> {
        let mut dt_out = dt_in;
        self.k_bundle.fill(0.0);
        // calculate the k_bundle entries
        for q in 0..self.bt.order {
            self.utilde.assign(u);
            for p in 0..q {
                self.utilde.assign_cons(
                    &(&self.utilde.cons() - dt_out * self.bt.a[[p, q]] * &self.k_bundle.index_axis(Axis(0), p)).view(),
                );
            }
            rhs.update::<E>(&mut self.utilde, mesh)
                .context("Calling rhs.update while calculating k_bundle in RungeKuttaFehlberg::calc_rkf_solution")?;
            let mut k_bundle_q = self.k_bundle.index_axis_mut(Axis(0), q);
            k_bundle_q.assign(&rhs.full_rhs);
        }

        // calculate high order solution
        self.utilde.assign_cons(&u.cons());
        for q in 0..self.bt.order {
            self.utilde.assign_cons(
                &(&self.utilde.cons() - dt_out * self.bt.b_high[q] * &self.k_bundle.index_axis(Axis(0), q)).view(),
            );
        }

        if !self.solution_accepted {
            // calculate low order solution
            self.u_cons_low.assign(&u.cons());
            for q in 0..self.bt.order {
                self.utilde.assign_cons(
                    &(&self.u_cons_low - dt_out * self.bt.b_low[q] * &self.k_bundle.index_axis(Axis(0), q)).view(),
                );
            }

            // calc err_new
            self.err_new = {
                let mut err_max = 0.0f64;
                for j in 0..E {
                    for i in 2..S - 2 {
                        err_max = err_max.max(
                            ((self.utilde.cons_entry(j, i) - self.u_cons_low[[j, i]])
                                / (self.asc_relative_tolerance
                                    * (self.utilde.cons_entry(j, i) + self.asc_absolute_tolerance)))
                                .abs(),
                        );
                    }
                }
                err_max
            };

            // adjust dt
            dt_out = {
                // Step 1
                // P-Controller step: New time step estimate based on the error between the high and low order
                // solutions, and the order of the high order scheme. If the error is > 0, this suggests that dt
                // needs to be reduced. Otherwise, increase it.
                let dt_temp = if self.err_new > 0.0 {
                    0.9 * dt_out * (-1.0 * self.err_new.ln() / self.bt.order as f64).exp()
                } else {
                    4.0 * dt_out
                };

                // Step 2
                // PI-Controller step: Controlled increase or decrease of the time step based on the error of the
                // previous try-out. If the error was also > 0 in the last try-out, reduce it further and save
                // in dt_new. Otherwise, leave it as is and save in dt_new.
                let mut dt_new = if self.err_old > 0.0 {
                    dt_temp * (-1.0 * self.asc_timestep_friction * (self.err_old / self.err_new).ln()).exp()
                } else {
                    dt_temp
                };

                // Step 3
                // Final adjustments
                // If the error is small, and the new time step is much larger than the previous dt_out, cap dt_new at
                // 4.0 * dt_old.
                if self.err_new < 1.0 {
                    dt_new = dt_new.min(4.0 * dt_out);
                } else {
                    // in case err_old is significantly larger than err_new, dt_new might be >
                    // dt_out.
                    // This sets it back to the P-Controller value.
                    if dt_new > dt_out {
                        dt_new = dt_temp;
                    }
                    // either way, decrease dt by setting it to either dt_new or, at most, 0.25 *
                    // dt_out
                    dt_new = dt_new.max(0.25 * dt_out);
                }

                dt_new
            };
            self.err_old = self.err_new;
        }

        self.validate()
            .context("Validating RungeKuttaFehlberg at the end of RungeKuttaFehlberg::calc_rkf_solution")?;
        return Ok(dt_out);
    }
}

impl<P: Physics + Validation> Validation for RungeKuttaFehlberg<P> {
    fn validate(&self) -> Result<()> {
        ensure!(
            self.u_cons_low.fold(true, |acc, x| acc && x.is_finite()),
            "RungeKuttaFehlberg::u_cons_low must be finite! Got: {}",
            self.u_cons_low
        );
        self.utilde
            .validate()
            .context("Validating RungeKuttaFehlberg::utilde in RungeKuttaFehlberg::validate()")?;
        return Ok(());
    }
}
