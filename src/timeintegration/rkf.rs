// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! TODO

use color_eyre::Result;
use ndarray::{Array2, Array3, Axis};

use crate::{
    config::{numericsconfig::RkfConfig, physicsconfig::PhysicsConfig},
    mesh::Mesh,
    physics::{init_physics, Physics},
    rhs::Rhs,
};

use self::butchertableau::ButcherTableau;

use super::{timestep::TimeStep, TimeSolver};

mod butchertableau;

pub struct RungeKuttaFehlberg<const S: usize, const EQ: usize> {
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
    utilde: Physics<S, EQ>,

    /// Temporary holder for dt, needed for embedded methods
    dt_temp: f64,
}

impl<const S: usize, const EQ: usize> TimeSolver<S, EQ> for RungeKuttaFehlberg<S, EQ> {
    fn next_solution(
        &mut self,
        time: &mut TimeStep,
        u: &mut Physics<S, EQ>,
        rhs: &mut Rhs<S, EQ>,
        mesh: &Mesh<S>,
    ) -> Result<()> {
        time.iter += 1;
        self.n_asc = 0;
        time.calc_dt_expl(u, mesh)?;
        time.cap_dt();

        if self.bt.asc {
            self.u_prim_old.assign(&u.prim);
            self.solution_accepted = false;
            self.dt_temp = time.dt;
        } else {
            self.solution_accepted = true;
        }

        while !self.solution_accepted {
            self.dt_temp = self.calc_rkf_solution(time.dt, u, rhs, mesh)?;
            if self.err_new < 1.0 {
                time.dt = self.dt_temp;
                time.cap_dt();
                self.solution_accepted = true;
            }

            // reset before calculating the final rkf solution or trying again depending on
            // solution_accepted
            u.prim.assign(&self.u_prim_old);
            u.update_cons();
            u.update_derived_variables();
            self.n_asc += 1;
        }

        // Calculate the final rkf solution.
        // TODO:
        // For embedded methods, this might actually be a redundant call to this method, since u
        // should already have that solution from the last iteration of the while loop.
        // In an earlier iteration skipping this call for embedded methods did not work, so I
        // should play around with this again, since this call can be quite expensive.
        let _ = self.calc_rkf_solution(time.dt, u, rhs, mesh)?;
        time.t += time.dt;
        return Ok(());
    }
}

impl<const S: usize, const EQ: usize> RungeKuttaFehlberg<S, EQ> {
    pub fn new(rkfconfig: &RkfConfig, physicsconfig: &PhysicsConfig) -> Self {
        let bt = ButcherTableau::new(rkfconfig);
        let order = bt.order;
        return Self {
            bt,
            k_bundle: Array3::zeros([order, EQ, S]),
            err_new: 0.0,
            err_old: 0.0,
            n_asc: 0,
            asc_relative_tolerance: rkfconfig.asc_relative_tolerance,
            asc_absolute_tolerance: rkfconfig.asc_absolute_tolerance,
            asc_timestep_friction: rkfconfig.asc_timestep_friction,
            solution_accepted: false,
            u_cons_low: Array2::zeros((EQ, S)),
            u_prim_old: Array2::zeros((EQ, S)),
            utilde: init_physics(physicsconfig),
            dt_temp: 0.0,
        };
    }

    fn calc_rkf_solution(
        &mut self,
        dt_in: f64,
        u: &mut Physics<S, EQ>,
        rhs: &mut Rhs<S, EQ>,
        mesh: &Mesh<S>,
    ) -> Result<f64> {
        let mut dt_out = dt_in;
        self.k_bundle.fill(0.0);
        // calculate the k_bundle entries
        for q in 0..self.bt.order {
            self.utilde.cons.assign(&u.cons);
            for p in 0..q {
                self.utilde
                    .cons
                    .assign(&(&self.utilde.cons - dt_out * &self.bt.a[[p, q]] * &self.k_bundle.index_axis(Axis(0), p)));
            }
            rhs.update(u, mesh);
            let mut k_bundle_q = self.k_bundle.index_axis_mut(Axis(0), q);
            k_bundle_q.assign(&rhs.full_rhs);
        }

        // calculate high order solution
        self.utilde.cons.assign(&u.cons);
        for q in 0..self.bt.order {
            self.utilde
                .cons
                .assign(&(&self.utilde.cons - dt_out * &self.bt.b_high[q] * &self.k_bundle.index_axis(Axis(0), q)));
        }

        if !self.solution_accepted {
            // calculate low order solution
            self.u_cons_low.assign(&u.cons);
            for q in 0..self.bt.order {
                self.utilde
                    .cons
                    .assign(&(&self.u_cons_low - dt_out * &self.bt.b_high[q] * &self.k_bundle.index_axis(Axis(0), q)));
            }

            // calc err_new
            self.err_new = {
                let mut err_max = 0.0f64;
                for j in 0..EQ {
                    for i in 2..S - 2 {
                        err_max = err_max.max(
                            (self.utilde.cons[[j, i]]
                                - self.u_cons_low[[j, i]]
                                    / (self.asc_relative_tolerance
                                        * (self.utilde.cons[[j, i]] + self.asc_absolute_tolerance).abs()))
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
            // NOTE: is this necessary?
            // self.u_cons_low.fill(0.0);
        }

        u.cons.assign(&self.utilde.cons);
        rhs.update_physics(u, mesh);
        return Ok(dt_out);
    }
}
