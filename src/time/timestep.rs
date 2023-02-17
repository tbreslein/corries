// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [TimeStep]

use super::DtKind;
use crate::{
    config::numericsconfig::NumericsConfig,
    data::{Data, DataName, StructAssociation},
    mesh::Mesh,
    state::Physics,
    Collectable, DataValue, State,
};
use color_eyre::{eyre::bail, Result};

/// Information about the time coordinate and related data
#[derive(Debug, Clone, Default, PartialEq)]
pub struct TimeStep {
    /// Current iteration counter for the simulation
    pub iter: usize,

    /// Maximum value for `iter`; the simulation will error out if `iter` hits this value
    pub iter_max: usize,

    /// Current time coordinate
    pub t: f64,

    /// Time coordinate at which the simulation terminates
    pub t_end: f64,

    /// Current time step width
    pub dt: f64,

    /// Minimum value for `dt`; the simulation will error out if `dt` dips below this value
    pub dt_min: f64,

    /// Maximum value for `dt`; simply a safety measure if a good value for `dt_cfl_param` is not
    /// known yet for this simulation
    dt_max: f64,

    /// CFL safety parameter
    dt_cfl_param: f64,

    /// The limiting factor this iteration's time step width
    pub dt_kind: DtKind,

    /// Time step width between outputs
    pub dt_output: f64,

    /// At which time to perform the next output
    pub t_next_output: f64,
}

unsafe impl Send for TimeStep {}
unsafe impl Sync for TimeStep {}

impl TimeStep {
    /// Constructs a new [TimeStep] object
    ///
    /// # Arguments
    ///
    /// * `numericsconfig` - Configuration for everything regarding numerics
    /// * `output_counter_max` - The number of outputs this simulation produces
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    /// use corries::timestep::TimeStep;
    ///
    /// // set up constants
    /// set_Physics_and_E!(Euler1DAdiabatic);
    /// const S: usize = 100;
    /// type N = Kt<E,S>;
    ///
    /// let t_end = 1.0;
    ///
    /// // define the config instance
    /// let numerics_config = NumericsConfig::default_riemann_test::<N, E, S>(t_end);
    ///
    /// let time_step = TimeStep::new(&numerics_config, 2);
    /// ```
    pub fn new(numericsconfig: &NumericsConfig, output_counter_max: usize) -> Self {
        Self {
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
        }
    }

    /// Calculates the explicit time step width and updates it in `self.dt` and `self.dt_kind`.
    ///
    /// # Arguments
    ///
    /// * `u` - The current [State]
    /// * `mesh` - Information about spatial properties
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use corries::prelude::*;
    /// use corries::timestep::TimeStep;
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
    /// let config = CorriesConfig::default_riemann_test::<N, E, S>(t_end, folder_name, file_name);
    /// let (mut u, mut solver, mesh, _) = config.init_corries::<P, N, T, E, S>(|_,_,_| Ok(())).unwrap();
    /// solver.timestep.calc_dt_expl(&mut u, &mesh).unwrap();
    /// ```
    pub fn calc_dt_expl<P: Physics<E, S>, const E: usize, const S: usize>(
        &mut self,
        u: &mut State<P, E, S>,
        // rhs: &Rhs<S, EQ>,
        mesh: &Mesh<S>,
    ) -> Result<()> {
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
        Ok(())
    }

    /// Caps the time step width to an upper limit, for example so that the time coordinate does
    /// not exceed the point where the next output should be produced.
    pub fn cap_dt(&mut self) {
        self.dt = self.dt.min(self.dt_max);
        if (self.dt + self.t) / self.t_next_output > 1.0 {
            self.dt = self.t_next_output - self.t;
        }
    }
}

impl Collectable for TimeStep {
    fn collect_data(&self, data: &mut Data, _: usize) -> Result<()> {
        match (data.association, data.name) {
            (StructAssociation::TimeStep, DataName::Iter) => data.payload = DataValue::Usize(self.iter),
            (StructAssociation::TimeStep, DataName::T) => data.payload = DataValue::Float(self.t),
            (StructAssociation::TimeStep, DataName::Dt) => data.payload = DataValue::Float(self.dt),
            (StructAssociation::TimeStep, DataName::DtKind) => {
                data.payload = DataValue::String(format!("{}", self.dt_kind))
            },
            (StructAssociation::TimeStep, x) => bail!("Tried associating {:?} with Time!", x),
            (StructAssociation::Mesh, x) | (StructAssociation::Physics, x) => {
                bail!("name.association() for {:?} returned {:?}", x, data.association)
            },
        };
        Ok(())
    }
}
