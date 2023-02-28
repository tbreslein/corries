// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

//! Exports the [CorriesConfig] structs and its nested structs for configuring Corries simulations.

use crate::{components::*, errorhandling::Validation, initfuncs::InitFn, prelude::*};
use color_eyre::{eyre::Context, Result};
pub use meshconfig::*;
pub use numericsconfig::*;
pub use outputconfig::*;
pub use physicsconfig::*;
use serde::Serialize;

pub mod meshconfig;
pub mod numericsconfig;
pub mod outputconfig;
pub mod physicsconfig;

/// Enumerates the different boundary conditions.
///
/// Defaults to [NoGradients](BoundaryMode::NoGradients)
#[derive(Debug, Serialize, Clone, Default, PartialEq)]
pub enum BoundaryMode {
    /// Set of custom boundary conditions applied to seperate equations identified by their
    /// equation indexes
    Custom(Vec<(usize, CustomBoundaryMode)>),

    /// Sets no-gradients boundaries for all equations (default)
    #[default]
    NoGradients,
}

unsafe impl Send for BoundaryMode {}
unsafe impl Sync for BoundaryMode {}

/// Enumerates the possible custom boundary conditions.
///
/// Defaults to [NoGradients](CustomBoundaryMode::NoGradients)
#[derive(Debug, Serialize, Clone, Copy, Default, PartialEq, Eq)]
pub enum CustomBoundaryMode {
    /// Extrapolate the values near the boundary into the ghost cells
    Extrapolate,

    /// Specialised version of Extrapolate for mass density in the Kepler case
    ExtrapolateDensityKepler,

    /// Specialised version of Extrapolate for eta velocity in the Kepler case
    ExtrapolateEtaVelocityKepler,

    /// Like NoGradients, but multiplies the value in the ghost cell with a very small number
    NearZero,

    /// No gradients condition; copies the value closest to the boundary into the ghost cells.
    #[default]
    NoGradients,

    /// Reflects outgoing flow, and applies NoGradients to incoming flow
    OutflowNoGradients,

    /// Reflects outgoing flow, and applies Extrapolate to incoming flow
    OutFlowExtrapolate,

    /// Like NoGradients, but switches the sign of the values in the ghost cells
    Reflecting,
}

unsafe impl Send for CustomBoundaryMode {}
unsafe impl Sync for CustomBoundaryMode {}

/// Struct that carries the full configuration info for a simulation.
///
/// This struct is used in the beginning of a run to initialise all the runtime-objects that are
/// used throughout the simulation.
#[derive(Debug, Serialize, Clone)]
pub struct CorriesConfig {
    /// Whether to print the Corries banner to stdout
    pub print_banner: bool,

    /// Config for [Mesh](crate::mesh::Mesh) objects
    pub mesh_config: MeshConfig,

    /// Config for [Physics](crate::state::Physics) objects
    pub physics_config: PhysicsConfig,

    /// boundary condition on the west border of the computational area
    pub boundary_condition_west: BoundaryMode,

    /// boundary condition on the east border of the computational area
    pub boundary_condition_east: BoundaryMode,

    /// Config for everything related to numerics
    pub numerics_config: NumericsConfig,

    /// The number of outputs to write during the simulation, not counting output for the initial
    /// state
    pub output_counter_max: usize,

    /// Config for [Writer](crate::writer::Writer) objects
    pub writer_config: Vec<OutputConfig>,
}

unsafe impl Send for CorriesConfig {}
unsafe impl Sync for CorriesConfig {}

impl CorriesConfig {
    /// Sets up a default [CorriesConfig] that can be used by most Riemann tests.
    ///
    /// Check the asserts in the example, as well as the docs for the following methods to see the
    /// full config:
    ///
    /// * [MeshConfig::default_riemann_test()]
    /// * [PhysicsConfig::default()]
    /// * [NumericsConfig::default_riemann_test()]
    /// * [OutputConfig::default_stdout()]
    /// * [OutputConfig::default_file()]
    ///
    /// # Arguments
    ///
    /// * `t_end` - The time coordinate at which the simulation is terminated
    /// * `folder_name` - The folder to write the file output to
    /// * `file_name` - The base name of the files that output is being written to
    ///
    /// # Examples
    ///
    /// ```
    /// use corries::prelude::*;
    ///
    /// // set up constants
    /// set_Physics_and_E!(Euler1DAdiabatic);
    /// const S: usize = 100;
    /// type N = Hll<E,S>;
    ///
    /// let t_end = 0.5;
    /// let folder_name = "results";
    /// let file_name = "noh";
    ///
    /// // define the config instance
    /// let config = CorriesConfig::default_riemann_test::<N, E, S>(t_end, folder_name, file_name);
    /// assert_eq!(config.print_banner, false);
    /// assert_eq!(config.boundary_condition_west, BoundaryMode::NoGradients);
    /// assert_eq!(config.boundary_condition_east, BoundaryMode::NoGradients);
    /// assert_eq!(config.output_counter_max, 1);
    /// ```
    pub fn default_riemann_test<N: NumFlux<E, S> + 'static, const E: usize, const S: usize>(
        t_end: f64,
        folder_name: &str,
        file_name: &str,
    ) -> Self {
        Self {
            print_banner: false,
            mesh_config: MeshConfig::default_riemann_test(),
            physics_config: PhysicsConfig::default(),
            boundary_condition_west: BoundaryMode::NoGradients,
            boundary_condition_east: BoundaryMode::NoGradients,
            numerics_config: NumericsConfig::default_riemann_test::<N, E, S>(t_end),
            output_counter_max: 1,
            writer_config: vec![
                OutputConfig::default_stdout(),
                OutputConfig::default_file(folder_name, file_name, E),
            ],
        }
    }

    /// Initialises all objects needed to run a corries simulation.
    ///
    /// Apart from the `self` argument, the important bits that also help configuring the simulation
    /// are the template Parameters. For example, the type you pass as the first template argument
    /// determines the type of [Physics] used throughout the whole simulation!
    ///
    /// In addition to that, you also need to pass a function that takes a [State], a [Solver], and
    /// a [Mesh], and returns a [color_eyre::Result].
    /// That function is used to apply the initial conditions to the [State] object.
    ///
    /// In that `init_fn` you only need to concern yourself with setting initial values for the
    /// [State::cent.prim] field for the [State] object.
    /// After calling `init_fn`, this `init_corries` method will make sure that [State::west] and
    /// [State::east] are initialised, as well as that the boundary conditions are applied.
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
    /// CorriesConfig::default_riemann_test::<N, E, S>(t_end, folder_name, file_name)
    ///     .init_corries::<P, N, T, E, S>(|u, solver, mesh| {
    ///         /* ... some complex initial conditions for u ... */
    ///         Ok(())
    ///     }).unwrap();
    /// ```
    pub fn init_corries<P, N, T, const E: usize, const S: usize>(
        &self,
        init_fn: InitFn<P, N, T, E, S>,
    ) -> Result<CorriesComponents<P, N, T, E, S>>
    where
        P: Physics<E, S> + 'static,
        N: NumFlux<E, S>,
        T: TimeSolver<P, E, S>,
    {
        self.validate().context("Validating CorriesConfig")?;

        let mesh = Mesh::<S>::new(&self.mesh_config).context("Constructing Mesh")?;
        let mut u = State::<P, E, S>::new(&self.physics_config);
        let mut solver = Solver::<P, N, T, E, S>::new(self, &mesh).context("Constructing Solver")?;
        let mut writer = Writer::new::<S>(self, &mesh).context("Constructing Writer")?;

        if writer.print_banner {
            println!("# ****************************************");
            println!("# Corries - corrosive Riemann solver ");
            println!("# ");
            println!("# Version: {}", env!("CARGO_PKG_VERSION"));
            println!("# Copyright (c) 2022-2023");
            println!("# Author: tbreslein <github.com/tbreslein>");
            println!("# License: MIT");
            println!("# ****************************************");
        }
        writer
            .write_metadata::<S>()
            .context("Calling writer.write_metadata in run_corries")?;

        init_fn(&mut u, &mut solver, &mesh).context("Calling init_fn in CorriesConfig::init_corries")?;
        u.update_vars_from_prim(&mut solver.rhs.boundary_west, &mut solver.rhs.boundary_east, &mesh);
        u.init_west_east();
        Ok((u, solver, mesh, writer))
    }
}

impl Validation for CorriesConfig {
    fn validate(&self) -> Result<()> {
        self.mesh_config.validate().context("Validating config.meshconfig")?;
        self.physics_config
            .validate()
            .context("Validating config.physicsconfig")?;
        self.numerics_config
            .validate()
            .context("Validating config.numericsconfig")?;
        for outputconf in self.writer_config.iter() {
            outputconf.validate().context("Validating config.writerconf")?;
        }
        Ok(())
    }
}
