// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

#![warn(missing_docs)]

//! Corries - CORrosive RIEman Solver
//!
//! A library/framework (honestly what's the difference nowadays anyways) to run 1D-hydrodynamics
//! simulations solved with Riemann solvers.
//! As such, corries is grid-based (in opposition to doing smooth particle hydrodynamics).
//!
//! Corries should be ready to use for simple shocktube-like simulations at this point.
//!
//! # Disclaimer about nomenclature
//!
//! Just a couple of notes about how things are named:
//!
//! * corries uses generalised coordinates, and names them `xi`, `eta`, and `Phi`
//! * the coordinate systems in corries are always right-handed
//! * corries may be specialised for 1D simulations, but the other dimensions are still important,
//! even though they are only 1 cell wide
//! * `xi` is the primary coordinate, the one that has more than one cell. That's why the
//! coordinate vectors in [Mesh] are called xi_*, for example.
//!
//! # Dependencies
//!
//! Corries has one dependency that you need to use in your apps, which is `color_eyre`.
//! Thus, in order to use corries, you need add to lines to the dependencies in your `Cargo.toml`:
//!
//! ```toml
//! [dependencies]
//! color-eyre = { version = "0.6", default-features = false }
//! corries = { version = "0.4" }
//! ```
//!
//! corries also has one default feature: `validation`.
//! This feature performs run-time checks on the state of your simulation, like checking for finite
//! values, positive mass densities and pressures, and the like.
//!
//! As useful as they are, they do slow down the simulation.
//! Such you can disable it by turning off the the default-features by modifying the corries entry
//! like this:
//!
//! ```toml
//! corries = { versoin = "0.4", default-features = false }
//! ```
//!
//! # Usage
//!
//! The main building blocks of preparing and running a simulation with corries are:
//!
//! * Choosing a couple of type parameters, namely for the traits: [Physics], [NumFlux], and
//! [TimeSolver], as well as the [Mesh] size (also through a compile time constant)
//! * Building a [CorriesConfig], which allows for fine-grained control over the simulation
//! * Constructing and applying initial conditions to the necessary objects through the
//! [CorriesConfig::init_corries]
//! * Running the simulation with the [CorriesComponents]'s [Runner::run_corries] method
//!
//! You can see examples of this in the docs for [CorriesConfig::init_corries] and
//! [Runner::run_corries], as well as when looking through the source code for the integration
//! tests in the `tests/` folder.
//! In those examples I usually use a couple of default constructors for configuration structs, but
//! in `examples/template.rs` you can also see how a full-blown [CorriesConfig] can be set up.
//!
//! # Building a Sod test
//!
//! Let's walk through building a Sod test.
//!
//! The first thing we should do is setting some compile constants and type definitions.
//! The type defs are not necessary, but I highly recommend adding them to make your code more
//! flexible when you decide to switch out some of these type defs.
//!
//! ```
//! use color_eyre::Result;
//! use corries::prelude::*;
//!
//! // The number of grid cells in the Mesh
//! const S: usize = 100;
//!
//! // Expands to:
//! // type P = Euler1DAdiabatic<S>;
//! // const E: usize = P::NUM_EQ; // which in turn equals 3
//! set_Physics_and_E!(Euler1DAdiabatic);
//!
//! // The type for the numerical flux scheme
//! type N = Kt<E, S>;
//!
//! // The type for the time discretisation scheme
//! type T = RungeKuttaFehlberg<P, E, S>;
//! ```
//!
//! To break this down, first we have a couple of `use` statements.
//! corries features a prelude module that imports everything you need, and then you also need to
//! use `color_eyre::Result`, because that is the return type for a lot of functions in `corries`.
//!
//! Then we need to set the size of your mesh, by definining `const S: usize`.
//! In this case, we set the mesh to have 100 grid cells (including ghost cells).
//!
//! Next, we set up our [Physics] type.
//! This determines the differential equations that are being solved in the simulation.
//! You can check out the current options in the `Implementors` section of [Physics].
//! This macro also sets up `const E: usize` to be the number of equations for the [Physics]
//! implementor you chose.
//!
//! Then we set `N` to the [NumFlux] implementor we want, and `T` to the [TimeSolver] implementor.
//! Again, check the docs for those traits to see the options.
//!
//! Next up, we need to set up our [CorriesConfig].
//! This is probably the most complicated part, because there are so many options.
//! In case this piece of the docs ever gets outdated, I would implore you to the [CorriesConfig]
//! docs to check out the options you have, as well as the default constructors for that struct as
//! well as some fields in it.
//!
//! Either way, a hand-rolled [CorriesConfig] for this type of simulation would look like this:
//!
//! ```
//! # use color_eyre::Result;
//! # use corries::prelude::*;
//! # const S: usize = 100;
//! # set_Physics_and_E!(Euler1DAdiabatic);
//! # type N = Kt<E, S>;
//! # type T = RungeKuttaFehlberg<P, E, S>;
//! let config: CorriesConfig = CorriesConfig {
//!     // prints a welcome banner upon running [init_corries]
//!     print_banner: true,
//!
//!     // Sets up a cartesian mesh, with the edges of the computational area (not counting ghost
//!     // cells) are at the coordinates `xi = 1.0` and `xi = 2.0`. Remember that `xi` is the name
//!     // of the primary coordinate in corries.
//!     mesh_config: MeshConfig {
//!         mode: MeshMode::Cartesian,
//!         xi_in: 1.0,
//!         xi_out: 2.0,
//!     },
//!
//!     // Since we want to run an adiabatic simulaion, we need to set the adiabatic index.
//!     // For the Sod test, let's set it to `1.4`.
//!     // We can also set the type of units we want, but there is not much reason for it in
//!     // corries right now (yes, I should just remove it for now...). Just set it to
//!     // [UnitsMode::SI].
//!     physics_config: PhysicsConfig {
//!         adiabatic_index: 1.4,
//!         units_mode: UnitsMode::SI,
//!     },
//!
//!     // Sets the boundary conditions for the west (inner) edge of the computational area.
//!     // [BoundaryMode::Custom] allows you to set individual boundary conditions for
//!     // every equation by putting pairs of the equation index the the custom boundary
//!     // mode in the embedded vector. Check out the documentation for the [Physics]
//!     // implementor you chose to see which indexes correspond to which equations, and check out
//!     // the documentation for [CustomBoundaryMode] to see what the available conditions are.
//!     //
//!     // Alternatively, this field could have been set to:
//!     // `boundary_condition_west: BoundaryMode::NoGradients`
//!     // This would set the `NoGradients` condition for all equations, without setting up this
//!     // vector. See `boundary_condition_east` to see how that would look like.
//!     boundary_condition_west: BoundaryMode::Custom(vec![
//!         (0, CustomBoundaryMode::NoGradients),
//!         (1, CustomBoundaryMode::NoGradients),
//!         (2, CustomBoundaryMode::NoGradients),
//!     ]),
//!
//!     // Sets the boundary conditions for the east (outer) edge of the computational area.
//!     boundary_condition_east: BoundaryMode::NoGradients,
//!
//!     // Sets up everything regarding numerics.
//!     numerics_config: NumericsConfig {
//!
//!         // Sets up configuration for the numerical flux scheme.
//!         // [NumFluxConfig] is an enum that needs to correspond to the [NumFlux] implementor you
//!         // chose up top, otherwise the constructor for the [NumFlux] object will panic.
//!         //
//!         // `Kt` should be your go-to. It's fast and accurate, without running into Godunov's
//!         // order barrier since we only use linear limiting functions.
//!         numflux_config: NumFluxConfig::Kt {
//!             // When choosing `Kt`, you also need to set the limiter function here. Your go-to
//!             // should be the `VanLeer` limiter.
//!             limiter_mode: LimiterMode::VanLeer,
//!         },
//!
//!         // Sets up the time integration scheme.
//!         // Currently, corries only supports Runge-Kutta-Fehlberg schemes, which are set here.
//!         time_integration_config: TimeIntegrationConfig::Rkf(RkfConfig {
//!             // The only important bit about this config is the exact scheme you want to use.
//!             // SSPRK5 and RKF4 are the go-to choices, though I would recommend the first
//!             // usually.
//!             rkf_mode: RKFMode::SSPRK5,
//!
//!             // Whether or not to use automatic step control. Not really needed for this test,
//!             // so we set this to false.
//!             asc: false,
//!
//!             // The relative tolerance when calculating the error between high and low order
//!             // solutions, used by the automatic step control.
//!             asc_relative_tolerance: 0.001,
//!
//!             // The absolute tolerance when calculating the error between high and low order
//!             // solutions, used by the automatic step control.
//!             asc_absolute_tolerance: 0.1,
//!
//!             // The timestep friction factor. This dampens new time steps widths calculated by
//!             // the automatic step control.
//!             asc_timestep_friction: 0.08,
//!         }),
//!
//!         // How many full loop iterations corries is allowed to run before aborting the
//!         // simulation, if it has not ended due to other break conditions yet.
//!         iter_max: usize::MAX - 2,
//!
//!         // The time coordinate at which the simulation starts.
//!         t0: 0.0,
//!
//!         // The time coordinate at which the simulation ends. We set this to 0.25, so that the
//!         // the simulation stays contained within our spatial bounds.
//!         t_end: 0.25,
//!
//!         // If the calculated time step width falls below this value, the simulation ends with
//!         // an error. This is to prevent simulations that would run on too long because
//!         // something is driving the time step width into the ground.
//!         dt_min: 1.0e-12,
//!
//!         // The maximum value for the time step width. If the calculated time step width exceeds
//!         // this value, former will simply be capped to this.
//!         dt_max: f64::MAX,
//!
//!         /// The CFL (Courant-Friedrichs-Lewy) condition parameter. Basically a safety factor to
//!         // dampen the time step widths calculated by the CFL condition.
//!         dt_cfl_param: 0.4,
//!     },
//!
//!     // How many times should [Writer] write outputs during the simulation (not counting the
//!     // initial output).
//!     // These are evenly distributed throughout the simulation time, so for example, if you set
//!     // `numerics_config.t0 = 0.0` and `numerics_config.t_end = 10.0`, then setting
//!     // `output_counter_max = 10` would mean that corries writes output at
//!     // `t = [1.0, 2.0, ..., 10.0]` in addition to the extra output at `t = 0.0`.
//!     output_counter_max: 10,
//!
//!     // Sets up how outputs are generated.
//!     // This is a vector, where every entry is a instance of [OutputConfig], and each of these
//!     // instances sets up one output stream. In this example, we set up to streams: one for
//!     // stdout (the first entry), and one for files (the second entry). You are free to mix and
//!     // match different output streams as you like, but this is probably a good go-to.
//!     writer_config: vec![
//!
//!         // Sets up a stdout output stream, with TSV formatting (tab separated values) for
//!         // scalar values, i.e. values that do not differ for different parts of the
//!         // computational area.
//!         //
//!         // In most cases, setting this stream up like this is probably too verbose, so
//!         // OutputConfig provides constructor that you only need to pass the data_names field
//!         // to.
//!         // This same setup could have been achieved by replacing this initialiser with:
//!         // OutputConfig::default_stdout_with_names(
//!         //     vec![DataName::Iter, DataName::T, DataName::Dt, DataName::DtKind]
//!         // ),
//!         //
//!         // You can also replicate this setup with the default_stout constructor:
//!         //
//!         // OutputConfig::default_stdout(),
//!         //
//!         // that sets up data_names exactly like this too.
//!         OutputConfig {
//!             // Sets up that this stream writes to stdout
//!             stream_mode: StreamMode::Stdout,
//!
//!             // Formats the output vales such that the values are separated by tabs
//!             formatting_mode: FormattingMode::TSV,
//!
//!             // Tells the stream to only print "global" data, i.e. values that are the same
//!             // everywhere in the computational area.
//!             string_conversion_mode: ToStringConversionMode::Scalar,
//!
//!             // Ignored for stdout output; for file output this field sets the folder the
//!             // simulation data should be written to.
//!             folder_name: "".to_string(),
//!
//!             // Ignored for stdout output; Whether to clear the contents of that folder during
//!             // initialisation.
//!             // THIS CAN LEAD TO LOSS OF DATA SO BE VERY CAREFUL WHAT FOLDER_NAME POINTS TO WHEN
//!             // SETTING THIS VALUE TO true!
//!             should_clear_out_folder: true,
//!
//!             // Irrelevant for stdout output; for file output this field sets the base file name
//!             // for the simulation data.
//!             file_name: "".to_string(),
//!
//!             // Floating point precision when printing floating point numbers
//!             precision: 3,
//!
//!             // Ignored for scalar output; for vector output, this sets whether to include the
//!             // ghost cells in the output.
//!             should_print_ghostcells: false,
//!
//!             // Whether to print the configuration metadata. For file based output, this
//!             // generates a <file_name>__metadata.json next to your simulation data, which is the
//!             // deserialised [CorriesConfig] for this simulation.
//!             should_print_metadata: true,
//!
//!             // Defines which values should be printed to the output.
//!             // These values will be printed from left to right in the output. Check out the
//!             // docs for [DataName] for the options you have.
//!             data_names: vec![DataName::Iter, DataName::T, DataName::Dt, DataName::DtKind],
//!         },
//!
//!         // Sets up a file output stream.
//!         //
//!         // There is also a handy constructor that sets up a default file output.
//!         // You could also define this exact same setup by constructing the OutputConfig like
//!         // this:
//!         // OutputConfig::default_file_with_names(
//!         //     "results/sod",
//!         //     "sod",
//!         //     vec![DataName::T, DataName::XiCent, DataName::Prim(0), DataName::Prim(1), DataName::Prim(2)],
//!         // ),
//!         //
//!         // There is also OutputConfig::default_file() that sets the data_names field like this
//!         // automatically. You would call it like this:
//!         //
//!         // OutputConfig::default_file("results/sod", "sod", E),
//!         //
//!         // where E is the number of equations (we set this through the set_Physics_and_E macro
//!         // up top. This function would also replicate this exact setup.
//!         OutputConfig {
//!             // Configures this object to write to a file
//!             stream_mode: StreamMode::File,
//!
//!             // Formats the data into csv files
//!             formatting_mode: FormattingMode::CSV,
//!
//!             // Defines that the data is written as vectors. This allows the output to write
//!             // data like DataName::XiCent or DataName::Prim(n), which are vector-like values.
//!             string_conversion_mode: ToStringConversionMode::Vector,
//!
//!             // Which folder to write the data to, relative to the cwd when calling the
//!             // executable.
//!             folder_name: "results/sod".to_string(),
//!
//!             // Whether to clear out the folder defined in `folder_name` during initialisation.
//!             // THIS CAN LEAD TO LOSS OF DATA SO BE VERY CAREFUL WHAT FOLDER_NAME POINTS TO WHEN
//!             // SETTING THIS VALUE TO true!
//!             should_clear_out_folder: true,
//!
//!             // The base file name for output files.
//!             // In this setup, your output files would be named:
//!             // `results/sod/sod_{00,01,02,..,10}.csv`
//!             file_name: "sod".to_string(),
//!
//!             // Floating point precision for floating point numbers
//!             precision: 7,
//!
//!             // Whether to include ghostcells when writing vector-like values.
//!             should_print_ghostcells: true,
//!
//!             // Whether to write the metadata dump. In this setup, this file would be called:
//!             // `results/sod/sod__metadata.json`
//!             should_print_metadata: false,
//!
//!             // Which values to include in this output, read from left to right.
//!             data_names: vec![
//!                 DataName::T, DataName::XiCent, DataName::Prim(0), DataName::Prim(1), DataName::Prim(2)
//!             ],
//!         },
//!     ],
//! };
//! ```
//!
//! Now with our configuration done, we can initialise the objects we need.
//! The [CorriesConfig::init_corries()] method does exactly what you need here, and all you need to
//! pass it are the generic parameters we set up top, as well as a function for your initial
//! conditions.
//!
//! That function has the signature (barring generic parameters):
//!
//! ```ignore
//! fn(&mut State, &mut Solver, &Mesh) -> Result<()>
//! ```
//!
//! As was already stated, the purpose of it is to apply your initial conditions to the [State]
//! object.
//!
//! The Sod test for this set of equations is set up like this:
//!
//! * net zero velocities
//! * a density and pressure jump in the horizontal centre of the shocktube:
//!   * To the left of that shock the mass density and pressure are set to 1.0
//!   * To the right of that shock the mass density is set to 0.125 and the pressure is set to 1.0
//!
//! When setting up initial conditions we generally only have to set up the primitive variables for
//! the cell centric variable set.
//! The [CorriesConfig::init_corries()] method handles applying the boundary conditions and making
//! sure that [State::west] and [State::east] are set up too.
//!
//! ```no_run
//! # use color_eyre::{eyre::Context, Result};
//! # use corries::prelude::*;
//! # const S: usize = 100;
//! # set_Physics_and_E!(Euler1DAdiabatic);
//! # type N = Kt<E, S>;
//! # type T = RungeKuttaFehlberg<P, E, S>;
//! # fn get_config<N: NumFlux<E, S> + 'static, const E: usize>(folder_name: &str, file_name: &str) -> CorriesConfig {
//! #     return CorriesConfig {
//! #         print_banner: false,
//! #         mesh_config: MeshConfig::default_riemann_test(),
//! #         physics_config: PhysicsConfig {
//! #             units_mode: UnitsMode::SI,
//! #             adiabatic_index: 1.4,
//! #         },
//! #         boundary_condition_west: BoundaryMode::NoGradients,
//! #         boundary_condition_east: BoundaryMode::NoGradients,
//! #         numerics_config: NumericsConfig::default_riemann_test::<N, E, S>(0.25),
//! #         output_counter_max: 1,
//! #         writer_config: vec![
//! #             OutputConfig::default_stdout(),
//! #             OutputConfig::default_file(folder_name, file_name, E),
//! #         ],
//! #     };
//! # }
//! # let config = get_config::<N,E>("results/sod", "sod");
//! config.init_corries::<P, N, T, E, S>(|u, _, _| {
//!     // this function does need the [Solver] and [Mesh] objects, hence the two `_` in the
//!     // arguments.
//!     let breakpoint_index = (S as f64 * 0.5) as usize;
//!     for i in 0..breakpoint_index {
//!         // set mass density (rho) and pressure of the cell central values left of the shock to
//!         // 1.0
//!         u.cent.prim[[P::JRHO, i]] = 1.0;
//!         u.cent.prim[[P::JPRESSURE, i]] = 1.0;
//!
//!         // set the xi velocity to 0.0
//!         u.cent.prim[[P::JXI, i]] = 0.0;
//!     }
//!     for i in breakpoint_index..S {
//!         // set the mass density to the right of the shock to 0.125
//!         u.cent.prim[[P::JRHO, i]] = 0.125;
//!
//!         // set the pressure to the right of the shock to 0.1
//!         u.cent.prim[[P::JPRESSURE, i]] = 0.1;
//!
//!         // set the xi velocity to 0.0
//!         u.cent.prim[[P::JXI, i]] = 0.0;
//!     }
//!     Ok(())
//! }).unwrap();
//! ```
//!
//! Obviously, you should be using the ? operator instead of unwrap, but then this doc test does
//! not run correctly... ¯\\_(ツ)_/¯
//!
//! The object returned here is a [CorriesComponents], which is simply a tuple around:
//!
//! * `u`: the [State] object that holds the state of the physical variables
//! * `solver`: the [Solver] object that generates new solutions
//! * `mesh`: the [Mesh] object that models the grid the simulation runs on
//! * `writer`: the [Writer] object that handles writing output
//!
//! Almost done.
//! Now we can call [Runner::run_corries] to run the simulation!
//!
//! ```no_run
//! # use color_eyre::{eyre::Context, Result};
//! # use corries::prelude::*;
//! # const S: usize = 100;
//! # set_Physics_and_E!(Euler1DAdiabatic);
//! # type N = Kt<E, S>;
//! # type T = RungeKuttaFehlberg<P, E, S>;
//! # fn get_config<N: NumFlux<E, S> + 'static, const E: usize>(folder_name: &str, file_name: &str) -> CorriesConfig {
//! #     return CorriesConfig {
//! #         print_banner: false,
//! #         mesh_config: MeshConfig::default_riemann_test(),
//! #         physics_config: PhysicsConfig {
//! #             units_mode: UnitsMode::SI,
//! #             adiabatic_index: 1.4,
//! #         },
//! #         boundary_condition_west: BoundaryMode::NoGradients,
//! #         boundary_condition_east: BoundaryMode::NoGradients,
//! #         numerics_config: NumericsConfig::default_riemann_test::<N, E, S>(0.25),
//! #         output_counter_max: 1,
//! #         writer_config: vec![
//! #             OutputConfig::default_stdout(),
//! #             OutputConfig::default_file(folder_name, file_name, E),
//! #         ],
//! #     };
//! # }
//! # let config = get_config::<N,E>("results/sod", "sod");
//! config.init_corries::<P, N, T, E, S>(|u, _, _| {
//!     // this function does need the [Solver] and [Mesh] objects, hence the two `_` in the
//!     // arguments.
//!     let breakpoint_index = (S as f64 * 0.5) as usize;
//!     for i in 0..breakpoint_index {
//!         // set mass density (rho) and pressure of the cell central values left of the shock to
//!         // 1.0
//!         u.cent.prim[[P::JRHO, i]] = 1.0;
//!         u.cent.prim[[P::JPRESSURE, i]] = 1.0;
//!
//!         // set the xi velocity to 0.0
//!         u.cent.prim[[P::JXI, i]] = 0.0;
//!     }
//!     for i in breakpoint_index..S {
//!         // set the mass density to the right of the shock to 0.125
//!         u.cent.prim[[P::JRHO, i]] = 0.125;
//!
//!         // set the pressure to the right of the shock to 0.1
//!         u.cent.prim[[P::JPRESSURE, i]] = 0.1;
//!
//!         // set the xi velocity to 0.0
//!         u.cent.prim[[P::JXI, i]] = 0.0;
//!     }
//!     Ok(())
//! }).unwrap().run_corries().unwrap();
//! ```
//!
//! # Plans
//!
//! Currently, corries only supports cartesian meshes, though non-cartesian meshes are coming very
//! soon.
//!
//! * Add cylindrical geometry + Source trait + GeometricSource + Sedov test
//! * Add spherical geometry + new Sedov case
//! * Add 2D adiabatic Euler and 2D isothermal Euler
//! * Add adiabatic and isothermal Navier Stokes physics
//! * Add logcylindrical geometry + Vortex test
//! * Add gravitational source + Bondi test
//! * Add viscosity source + Pringle test
//! * investigate the possibility of putting boundary conditions into [State]
//!   * not really possible currently, because I cannot split the borrows between the boundary
//!   conditions and the variables without using Cell for interior mutability.
//!
//! Most importantly exports the module [prelude].

mod boundaryconditions;
pub mod components;
pub mod config;
pub mod directions;
mod errorhandling;
#[macro_use]
pub mod macros;
pub mod initfuncs;
pub mod mesh;
pub mod rhs;
pub mod solver;
pub mod state;
pub mod time;
pub mod units;
pub mod writer;

/// Exports everything you need to run a corries simulation. This includes the following modules
///
/// * [corries::components](crate::components)
/// * [corries::config](crate::config)
/// * [corries::directions](crate::directions)
/// * [corries::initfuncs](crate::initfuncs)
/// * [corries::mesh](crate::mesh)
/// * [corries::rhs](crate::rhs)
/// * [corries::solver](crate::solver)
/// * [corries::state](crate::state)
/// * [corries::time](crate::time)
/// * [corries::units](crate::units)
/// * [corries::writer](crate::writer)
///
/// as well as the [set_Physics_and_E] macro
pub mod prelude {
    pub use crate::components::*;
    pub use crate::config::*;
    pub use crate::directions::*;
    pub use crate::mesh::*;
    pub use crate::rhs::*;
    pub use crate::set_Physics_and_E;
    pub use crate::solver::*;
    pub use crate::state::*;
    pub use crate::time::*;
    pub use crate::units::*;
    pub use crate::writer::*;
}
pub use prelude::*;
