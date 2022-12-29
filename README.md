# corries

CORrosive RIEmann Solver (corrosive because it's written in Rust...) for 1d hydrodynamics problems specialising, planned to be specialised for non-relativistic astrophysics.

## What does this thing do?

This library allows one to write setups for hydrodynamics simulations, and then let those simulations run.
Since this library is supposed to be used for long-running simulations, depending on the setup running for days or weeks, many of the optimisations happen at compile time.
This means that even your parameter studies should be setup such that the parameter study runs from a single executable, either by passing the paramters over the command line, or backing the study's parameter space into the main function.

The things Corries can simulate are fairly limited at the time as it only feature cartesian meshes and no source terms.
You can check the TODO section for what's coming up, if you know what those terms mean.

The plan with this library is to:

- add at least two more numerical flux solvers (a Kurganov-Tadmor solver, as well as a simple finite differences solver)
- add equation systems with additional momentum components (at least 2D)
  - Note that Corries still stays 1D even with those 2D systems, as they assume symmetry along their axes; for example in a disk we would resolve the radius, but assume rotational symmetry
- add non-cartesian symmetry, most notable cylindrical and logcylindrical meshes
- add capabilities for source terms like gravity and viscosity
- add more time integration schemes

## Preparing a simulation

Currently, you can the integration tests as examples for how to set up a simulation, as well as the `examples/template.rs` file.
The base steps are always:

- write a `Corries::Config` struct
- initialise all objects that `Corries::run_corries` expects; you can usually use the `Corries::init_corries` function for that
- run `Corries::run_corries`

As you can see in those examples, the configuration struct is fairly big, and you also need to set up a couple of compile time values that the generic functions need.

### `Corries::Config`

TODO

### Compile time configuration

TODO

## Building

TODO

## Examples ++ Tests

TODO

## TODO

- more detailed README
- Add Sod test
- Add KT solver
- Add Euler2DAdiabatic
- Add NavStoIsot + NavStoAdiabatic
- Add cylindrical geometry + Source + GeometricSource + Sedov test
- Add spherical geometry + new Sedov case
- Add logcylindrical geometry + Vortex test
- Add gravitational source + Bondi test
- Add viscosity source + Pringle test
- big docs
