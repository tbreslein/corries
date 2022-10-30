# corries

CORrosive RIEmann Solver (corrosive because it's written in Rust...)

## Ideas

### Threads

Deps:

- Add rayon for parallel iterators ++ rayon feature for ndarray
- Add crossbeam for scoped threads

Things I can probably parallelise:

- The different outputs running their pipelines starting with converting data_matrix to stream_strings
- computing numflux and sources in parallel (and those sources in parallel to one another)
- computations that iterate over the set of equations

Things to read:

- [rust-lang-nursery](https://rust-lang-nursery.github.io/rust-cookbook/concurrency/threads.html)
- [rayon](https://docs.rs/rayon/latest/rayon/)
- [crossbeam](https://docs.rs/crossbeam/latest/crossbeam/)

## TODO

### Small things

- [ ] Fix Mesh unit tests
- [ ] add approx features for ndarray?
- [ ] Add unit tests for Units and Variables
- [ ] more linting?
- [ ] can I macro metadata_dump?

### Big Steps

- [ ] Implement Physics and Euler1DIsot
- [ ] Implement Numflux and HLL
- [ ] Implement TimeInteg and RKF
- [ ] Add boundary conditions
- [ ] Add general initial conditions structure + Noh integration test
- [ ] Add Euler1DAdiabatic + Sod test
- [ ] Add KT solver
- [ ] Add Euler2DIsot + Euler2DAdiabatic
- [ ] Add NavStoIsot + NavStoAdiabatic
- [ ] Play around with threading + add benchmarking pipeline
- [ ] Add cylindrical geometry + Source + GeometricSource + Sedov test
- [ ] Add spherical geometry + new Sedov case
- [ ] Add logcylindrical geometry + Vortex test
- [ ] Add gravitational source + Bondi test
- [ ] Add viscosity source + Pringle test
- [ ] Add Henyey solver
- [ ] big docs
