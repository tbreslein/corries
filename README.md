# corries

CORrosive RIEmann Solver (corrosive because it's written in Rust...)

## TODO

### Small things

- [ ] add `prelude.rs` to define common API exports
- [ ] can I macro metadata_dump?

### Big Steps

- [ ] fix doc string for Output::width
- [ ] Can I zip over data_names and data_matrix in Output::update_data_matrix?
- [ ] turn Physics into a trait and make the functions generic
- [ ] same with Numflux
- [ ] same with TimeSolver
- [ ] Add Sod test
- [ ] Add KT solver
- [ ] Add Euler2DAdiabatic
- [ ] Add NavStoIsot + NavStoAdiabatic
- [ ] Play around with threading + add benchmarking pipeline
- [ ] Add cylindrical geometry + Source + GeometricSource + Sedov test
- [ ] Add spherical geometry + new Sedov case
- [ ] Add logcylindrical geometry + Vortex test
- [ ] Add gravitational source + Bondi test
- [ ] Add viscosity source + Pringle test
- [ ] Add Henyey solver
- [ ] big docs
