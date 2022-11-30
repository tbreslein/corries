# corries

CORrosive RIEmann Solver (corrosive because it's written in Rust...)

## TODO

### Small things

- [ ] add `prelude.rs` to define common API exports

### Big Steps

- [ ] Get rid of the BoundaryConditionsContainer and just store 2 BoundaryCond instances
- [ ] turn Physics into a trait and make the functions generic
- [ ] same with Numflux
- [ ] same with TimeSolver
- [ ] turn calc_dt_expl into an external generic function and make it return `(double, DtKind)`
- [ ] find more opportunities to use `ensure!`, like in calc_dt_expl
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
