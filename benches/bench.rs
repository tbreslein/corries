// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use corries::prelude::*;
use criterion::{criterion_group, criterion_main, Criterion};
use ndarray::{Array1, Array2};

const S: usize = 500;

fn init<P: Physics<E, S>, N: NumFlux<E, S>, T: TimeSolver<P, E, S>, const E: usize, const S: usize>(
    u: &mut State<P, E, S>,
    solver: &mut Solver<P, N, T, E, S>,
    mesh: &Mesh<S>,
) {
    let breakpoint_index = (S as f64 * 0.5) as usize;
    let mut prim = Array2::ones((E, S));
    for i in breakpoint_index..S {
        prim[[1, i]] = -1.0;
    }
    if u.is_adiabatic() {
        prim.row_mut(E - 1).fill(1.0E-5);
    } else {
        u.cent.c_sound.assign(&Array1::ones(S).view());
    }
    u.cent.prim.assign(&prim.view());
    u.update_vars_from_prim(&mut solver.rhs.boundary_west, &mut solver.rhs.boundary_east, &mesh);
    u.init_west_east();
    return;
}

pub fn noh_hll_run(c: &mut Criterion) {
    let mut group = c.benchmark_group("noh_hll_run");

    set_Physics_and_E!(Euler1DAdiabatic);
    type N = Hll<E, S>;
    type T = RungeKuttaFehlberg<P, E, S>;
    let config = CorriesConfig::default_riemann_test::<N, E, S>(0.5, "results/bench/bench_noh", "bench_noh");
    let (mut u, mut solver, mesh, mut writer) = init_corries::<P, N, T, E, S>(&config).unwrap();
    init(&mut u, &mut solver, &mesh);

    group.bench_function("noh_hll", |b| {
        b.iter(|| {
            run_corries(&mut u, &mut solver, &mesh, &mut writer).unwrap();
        })
    });
    group.finish();
}

pub fn noh_kt_run(c: &mut Criterion) {
    let mut group = c.benchmark_group("noh_kt_run");

    set_Physics_and_E!(Euler1DAdiabatic);
    type N = Kt<E, S>;
    type T = RungeKuttaFehlberg<P, E, S>;
    let config = CorriesConfig::default_riemann_test::<N, E, S>(0.5, "results/bench/bench_noh", "bench_noh");
    let (mut u, mut solver, mesh, mut writer) = init_corries::<P, N, T, E, S>(&config).unwrap();
    init(&mut u, &mut solver, &mesh);

    group.bench_function("noh_kt", |b| {
        b.iter(|| {
            run_corries(&mut u, &mut solver, &mesh, &mut writer).unwrap();
        })
    });
    group.finish();
}

pub fn euler1d_isot_conversions(c: &mut Criterion) {
    const PHYSICS_CONFIG: PhysicsConfig = PhysicsConfig {
        adiabatic_index: 1.4,
        units_mode: UnitsMode::SI,
    };
    let mut group = c.benchmark_group("physics_conversions");
    group.sample_size(1000);

    let mut u = State::<Euler1DIsot<S>, 2, S>::new(&PHYSICS_CONFIG);
    u.cent.prim.row_mut(0).fill(2.0);
    u.cent.prim.row_mut(1).fill(3.0);
    group.bench_function("from prim to cons and back; euler 1d isot", |b| {
        b.iter(|| {
            u.update_cons();
            u.update_prim();
        })
    });

    let mut u = State::<Euler1DAdiabatic<S>, 3, S>::new(&PHYSICS_CONFIG);
    u.cent.prim.row_mut(0).fill(2.0);
    u.cent.prim.row_mut(1).fill(3.0);
    u.cent.prim.row_mut(2).fill(4.0);
    group.bench_function("from prim to cons and back; euler 1d adiabatic", |b| {
        b.iter(|| {
            u.update_cons();
            u.update_prim();
        })
    });

    group.finish();
}

criterion_group!(benches, euler1d_isot_conversions, noh_hll_run, noh_kt_run);
criterion_main!(benches);
