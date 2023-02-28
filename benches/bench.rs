// Copyright (c) 2022-2023
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use corries::{initfuncs::init_noh, prelude::*};
use criterion::{criterion_group, criterion_main, Criterion};

const S: usize = 500;

pub fn noh_hll_run(c: &mut Criterion) {
    let mut group = c.benchmark_group("noh_hll_run");

    set_Physics_and_E!(Euler1DAdiabatic);
    type N = Hll<E, S>;
    type T = RungeKuttaFehlberg<P, E, S>;
    let mut comps = CorriesConfig::default_riemann_test::<N, E, S>(0.5, "results/bench/bench_noh", "bench_noh")
        .init_corries::<P, N, T, E, S>(init_noh)
        .unwrap();

    group.bench_function("noh_hll", |b| {
        b.iter(|| {
            comps.run_corries().unwrap();
        })
    });
    group.finish();
}

pub fn noh_kt_run(c: &mut Criterion) {
    let mut group = c.benchmark_group("noh_kt_run");

    set_Physics_and_E!(Euler1DAdiabatic);
    type N = Kt<E, S>;
    type T = RungeKuttaFehlberg<P, E, S>;
    let mut comps = CorriesConfig::default_riemann_test::<N, E, S>(0.5, "results/bench/bench_noh", "bench_noh")
        .init_corries::<P, N, T, E, S>(init_noh)
        .unwrap();

    group.bench_function("noh_kt", |b| {
        b.iter(|| {
            comps.run_corries().unwrap();
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
