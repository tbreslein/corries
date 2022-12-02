// Copyright (c) 2022
// Author: Tommy Breslein (github.com/tbreslein)
// License: MIT

use corries::prelude::*;
use criterion::{criterion_group, criterion_main, Criterion};

const S: usize = 500;
const PHYSICS_CONFIG: PhysicsConfig = PhysicsConfig { adiabatic_index: 1.4 };

pub fn euler1d_isot_conversions(c: &mut Criterion) {
    let mut group = c.benchmark_group("physics_conversions");
    group.sample_size(1000);

    let mut u = Euler1DIsot::<S>::new(&PHYSICS_CONFIG);
    u.prim.row_mut(0).fill(2.0);
    u.prim.row_mut(1).fill(3.0);
    group.bench_function("from prim to cons and back; euler 1d isot", |b| {
        b.iter(|| {
            u.update_cons();
            u.update_prim();
        })
    });

    let mut u = Euler1DAdiabatic::<S>::new(&PHYSICS_CONFIG);
    u.prim.row_mut(0).fill(2.0);
    u.prim.row_mut(1).fill(3.0);
    u.prim.row_mut(2).fill(4.0);
    group.bench_function("from prim to cons and back; euler 1d adiabatic", |b| {
        b.iter(|| {
            u.update_cons();
            u.update_prim();
        })
    });

    group.finish();
}

criterion_group!(benches, euler1d_isot_conversions);
criterion_main!(benches);
