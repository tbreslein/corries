use approx::assert_relative_eq;
use corries::{
    config::physicsconfig::{PhysicsConfig, PhysicsMode},
    physics::{self, init_physics},
};
use criterion::{criterion_group, criterion_main, Criterion};

const S: usize = 10_000;
const PHYSICS_MODE: PhysicsMode = PhysicsMode::Euler1DAdiabatic;
const EQ: usize = corries::get_n_equations(PHYSICS_MODE);

pub fn criterion_benchmark(c: &mut Criterion) {
    let physicsconf = PhysicsConfig {
        adiabatic_index: 0.5,
        c_sound_0: 1.0,
        mode: PHYSICS_MODE,
        units_mode: corries::units::UnitsMode::SI,
    };
    let mut group = c.benchmark_group("Physics implementations");

    let u0: physics::Physics<S, EQ> = init_physics(&physicsconf);
    let mut u: physics::Physics<S, EQ> = init_physics(&physicsconf);
    u.update_cons();
    u.update_prim();
    u.update_cons();
    u.update_prim();
    u.update_cons();
    u.update_prim();
    assert_relative_eq!(u.prim, u0.prim, epsilon = 10.0f64.powi(-12));

    group.sample_size(1000);
    group.bench_function("update", |b| {
        b.iter(|| {
            u.update_cons();
            u.update_prim();
        })
    });
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
