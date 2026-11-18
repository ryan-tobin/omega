//! Runs with forces from the GPU (gpu/forces.rs): they hold energy and
//! momentum, retrace their path, and follow a run on the CPU until
//! rounding differences have had time to grow.
//!
//! Each test starts the 108 atoms of tests/conservation.rs: an fcc lattice
//! at Rahman's density, velocities for his temperature from seed 1, his
//! cutoff and our time step. They need a GPU and the `gpu` feature:
//! cargo test --features gpu. Every tolerance is our own choice.
#![cfg(feature = "gpu")]

use omega::forces::Method;
use omega::gpu::Gpu;
use omega::gpu::forces::GpuForces;
use omega::lattice;
use omega::periodic::PeriodicBox;
use omega::rahman::{CUTOFF, DT, SIDE, TEMPERATURE};
use omega::rng::Rng;
use omega::system::System;
use omega::vec3::Vec3;
use omega::velocities;
use std::sync::OnceLock;

/// Forces computed on this machine's GPU, shared by all tests.
fn on_gpu() -> Method {
    static GPU: OnceLock<Gpu> = OnceLock::new();
    let gpu = GPU.get_or_init(|| Gpu::new().expect("these tests need a GPU"));
    Method::Gpu(Box::new(GpuForces::new(gpu)))
}

/// 108 atoms on an fcc lattice at Rahman's density, with velocities for
/// his temperature, and forces computed by `method`.
fn argon(method: Method) -> System {
    let space = PeriodicBox { side: SIDE / 2.0 };
    let positions = lattice::fcc(3, space.side);
    let mut rng = Rng::new(1);
    let velocities = velocities::maxwell_boltzmann(positions.len(), TEMPERATURE, &mut rng);
    System::with_forces(positions, velocities, space, CUTOFF, method)
}

/// The farthest any entry of `a` is from the same entry of `b`.
fn apart(a: &[Vec3], b: &[Vec3]) -> f64 {
    let pairs = a.iter().zip(b);
    pairs.map(|(&p, &q)| (p - q).norm()).fold(0.0, f64::max)
}

#[test]
fn gpu_run_holds_energy_and_momentum_over_10000_steps() {
    // The 108 atoms and seed of tests/conservation.rs, where all pairs
    // gives a largest energy change of 2.29e-4 and a largest momentum of
    // 7.6e-14. Measured on a software GPU: 2.24e-4 and 7.1e-6. The energy
    // line is the same 5e-4. The momentum line is 1e-4 where the CPU has
    // 1e-10: each atom's force is rounded on its own, so the forces no
    // longer cancel exactly. Without exactly opposite pair forces the
    // momentum reached 4e-4 (measured, see forces.wgsl).
    let mut system = argon(on_gpu());
    let e0 = system.total_energy();
    let (mut worst, mut drift): (f64, f64) = (0.0, 0.0);
    for _ in 0..10_000 {
        system.step(DT);
        worst = worst.max((system.total_energy() - e0).abs() / e0.abs());
        drift = drift.max(system.momentum().norm());
    }
    assert!(worst <= 5e-4, "largest relative change {worst:e}");
    assert!(drift <= 1e-4, "largest momentum {drift:e}");
}

#[test]
fn gpu_run_retraces_its_path() {
    // 1000 steps forward, flip every velocity, 1000 steps back, as in
    // tests/conservation.rs. Measured on a software GPU over 80 seeds: 77
    // runs came back within 9e-15 sigma, closer than the CPU's 1.5e-10.
    // The forces see each position only to the nearest tick (sorted.rs),
    // so they repeat bit for bit on the way back. In the other 3 runs a
    // coordinate fell on the other side of a tick, a force changed by one
    // rounding, and the miss grew to at most 4.0e-7. Line: 1e-5, 25 times
    // that. With a skin of 0.3 sigma on the cells, which makes the forces
    // depend on when the atoms were last sorted, the miss was 1e-2.
    let mut system = argon(on_gpu());
    let start = system.positions.clone();
    for _ in 0..1000 {
        system.step(DT);
    }
    for v in &mut system.velocities {
        *v = -*v;
    }
    for _ in 0..1000 {
        system.step(DT);
    }
    let miss = apart(&system.positions, &start);
    assert!(miss <= 1e-5, "{miss:e} sigma from the start");
}

#[test]
fn gpu_and_all_pairs_runs_stay_together_for_100_steps() {
    // The same start on the GPU's forces and on all pairs. Rounding
    // differences grow as the runs go on. Measured on a software GPU:
    // 1.7e-7 sigma apart after 100 steps, 2.4e-5 after 200, 4.3e-4 after
    // 500 and 0.2 after 1000. Line at 100 steps: 1e-4.
    let mut all = argon(Method::AllPairs);
    let mut gpu = argon(on_gpu());
    for _ in 0..100 {
        all.step(DT);
        gpu.step(DT);
    }
    let gap = apart(&all.positions, &gpu.positions);
    assert!(gap <= 1e-4, "{gap:e} sigma apart");
}
