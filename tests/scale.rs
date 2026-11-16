//! Forces from neighbor lists and threads in running systems (forces.rs):
//! they match all pairs at every step, hold energy and momentum, retrace
//! their path, give the same run for any number of threads, and follow
//! the all-pairs run until rounding errors have had time to grow.
//!
//! Each test starts an fcc lattice at Rahman's density, with velocities
//! for his temperature, his cutoff and a skin of 0.3 sigma. Every
//! tolerance is our own choice and says why.

use omega::forces::Method;
use omega::lattice;
use omega::pairs;
use omega::periodic::PeriodicBox;
use omega::rahman::{CUTOFF, DT, SIDE, TEMPERATURE};
use omega::rng::Rng;
use omega::system::System;
use omega::vec3::Vec3;
use omega::velocities;

const SKIN: f64 = 0.3;

/// 4 * cells^3 atoms on an fcc lattice at Rahman's density, with
/// velocities for his temperature drawn from `seed`.
fn argon(cells: usize, seed: u64, method: Method) -> System {
    let space = PeriodicBox {
        side: SIDE / 6.0 * cells as f64,
    };
    let positions = lattice::fcc(cells, space.side);
    let mut rng = Rng::new(seed);
    let velocities = velocities::maxwell_boltzmann(positions.len(), TEMPERATURE, &mut rng);
    System::with_forces(positions, velocities, space, CUTOFF, method)
}

/// The farthest any atom of `a` is from the same atom of `b`, in sigma.
fn apart(a: &[Vec3], b: &[Vec3]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(&p, &q)| (p - q).norm())
        .fold(0.0, f64::max)
}

#[test]
fn listed_forces_match_all_pairs_at_every_step() {
    // 500 atoms, 3 cells a side. The crystal does not melt (T falls to
    // about 0.43), but its atoms move far enough that the list is built
    // again every 12 steps or so. Every 10 steps the forces, energy and
    // virial must match all pairs at the same positions. Measured over
    // 400 steps: forces within 2.2e-14 eps / sigma, on forces up to 74,
    // energy and virial within 7e-15 relative, and 34 builds. Tolerances
    // 1e-12, as in tests/forces.rs.
    let mut system = argon(5, 1964, Method::listed(SKIN, 2));
    let mut forces = vec![Vec3::new(0.0, 0.0, 0.0); 500];
    for step in 1..=400 {
        system.step(DT);
        if step % 10 != 0 {
            continue;
        }
        let (u, w) = pairs::compute(&system.positions, system.space, CUTOFF, &mut forces);
        assert!(apart(&system.forces, &forces) <= 1e-12, "step {step}");
        assert!((system.potential / u - 1.0).abs() <= 1e-12, "step {step}");
        assert!((system.virial / w - 1.0).abs() <= 1e-12, "step {step}");
    }
    assert!(
        system.method.builds() >= 10,
        "{} builds",
        system.method.builds()
    );
}

#[test]
fn listed_run_holds_energy_and_momentum_over_10000_steps() {
    // The 108 atoms and seed of tests/conservation.rs, where all pairs
    // gives a largest energy change of 2.29e-4. Measured here: 2.15e-4,
    // and momentum up to 1.4e-13, with the list built 752 times. The pass
    // lines are those of all pairs: 5e-4 and 1e-10.
    let mut system = argon(3, 1, Method::listed(SKIN, 2));
    let e0 = system.total_energy();
    let (mut worst, mut drift): (f64, f64) = (0.0, 0.0);
    for _ in 0..10_000 {
        system.step(DT);
        worst = worst.max((system.total_energy() - e0).abs() / e0.abs());
        drift = drift.max(system.momentum().norm());
    }
    assert!(worst <= 5e-4, "largest relative change {worst:e}");
    assert!(drift <= 1e-10, "largest momentum {drift:e}");
}

#[test]
fn listed_run_retraces_its_path() {
    // As in tests/conservation.rs: 1000 steps forward, flip every
    // velocity, 1000 steps back. Measured miss: 2.1e-10 sigma for the
    // positions and 1.1e-9 sigma / tau for the velocities (all pairs:
    // 1.5e-10 and 5.4e-10). Tolerances 1e-8, the same as for all pairs.
    let mut system = argon(3, 1, Method::listed(SKIN, 2));
    let start = system.positions.clone();
    let start_velocities = system.velocities.clone();
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
    assert!(miss <= 1e-8, "{miss:e} sigma from the start");
    let back: Vec<Vec3> = system.velocities.iter().map(|&v| -v).collect();
    let miss = apart(&back, &start_velocities);
    assert!(miss <= 1e-8, "velocities off by {miss:e}");
}

#[test]
fn any_number_of_threads_gives_the_same_run() {
    // 500 atoms, 200 steps, on 1 and on 3 threads. Exact.
    let mut one = argon(5, 7, Method::listed(SKIN, 1));
    let mut three = argon(5, 7, Method::listed(SKIN, 3));
    for _ in 0..200 {
        one.step(DT);
        three.step(DT);
    }
    assert!(one.positions == three.positions);
    assert!(one.velocities == three.velocities);
    assert_eq!(one.potential, three.potential);
    assert_eq!(one.method.builds(), three.method.builds());
}

#[test]
fn listed_and_all_pairs_runs_stay_together_for_5_tau() {
    // The two differ only in the order of their sums, by rounding at
    // each step, and chaos makes such differences grow. Measured: 6.5e-11
    // sigma apart after 1000 steps (5 tau), 4.8e-5 after 2000, and 0.57
    // after 3000, by when they have become two separate runs.
    // Tolerance after 1000 steps 1e-8 sigma, our own choice.
    let mut all = argon(3, 1, Method::AllPairs);
    let mut listed = argon(3, 1, Method::listed(SKIN, 2));
    for _ in 0..1000 {
        all.step(DT);
        listed.step(DT);
    }
    let gap = apart(&all.positions, &listed.positions);
    assert!(gap <= 1e-8, "{gap:e} sigma apart");
}