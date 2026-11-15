//! Conservation tests for many atoms: energy, momentum and time reversal.
//!
//! Each test starts 108 atoms (4 * 3^3) on an fcc lattice at Rahman's
//! density, in a box of half his side, with velocities for his
//! temperature. The cutoff is his, 2.25 sigma, with the shifted potential.
//! The time step is 0.005 tau. Every tolerance is our own choice, set from
//! measured values with a margin, and each one says why.

use omega::lattice;
use omega::periodic::PeriodicBox;
use omega::rng::Rng;
use omega::system::System;
use omega::velocities;

const DT: f64 = 0.005;

fn argon_108(seed: u64) -> System {
    let space = PeriodicBox { side: 10.229 / 2.0 };
    let positions = lattice::fcc(3, space.side);
    let mut rng = Rng::new(seed);
    let velocities = velocities::maxwell_boltzmann(positions.len(), 0.7867, &mut rng);
    System::new(positions, velocities, space, 2.25)
}

#[test]
fn total_energy_holds_over_10000_steps() {
    // The total energy is not held exactly. Velocity Verlet makes it wobble
    // by an amount that grows as dt^2. On top of that, the force jumps by
    // 0.081 eps / sigma at the cutoff, so each pair that crosses it gets a
    // small kick of random sign. The kicks add up like a random walk, so
    // the energy drifts slowly. At dt = 0.005 the two are about the same
    // size. Relative to the total, both shrink as 1/sqrt(N), so 108 atoms
    // show about 2.8 times more than 864. Over four seeds the largest
    // change was 0.9e-4 to 2.4e-4 of the starting energy. This seed gives
    // 2.29e-4. Tolerance: 5e-4.
    let mut system = argon_108(1);
    let e0 = system.total_energy();
    let mut worst: f64 = 0.0;
    for _ in 0..10_000 {
        system.step(DT);
        worst = worst.max((system.total_energy() - e0).abs() / e0.abs());
    }
    assert!(worst <= 5e-4, "largest relative change {worst:e}");
}

#[test]
fn total_momentum_stays_zero_over_10000_steps() {
    // Every pair pushes its two atoms equally and oppositely, so the total
    // momentum cannot change. What is left is rounding. Measured largest
    // value: 7.6e-14. Tolerance: 1e-10.
    let mut system = argon_108(1);
    let mut worst: f64 = system.momentum().norm();
    for _ in 0..10_000 {
        system.step(DT);
        worst = worst.max(system.momentum().norm());
    }
    assert!(worst <= 1e-10, "largest momentum {worst:e}");
}

#[test]
fn running_forward_then_backward_returns_to_the_start() {
    // Velocity Verlet is time-reversible: flip every velocity and the same
    // steps retrace the path. Rounding leaves errors near 1e-16 at each
    // step, and the motion is chaotic, so they grow fast. Measured return
    // miss for positions: 3.4e-13 sigma after 500 steps each way, 1.5e-10
    // after 1000 (5 tau), 4.6e-4 after 2000. For velocities: 3.3e-12,
    // 5.4e-10 and 2.0e-3 sigma / tau. Tolerance after 1000 steps: 1e-8
    // sigma for positions and 1e-8 sigma / tau for velocities.
    let mut system = argon_108(1);
    let start_positions = system.positions.clone();
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
    for (i, r) in system.positions.iter().enumerate() {
        let miss = (*r - start_positions[i]).norm();
        assert!(miss <= 1e-8, "atom {i} is {miss:e} sigma from its start");
    }
    for (i, v) in system.velocities.iter().enumerate() {
        let miss = (*v + start_velocities[i]).norm();
        assert!(miss <= 1e-8, "atom {i} velocity off by {miss:e}");
    }
}
