//! Velocity Verlet on two atoms (system.rs): its energy error shrinks as
//! dt^2, and the momentum is the sum of the velocities. These tests were in
//! system.rs until M5, when the file grew too long to hold them.

use omega::periodic::PeriodicBox;
use omega::system::System;
use omega::vec3::Vec3;

/// Two atoms 1.05 sigma apart, at rest, in a box of side 8. They swing
/// between 1.05 and 1.26 sigma and never reach the cutoff at 2.5, so the
/// force they feel is smooth. Returns the largest change in total energy,
/// in eps, over `steps` steps of `dt`.
fn two_atom_energy_error(dt: f64, steps: usize) -> f64 {
    let positions = vec![Vec3::new(3.475, 4.0, 4.0), Vec3::new(4.525, 4.0, 4.0)];
    let velocities = vec![Vec3::new(0.0, 0.0, 0.0); 2];
    let mut system = System::new(positions, velocities, PeriodicBox { side: 8.0 }, 2.5);
    let e0 = system.total_energy();
    let mut worst: f64 = 0.0;
    for _ in 0..steps {
        system.step(dt);
        worst = worst.max((system.total_energy() - e0).abs());
    }
    worst
}

#[test]
fn energy_error_shrinks_four_times_when_dt_halves() {
    // Velocity Verlet is a second-order method, so its energy error grows
    // as dt^2: halving dt should cut it by a factor of 4. Both runs cover
    // 2 tau, almost three swings of 0.70 tau. Measured: 4.0008. The
    // tolerance, 0.05, is our own choice.
    let coarse = two_atom_energy_error(0.005, 400);
    let fine = two_atom_energy_error(0.0025, 800);
    let ratio = coarse / fine;
    assert!((ratio - 4.0).abs() <= 0.05, "ratio {ratio}");
}

#[test]
fn momentum_is_the_sum_of_the_velocities() {
    // Two atoms 4 sigma apart, beyond the cutoff. Exact sums.
    let positions = vec![Vec3::new(1.0, 1.0, 1.0), Vec3::new(5.0, 1.0, 1.0)];
    let velocities = vec![Vec3::new(1.0, 2.0, 3.0), Vec3::new(0.5, 0.0, -1.0)];
    let system = System::new(positions, velocities, PeriodicBox { side: 8.0 }, 2.5);
    assert_eq!(system.momentum(), Vec3::new(1.5, 2.0, 2.0));
}