//! The thermostat (thermostat.rs) at work on liquid argon: 108 atoms at
//! Rahman's density and temperature (tests/common/mod.rs), with a time
//! step of 0.005 tau. Every tolerance is our own choice, set from measured
//! values with a margin, and each one says why.

mod common;

use common::{T, liquid, step};
use omega::thermo;

#[test]
fn atoms_and_bath_together_keep_their_energy() {
    // The thermostat moves energy in and out of the atoms, so their own
    // energy E wanders. E plus the energy of the bath (thermostat.rs) is
    // conserved by the exact motion, and by our steps as well as E is
    // without a thermostat (tests/conservation.rs). Over 20,000 steps and
    // four seeds its largest change was 3.3e-4 to 1.2e-3 of its start,
    // while E itself moved by 8% to 24%. This seed: 1.2e-3 and 24%.
    // Tolerances: under 3e-3, and E moving by over 2%.
    let (mut system, mut thermostat) = liquid(1);
    let e0 = system.total_energy();
    let h0 = e0 + thermostat.energy(108);
    let (mut worst, mut wander): (f64, f64) = (0.0, 0.0);
    for _ in 0..20_000 {
        step(&mut system, &mut thermostat);
        let e = system.total_energy();
        worst = worst.max(((e + thermostat.energy(108) - h0) / h0).abs());
        wander = wander.max(((e - e0) / e0).abs());
    }
    assert!(worst <= 3e-3, "largest relative change {worst:e}");
    assert!(
        wander >= 0.02,
        "the energy of the atoms moved by only {wander:e}"
    );
}

#[test]
fn average_temperature_is_the_one_asked_for() {
    // The friction grows at (T_now / T - 1) / period^2, so over a time t
    // the average of T_now is T (1 + period^2 (xi_end - xi_start) / t):
    // the miss is set by how far the friction moved, and shrinks as 1 / t.
    // Over 20,000 steps and four seeds the average was 0.78620 to 0.78673
    // for T = 0.7867, a miss of 6.3e-4 of T or less, and it agreed with
    // the formula to 3e-6 of T. Tolerances: 2e-3 of T, and 2e-5 of T.
    let (mut system, mut thermostat) = liquid(1);
    let start = thermostat.friction;
    let mut sum = 0.0;
    for _ in 0..20_000 {
        step(&mut system, &mut thermostat);
        sum += thermo::temperature(&system.velocities);
    }
    let mean = sum / 20_000.0;
    assert!((mean / T - 1.0).abs() <= 2e-3, "mean temperature {mean}");
    let moved = thermostat.friction - start;
    // 20,000 steps of 0.005 tau are 100 tau.
    let want = T * (1.0 + thermostat.period * thermostat.period * moved / 100.0);
    assert!(
        (mean - want).abs() <= 2e-5 * T,
        "mean {mean}, formula {want}"
    );
}

#[test]
fn kinetic_energy_spreads_as_the_canonical_distribution_says() {
    // At a fixed temperature the kinetic energy K is not fixed: with g =
    // 3N - 3 free velocity components its variance is 2 / g of its mean
    // squared (a chi-squared distribution). Plain rescaling to T would
    // give zero. Over 200,000 steps and eight seeds the variance was 0.87
    // to 1.06 of that value. This seed: 1.06. Tolerance: 0.8 to 1.2.
    let (mut system, mut thermostat) = liquid(1);
    let (mut sum, mut squares) = (0.0, 0.0);
    for _ in 0..200_000 {
        step(&mut system, &mut thermostat);
        let k = system.kinetic_energy();
        sum += k;
        squares += k * k;
    }
    let mean = sum / 2e5;
    let ratio = (squares / 2e5 - mean * mean) / (mean * mean) / (2.0 / 321.0);
    assert!(
        (0.8..=1.2).contains(&ratio),
        "variance {ratio} of the canonical"
    );
}

#[test]
fn running_forward_then_backward_returns_to_the_start() {
    // With every velocity and the friction flipped, the same steps
    // retrace the path, thermostat included. Rounding grows as in
    // tests/conservation.rs. Measured over four seeds after 500 steps each
    // way: positions within 8.6e-12 sigma, velocities within 1.1e-10
    // sigma / tau, friction within 5.4e-13 and eta within 1.4e-13.
    // Tolerance: 1e-8 for each.
    let (mut system, mut thermostat) = liquid(1);
    let (r0, v0) = (system.positions.clone(), system.velocities.clone());
    let (friction, eta) = (thermostat.friction, thermostat.eta);
    for _ in 0..500 {
        step(&mut system, &mut thermostat);
    }
    assert!(
        (system.positions[0] - r0[0]).norm() > 0.1,
        "the atoms moved"
    );
    for v in &mut system.velocities {
        *v = -*v;
    }
    thermostat.friction = -thermostat.friction;
    for _ in 0..500 {
        step(&mut system, &mut thermostat);
    }
    for i in 0..108 {
        let miss = (system.positions[i] - r0[i]).norm();
        assert!(miss <= 1e-8, "atom {i} is {miss:e} sigma from its start");
        let miss = (system.velocities[i] + v0[i]).norm();
        assert!(miss <= 1e-8, "atom {i} velocity off by {miss:e}");
    }
    assert!((thermostat.friction + friction).abs() <= 1e-8);
    assert!((thermostat.eta - eta).abs() <= 1e-8);
}
