//! The Monte Carlo barostat (barostat.rs): the volume of atoms that never
//! touch against the ideal gas law, what a kept move changes, the smallest
//! box, carrying on from saved counts, and the pressure of liquid argon.

mod common;

use common::{T, liquid, step};
use omega::barostat::Barostat;
use omega::forces::Method;
use omega::pairs;
use omega::periodic::PeriodicBox;
use omega::rng::Rng;
use omega::system::System;
use omega::thermo;
use omega::vec3::Vec3;

const ZERO: Vec3 = Vec3::new(0.0, 0.0, 0.0);

/// `n` atoms at rest, spread at random in a box of side 2, with a
/// cutoff of 0.001 sigma, so that no two ever touch.
fn gas(n: usize, seed: u64) -> System {
    let mut rng = Rng::new(seed);
    let mut x = || 2.0 * rng.uniform();
    let positions = (0..n).map(|_| Vec3::new(x(), x(), x())).collect();
    System::new(positions, vec![ZERO; n], PeriodicBox { side: 2.0 }, 1e-3)
}

#[test]
fn atoms_that_never_touch_fill_the_volume_of_an_ideal_gas() {
    // With U = 0 the weights are L^(3N - 1) exp(-P L^3 / T), which
    // makes the volume a gamma variable with mean N T / P and
    // variance N (T / P)^2. Here T / P = 2, for 1 atom and for 5.
    // Over four seeds the mean of 100,000 moves was within 0.6% of
    // N T / P and the variance within 2%. Tolerances 3% and 10%, our
    // own choice. Leaving out the 2, or using N for N - 1, moves the
    // mean by 13% or more.
    for (n, max_step) in [(1, 1.5), (5, 0.6)] {
        let mut system = gas(n, 1);
        let mut barostat = Barostat::new(0.75, 1.5, max_step, 2, 0, 0);
        let (mut sum, mut squares) = (0.0, 0.0);
        for k in 0..110_000 {
            barostat.attempt(&mut system, Method::AllPairs);
            let volume = system.space.side.powi(3);
            if k >= 10_000 {
                sum += volume;
                squares += volume * volume;
            }
        }
        let mean = sum / 1e5;
        let variance = squares / 1e5 - mean * mean;
        let want = n as f64 * 2.0;
        assert!((mean / want - 1.0).abs() <= 0.03, "{n}: mean {mean}");
        let miss = variance / (2.0 * want) - 1.0;
        assert!(miss.abs() <= 0.1, "{n}: variance {variance}");
    }
}

#[test]
fn a_kept_move_scales_the_atoms_and_gives_them_their_new_forces() {
    // Two atoms 1.1 sigma apart. After the first kept move every
    // position is the old one times L' / L, and the forces, energy
    // and virial are those of the new positions. Exact.
    let positions = vec![Vec3::new(1.0, 1.0, 1.0), Vec3::new(2.1, 1.0, 1.0)];
    let space = PeriodicBox { side: 6.0 };
    let mut system = System::new(positions.clone(), vec![ZERO; 2], space, 2.5);
    let mut barostat = Barostat::new(1.0, 1.0, 0.1, 5, 0, 0);
    while !barostat.attempt(&mut system, Method::AllPairs) {
        assert_eq!(
            system.positions, positions,
            "a refused move changed the atoms"
        );
        assert_eq!(system.space, space);
    }
    let scale = system.space.side / 6.0;
    assert!(scale != 1.0 && (scale - 1.0).abs() <= 0.1 / 6.0);
    assert_eq!(
        system.positions,
        [positions[0] * scale, positions[1] * scale]
    );
    let mut forces = [ZERO; 2];
    let (u, w) = pairs::compute(&system.positions, system.space, 2.5, &mut forces);
    assert_eq!(system.forces, forces);
    assert_eq!((system.potential, system.virial), (u, w));
    assert_eq!(barostat.accepted, 1);
}

#[test]
fn box_never_shrinks_below_two_cutoffs() {
    // A pressure of 1000 squeezes two far-apart atoms as hard as it
    // can. The side must stop at 2 * 1.0, where pairs.rs would refuse.
    let positions = vec![Vec3::new(0.2, 0.2, 0.2), Vec3::new(1.5, 1.5, 1.5)];
    let space = PeriodicBox { side: 3.0 };
    let mut system = System::new(positions, vec![ZERO; 2], space, 1.0);
    let mut barostat = Barostat::new(1000.0, 1.0, 0.3, 9, 0, 0);
    for _ in 0..500 {
        barostat.attempt(&mut system, Method::AllPairs);
        assert!(system.space.side >= 2.0);
    }
    assert!(system.space.side < 2.1, "side {}", system.space.side);
}

#[test]
fn a_barostat_made_from_saved_counts_carries_on_the_same() {
    // 40 moves in one go, against 15 moves, a new barostat made from
    // the counts, and 25 more. Exact.
    let mut whole = gas(5, 3);
    let mut barostat = Barostat::new(0.75, 1.5, 0.6, 4, 0, 0);
    for _ in 0..40 {
        barostat.attempt(&mut whole, Method::AllPairs);
    }
    let mut parts = gas(5, 3);
    let mut first = Barostat::new(0.75, 1.5, 0.6, 4, 0, 0);
    for _ in 0..15 {
        first.attempt(&mut parts, Method::AllPairs);
    }
    let mut second = Barostat::new(0.75, 1.5, 0.6, 4, first.tries, first.accepted);
    for _ in 0..25 {
        second.attempt(&mut parts, Method::AllPairs);
    }
    assert_eq!(parts.space, whole.space);
    assert_eq!(parts.positions, whole.positions);
    assert_eq!((second.tries, second.accepted), (40, barostat.accepted));
}

#[test]
fn pressure_from_the_forces_averages_the_barostat_s_pressure() {
    // The barostat never looks at the pressure. It only compares
    // energies. The pressure measured from the forces (thermo.rs) must
    // still average the pressure it was given. A try every 10 steps for
    // 60,000 steps, the first 10,000 left out. Over four seeds the
    // average was 0.987 to 1.011 for P = 1, each with a standard error
    // of 0.010 to 0.019 from ten blocks, at a density of 0.81, and 58%
    // to 59% of the moves were kept. This seed: 1.002. Tolerance: 0.05.
    let (mut system, mut thermostat) = liquid(1);
    let mut barostat = Barostat::new(1.0, T, 0.05, 101, 0, 0);
    let mut sum = 0.0;
    for k in 1..=60_000 {
        step(&mut system, &mut thermostat);
        if k % 10 == 0 {
            barostat.attempt(&mut system, Method::AllPairs);
        }
        if k > 10_000 {
            sum += thermo::pressure(&system.velocities, system.virial, system.space);
        }
    }
    let mean = sum / 50_000.0;
    assert!((mean - 1.0).abs() <= 0.05, "mean pressure {mean}");
    let kept = barostat.accepted as f64 / barostat.tries as f64;
    assert!((0.4..=0.7).contains(&kept), "{kept} of the moves kept");
}
