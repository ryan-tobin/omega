//! Checks of the measurements through the library's public functions, the
//! way the measure program uses them: the mean squared displacement of
//! atoms that cross the faces of the box, the virial the system keeps,
//! g(r) of atoms placed at random, the slope behind the diffusion
//! coefficient and the wrapping of XYZ positions. Every tolerance is our
//! own choice and says why.

use omega::msd::{self, Msd};
use omega::pairs;
use omega::periodic::PeriodicBox;
use omega::rdf::Rdf;
use omega::rng::Rng;
use omega::system::System;
use omega::units::ARGON;
use omega::vec3::Vec3;
use omega::xyz;
use std::f64::consts::PI;

const ZERO: Vec3 = Vec3::new(0.0, 0.0, 0.0);

#[test]
fn msd_counts_the_whole_path_across_the_box() {
    // Two atoms far apart in a box of side 8, both moving at the same
    // velocity, so they never come within the cutoff and nothing pushes
    // them. After 1000 steps of 0.01 tau, t = 10 and each has gone
    // (10, 5, 2.5) sigma, out through a face. Positions are never wrapped,
    // so MSD = |v|^2 t^2 = 1.3125 * 100 = 131.25. Rounding in the 1000
    // position updates is measured at 2.5e-12. Tolerance: 1e-9.
    let v = Vec3::new(1.0, 0.5, 0.25);
    let positions = vec![Vec3::new(1.0, 1.0, 1.0), Vec3::new(5.0, 5.0, 5.0)];
    let mut system = System::new(positions, vec![v, v], PeriodicBox { side: 8.0 }, 2.5);
    let msd = Msd::new(&system.positions);
    for _ in 0..1000 {
        system.step(0.01);
    }
    let m = msd.measure(&system.positions);
    assert!((m - 131.25).abs() <= 1e-9, "MSD {m}");
}

#[test]
fn system_keeps_the_virial_of_its_current_positions() {
    // Two atoms at rest 1 sigma apart. There the force is 24 eps / sigma,
    // pushing them apart, so r . F = 24 exactly (pairs.rs). After one
    // step they are further apart, the push is weaker and the virial is
    // below 24. It must be the virial of the new positions, the same
    // number pairs::compute gives for them.
    let positions = vec![Vec3::new(3.5, 4.0, 4.0), Vec3::new(4.5, 4.0, 4.0)];
    let mut system = System::new(positions, vec![ZERO; 2], PeriodicBox { side: 8.0 }, 2.5);
    assert_eq!(system.virial, 24.0);
    system.step(0.005);
    let mut forces = vec![ZERO; 2];
    let (_, w) = pairs::compute(&system.positions, system.space, system.cutoff, &mut forces);
    assert_eq!(system.virial, w);
    assert!(w < 24.0, "virial {w}");
}

#[test]
#[should_panic(expected = "r_max is over half the box")]
fn rdf_refuses_r_max_over_half_the_box() {
    // Past L/2 a pair would also be counted at a farther copy, which the
    // nearest-image rule never gives, so g(r) would be wrong there.
    Rdf::new(PeriodicBox { side: 3.0 }, 1.6, 10);
}

#[test]
fn ten_atoms_placed_at_random_give_g_of_0_9() {
    // 10 atoms placed evenly at random in a box of side 5, 20,000 times.
    // The nearest image of a pair is then evenly spread over the box, so a
    // shell of volume V_k inside r <= L/2 holds each pair with chance
    // V_k / V. Bin k expects 20,000 * 45 V_k / 125 pairs and g = 1 - 1/10.
    // The bins from 0.5 to 2.5, the last one too, expect 2,700 to 54,000
    // pairs, a relative standard error of 1 / sqrt(expected), 1.9% at
    // most. The worst bin is measured 1.8 standard errors off.
    // Tolerance: 5 standard errors in each bin.
    let space = PeriodicBox { side: 5.0 };
    let mut rng = Rng::new(21);
    let mut rdf = Rdf::new(space, 2.5, 25);
    for _ in 0..20_000 {
        let mut positions = Vec::with_capacity(10);
        for _ in 0..10 {
            let x = rng.uniform();
            let y = rng.uniform();
            let z = rng.uniform();
            positions.push(Vec3::new(x, y, z) * space.side);
        }
        rdf.add(&positions);
    }
    for (r, g) in rdf.g().into_iter().skip(5) {
        let (lower, upper) = (r - 0.05, r + 0.05);
        let shell = 4.0 / 3.0 * PI * (upper * upper * upper - lower * lower * lower);
        let expected = 20_000.0 * 45.0 * shell / 125.0;
        let error = 1.0 / expected.sqrt();
        assert!((g / 0.9 - 1.0).abs() <= 5.0 * error, "g({r}) = {g}");
    }
}

#[test]
fn scattered_points_give_the_least_squares_slope() {
    // Through (0, 0), (1, 6), (2, 6), (3, 18) the least-squares line has
    // slope sum (t - 1.5)(m - 7.5) / sum (t - 1.5)^2 = 27 / 5 = 5.4, so
    // D = 0.9. A line through the first and last points alone would give
    // 1. Tolerance 1e-12, for rounding.
    let d = msd::diffusion(&[0.0, 1.0, 2.0, 3.0], &[0.0, 6.0, 6.0, 18.0]);
    assert!((d - 0.9).abs() <= 1e-12, "D {d}");
}

#[test]
fn xyz_wraps_every_axis_into_the_box() {
    // An atom at (-1, 9, 17) sigma in a box of side 8 is drawn at
    // (7, 1, 1) sigma, which is (23.8, 3.4, 3.4) angstrom.
    let mut out = Vec::new();
    let atom = [Vec3::new(-1.0, 9.0, 17.0)];
    xyz::write_frame(&mut out, "Ar", &atom, PeriodicBox { side: 8.0 }, 0.0, ARGON).unwrap();
    let text = String::from_utf8(out).unwrap();
    let line = text.lines().nth(2).unwrap();
    assert_eq!(line, "Ar 23.800000 3.400000 3.400000");
}
