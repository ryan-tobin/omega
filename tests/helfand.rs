//! The stress and heat currents of flux.rs as the rates of change of the
//! Helfand moments along a real motion: the sum of y v_x (and its turns,
//! z v_y and x v_z) for the shear, and the sum of r_i e_i for the heat.

use omega::flux;
use omega::forces::Method;
use omega::lattice;
use omega::lj;
use omega::periodic::PeriodicBox;
use omega::rng::Rng;
use omega::system::System;
use omega::vec3::Vec3;
use omega::velocities;

const ZERO: Vec3 = Vec3::new(0.0, 0.0, 0.0);
const CUTOFF: f64 = 2.5;

/// The Helfand moments of a cluster: sum of (y v_x, z v_y, x v_z), sum
/// of (x v_x, y v_y, z v_z), and the energy moment, sum of r_i e_i.
fn moments(system: &System) -> (Vec3, Vec3, Vec3) {
    let (r, v) = (&system.positions, &system.velocities);
    let shift = lj::energy(Vec3::new(CUTOFF, 0.0, 0.0));
    let (mut g, mut d, mut h) = (ZERO, ZERO, ZERO);
    for i in 0..r.len() {
        g = g + Vec3::new(r[i].y * v[i].x, r[i].z * v[i].y, r[i].x * v[i].z);
        d = d + Vec3::new(r[i].x * v[i].x, r[i].y * v[i].y, r[i].z * v[i].z);
        let mut e = 0.5 * v[i].norm_sq();
        for j in (0..r.len()).filter(|&j| j != i) {
            if (r[i] - r[j]).norm() < CUTOFF {
                e += 0.5 * (lj::energy(r[i] - r[j]) - shift);
            }
        }
        h = h + r[i] * e;
    }
    (g, d, h)
}

/// The pairs closer than the cutoff, in a box no pair reaches across.
fn inside(system: &System) -> Vec<(usize, usize)> {
    let r = &system.positions;
    let n = r.len();
    let close = |&(i, j): &(usize, usize)| (r[i] - r[j]).norm() < CUTOFF;
    (0..n)
        .flat_map(|i| (i + 1..n).map(move |j| (i, j)))
        .filter(close)
        .collect()
}

#[test]
fn currents_are_the_rates_of_the_helfand_moments() {
    // 32 atoms of a 2-cell lattice, spaced for density 0.8442, shaken
    // and heated to T = 1, in the middle of a box of side 40 so that no
    // pair reaches across its edge. Steps of 1e-5 tau: the change of a
    // moment over two steps, against the current at the step between.
    // Velocity Verlet and the difference both miss by about dt^2 (at
    // 2e-5 tau the shear miss grew 4 times). Where a pair crosses the
    // cutoff the force jumps, and so does the rate, so those steps are
    // skipped. Measured: at most 7.4e-8 of the size of the shear current,
    // 3.1e-8 of the diagonal stress and 2.3e-6 of the heat current, 2 of
    // 500 steps skipped. Tolerance 1e-5, our own choice; a missing half on
    // any term misses by over 0.1.
    let side = 2.0 * (4.0_f64 / 0.8442).cbrt();
    let mut rng = Rng::new(5);
    let mut positions = lattice::fcc(2, side);
    for r in &mut positions {
        let d = Vec3::new(rng.uniform(), rng.uniform(), rng.uniform());
        *r = *r + (d - Vec3::new(0.5, 0.5, 0.5)) * 0.2 + Vec3::new(18.0, 18.0, 18.0);
    }
    let v = velocities::maxwell_boltzmann(positions.len(), 1.0, &mut rng);
    let space = PeriodicBox { side: 40.0 };
    let mut system = System::with_forces(positions, v, space, CUTOFF, Method::listed(0.3, 1));
    let dt = 1e-5;
    let mut checked = 0;
    for k in 0..500 {
        let (before, pairs) = (moments(&system), inside(&system));
        system.step(dt);
        let Method::Listed { list, .. } = &system.method else {
            panic!("the system was built with a list");
        };
        let p = &system.positions;
        let now = flux::measure(p, &system.velocities, space, CUTOFF, list, 1);
        system.step(dt);
        if inside(&system) != pairs {
            continue;
        }
        let after = moments(&system);
        let rate_g = (after.0 - before.0) * (0.5 / dt);
        let rate_d = (after.1 - before.1) * (0.5 / dt);
        let rate_h = (after.2 - before.2) * (0.5 / dt);
        let miss_g = (rate_g - now.shear).norm() / now.shear.norm();
        let miss_d = (rate_d - now.normal).norm() / now.normal.norm();
        let miss_h = (rate_h - now.heat).norm() / now.heat.norm();
        assert!(miss_g <= 1e-5, "step {k}: shear {rate_g:?} {:?}", now.shear);
        assert!(miss_d <= 1e-5, "step {k}: diagonal {:?}", now.normal);
        assert!(miss_h <= 1e-5, "step {k}: heat {rate_h:?} {:?}", now.heat);
        checked += 1;
    }
    assert!(checked >= 450, "only {checked} steps checked");
}
