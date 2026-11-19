//! The stress and heat currents of flux.rs against a plain sum over all
//! pairs in the repeating box. tests/helfand.rs checks them as rates.

use omega::flux::{self, Flux};
use omega::lattice;
use omega::lj;
use omega::neighbors::NeighborList;
use omega::periodic::PeriodicBox;
use omega::rng::Rng;
use omega::vec3::Vec3;
use omega::velocities;

const ZERO: Vec3 = Vec3::new(0.0, 0.0, 0.0);
const CUTOFF: f64 = 2.5;

/// 108 lattice atoms at density 0.8442, each moved at random by up to
/// 0.1 sigma along each axis, with velocities for T = 0.722.
fn shaken(seed: u64) -> (Vec<Vec3>, Vec<Vec3>, PeriodicBox) {
    let space = PeriodicBox {
        side: (108.0_f64 / 0.8442).cbrt(),
    };
    let mut rng = Rng::new(seed);
    let mut positions = lattice::fcc(3, space.side);
    for r in &mut positions {
        let d = Vec3::new(rng.uniform(), rng.uniform(), rng.uniform());
        *r = *r + (d - Vec3::new(0.5, 0.5, 0.5)) * 0.2;
    }
    let v = velocities::maxwell_boltzmann(positions.len(), 0.722, &mut rng);
    (positions, v, space)
}

/// The currents straight from their definitions, pair by pair, with each
/// pair i < j visited once.
fn by_hand(positions: &[Vec3], v: &[Vec3], space: PeriodicBox) -> Flux {
    let shift = lj::energy(Vec3::new(CUTOFF, 0.0, 0.0));
    let mut e: Vec<f64> = v.iter().map(|v| 0.5 * v.norm_sq()).collect();
    let (mut shear, mut normal, mut passed, mut pairs) = (ZERO, ZERO, ZERO, 0);
    for a in v {
        shear = shear + Vec3::new(a.x * a.y, a.y * a.z, a.z * a.x);
        normal = normal + Vec3::new(a.x * a.x, a.y * a.y, a.z * a.z);
    }
    for i in 0..positions.len() {
        for j in (i + 1)..positions.len() {
            let r = space.nearest_image(positions[i] - positions[j]);
            if r.norm() < CUTOFF {
                let f = lj::force(r);
                shear = shear + Vec3::new(r.x * f.y, r.y * f.z, r.z * f.x);
                normal = normal + Vec3::new(r.x * f.x, r.y * f.y, r.z * f.z);
                let u = lj::energy(r) - shift;
                e[i] += 0.5 * u;
                e[j] += 0.5 * u;
                // Pair i, j passes r (F . v_i) / 2 and, seen from j,
                // (-r) (-F . v_j) / 2.
                passed = passed + r * (0.5 * f.dot(v[i] + v[j]));
                pairs += 1;
            }
        }
    }
    let carried = (0..v.len()).fold(ZERO, |sum, i| sum + v[i] * e[i]);
    Flux {
        shear,
        normal,
        heat: carried + passed,
        pairs,
    }
}

#[test]
fn currents_match_a_sum_over_all_pairs() {
    // Only the order of the sums differs. Measured: at most 5.4e-15 of
    // the size of each current. Tolerance 1e-12 relative, our own choice.
    for seed in [1, 2, 3] {
        let (positions, v, space) = shaken(seed);
        let mut list = NeighborList::new(0.3);
        list.build(&positions, space, CUTOFF, 2);
        let got = flux::measure(&positions, &v, space, CUTOFF, &list, 2);
        let want = by_hand(&positions, &v, space);
        assert_eq!(got.pairs, want.pairs);
        let miss = |a: Vec3, b: Vec3| (a - b).norm() / b.norm();
        assert!(miss(got.shear, want.shear) <= 1e-12, "{got:?} {want:?}");
        assert!(miss(got.normal, want.normal) <= 1e-12, "{got:?} {want:?}");
        assert!(miss(got.heat, want.heat) <= 1e-12, "{got:?} {want:?}");
    }
}

#[test]
fn currents_are_the_same_on_any_number_of_threads() {
    let (positions, v, space) = shaken(4);
    let mut list = NeighborList::new(0.3);
    list.build(&positions, space, CUTOFF, 1);
    let one = flux::measure(&positions, &v, space, CUTOFF, &list, 1);
    for threads in [2, 3, 7] {
        assert_eq!(
            flux::measure(&positions, &v, space, CUTOFF, &list, threads),
            one
        );
    }
}
