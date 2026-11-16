//! Forces from a neighbor list split across threads (forces.rs) against
//! forces from all pairs (pairs.rs), at one moment: the same to rounding,
//! the same to the last bit for any number of threads, and still right
//! after every atom has moved almost half the skin since the list was
//! built. Every tolerance is our own choice and says why.

use omega::forces::Method;
use omega::lattice;
use omega::pairs;
use omega::periodic::PeriodicBox;
use omega::rng::Rng;
use omega::vec3::Vec3;

const ZERO: Vec3 = Vec3::new(0.0, 0.0, 0.0);
const CUTOFF: f64 = 2.25;
const SKIN: f64 = 0.3;

/// An fcc lattice of `cells` cells a side at Rahman's density, each atom
/// moved at random by up to 0.1 sigma along each axis.
fn shaken(cells: usize) -> (Vec<Vec3>, PeriodicBox) {
    let space = PeriodicBox {
        side: 10.229 / 6.0 * cells as f64,
    };
    let mut rng = Rng::new(9);
    let mut positions = lattice::fcc(cells, space.side);
    for r in &mut positions {
        let (x, y, z) = (rng.uniform(), rng.uniform(), rng.uniform());
        *r = *r + (Vec3::new(x, y, z) - Vec3::new(0.5, 0.5, 0.5)) * 0.2;
    }
    (positions, space)
}

/// Forces from `method` against all pairs: the largest difference of a
/// force, in eps / sigma, and the relative differences of the potential
/// energy and of the virial.
fn compare(method: &mut Method, positions: &[Vec3], space: PeriodicBox) -> [f64; 3] {
    let mut want = vec![ZERO; positions.len()];
    let (u0, w0) = pairs::compute(positions, space, CUTOFF, &mut want);
    let mut got = vec![ZERO; positions.len()];
    let (u, w) = method.compute(positions, space, CUTOFF, &mut got);
    let pairs = got.iter().zip(&want);
    let force = pairs.map(|(&a, &b)| (a - b).norm()).fold(0.0, f64::max);
    [force, (u / u0 - 1.0).abs(), (w / w0 - 1.0).abs()]
}

#[test]
fn listed_forces_agree_with_all_pairs() {
    // 500 atoms (3 cells a side) and 108 (2 a side). The sums run in
    // another order, so they agree to rounding. Measured: forces within
    // 1.6e-14 eps / sigma, on forces up to 76, and energy and virial
    // within 2.3e-15 relative. Tolerances 1e-12, our own choice: well
    // above rounding, and far below the 0.081 eps / sigma that a single
    // pair missed just inside the cutoff would leave out.
    for cells in [5, 3] {
        let (positions, space) = shaken(cells);
        let misses = compare(&mut Method::listed(SKIN, 3), &positions, space);
        assert!(misses.iter().all(|&m| m <= 1e-12), "{misses:?}");
    }
}

#[test]
fn any_number_of_threads_gives_the_same_result_to_the_last_bit() {
    let (positions, space) = shaken(5);
    let run = |threads: usize| {
        let mut forces = vec![ZERO; positions.len()];
        let mut method = Method::listed(SKIN, threads);
        let (u, w) = method.compute(&positions, space, CUTOFF, &mut forces);
        (forces, u, w)
    };
    let one = run(1);
    for threads in [2, 3, 4, 7, 16] {
        assert!(run(threads) == one, "{threads} threads");
    }
}

#[test]
fn forces_stay_right_while_no_atom_has_moved_half_the_skin() {
    // Build the list, then move every atom 0.1499 sigma in a random
    // direction, just under half the skin. The list is not built again,
    // and every pair now within r_c must be on it, so the forces still
    // match all pairs. Measured: forces within 1.6e-13 eps / sigma, on
    // forces up to 1107, and energy and virial within 1.6e-15 relative.
    // Tolerances 1e-12.
    let (mut positions, space) = shaken(5);
    let mut method = Method::listed(SKIN, 2);
    compare(&mut method, &positions, space);
    let mut rng = Rng::new(7);
    for r in &mut positions {
        let d = Vec3::new(rng.normal(), rng.normal(), rng.normal());
        *r = *r + d * (0.1499 / d.norm());
    }
    let misses = compare(&mut method, &positions, space);
    assert_eq!(method.builds(), 1);
    assert!(misses.iter().all(|&m| m <= 1e-12), "{misses:?}");
}