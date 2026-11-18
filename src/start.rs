//! The atoms a scenario starts from (scenario.rs).
//!
//! They sit on an fcc lattice that fills the box, or with `fill` only the
//! part of it below a fraction of the side along x: a slab, with empty
//! space for its vapor. Their velocities are drawn for `temperature`
//! (velocities.rs). With `hot`, the atoms at or beyond a fraction of the
//! side become a liquid next to the crystal: only a share of them is
//! kept, picked at random, because a liquid is less dense than its
//! crystal, and they are drawn for a second, higher temperature. Left to
//! move freely for a short while, that part melts and the rest does not.

use crate::lattice;
use crate::rng::Rng;
use crate::scenario::Scenario;
use crate::vec3::Vec3;
use crate::velocities;

/// The starting positions, in sigma, and velocities, in sigma / tau.
pub fn atoms(s: &Scenario) -> (Vec<Vec3>, Vec<Vec3>) {
    // A lattice plane counts as inside a part of the box if it lies a
    // quarter of a cell inside it, so a plane on the border is never in doubt.
    let quarter = 0.25 * s.side / s.cells as f64;
    let inside = |r: &Vec3, part: f64| r.x + quarter < part * s.side;
    let mut positions = lattice::fcc(s.cells, s.side);
    positions.retain(|r| inside(r, s.fill));
    let mut rng = Rng::new(s.seed);
    if let Some((from, _, keep)) = s.hot {
        positions.retain(|r| inside(r, from) || rng.uniform() < keep);
    }
    let mut v = velocities::maxwell_boltzmann(positions.len(), s.temperature, &mut rng);
    if let Some((from, hot, _)) = s.hot {
        // Scaling a draw for T by sqrt(hot / T) makes it a draw for `hot`.
        let scale = (hot / s.temperature).sqrt();
        for (v, r) in v.iter_mut().zip(&positions) {
            if !inside(r, from) {
                *v = *v * scale;
            }
        }
        // The atoms as a whole must stay at rest.
        let sum = v.iter().fold(Vec3::new(0.0, 0.0, 0.0), |sum, &v| sum + v);
        let mean = sum * (1.0 / v.len() as f64);
        for v in &mut v {
            *v = *v - mean;
        }
    }
    (positions, v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scenario;

    /// A scenario of 4 cells in a box of side 8, so the lattice planes
    /// across x sit at x = 0, 1, ... 7 with 32 atoms each, plus `extra`.
    fn planes(extra: &str) -> Scenario {
        let text = "cells 4\nside 8\ntemperature 0.5\nseed 3\ncutoff 2.5\ndt 0.005\n";
        scenario::parse(&format!("{text}{extra}\nstage only\nsteps 1")).unwrap()
    }

    #[test]
    fn fill_keeps_the_planes_below_its_fraction_of_the_side() {
        // No fill: all 256 atoms. fill 0.5: the planes at x = 0 to 3, 128
        // atoms, the one at x = 4 on the border left out. fill 0.3: below
        // 2.4, the planes at 0 and 1, 64 atoms. Exact.
        assert_eq!(atoms(&planes("")).0.len(), 256);
        let (half, velocities) = atoms(&planes("fill 0.5"));
        assert_eq!((half.len(), velocities.len()), (128, 128));
        assert!(half.iter().all(|r| r.x <= 3.0));
        assert_eq!(atoms(&planes("fill 0.3")).0.len(), 64);
    }

    #[test]
    fn hot_atoms_start_hot_and_the_whole_stays_at_rest() {
        // hot 0.5 8 1: the 128 atoms at x >= 4 are all kept and drawn for
        // T = 8, the other 128 for 0.5. The mean of v^2 / 3 over 128 atoms
        // scatters by sqrt(2 / 384) = 7% around T. Measured for this seed:
        // 0.510 and 7.88. Tolerance: 20%. The total momentum is zero to
        // rounding: measured 3.4e-14, tolerance 1e-12.
        let (positions, velocities) = atoms(&planes("hot 0.5 8 1"));
        let mut sums = [0.0, 0.0];
        let mut total = Vec3::new(0.0, 0.0, 0.0);
        for (r, v) in positions.iter().zip(&velocities) {
            sums[usize::from(r.x >= 4.0)] += v.norm_sq() / (3.0 * 128.0);
            total = total + *v;
        }
        assert!(
            (sums[0] / 0.5 - 1.0).abs() <= 0.2,
            "cold half at {}",
            sums[0]
        );
        assert!(
            (sums[1] / 8.0 - 1.0).abs() <= 0.2,
            "hot half at {}",
            sums[1]
        );
        assert!(total.norm() <= 1e-12, "total momentum {:e}", total.norm());
    }

    #[test]
    fn hot_part_keeps_its_share_of_the_atoms() {
        // hot 0.5 8 0.75: each of the 128 atoms at x >= 4 stays with a
        // chance of 3 in 4, so about 96 stay, give or take 5 (the spread
        // of 128 coin flips). This seed: 102. All 128 below x = 4 stay.
        let (positions, _) = atoms(&planes("hot 0.5 8 0.75"));
        let hot = positions.iter().filter(|r| r.x >= 4.0).count();
        assert_eq!(positions.len() - hot, 128);
        assert_eq!(hot, 102);
    }
}
