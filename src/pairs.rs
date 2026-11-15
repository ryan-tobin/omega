//! Forces and potential energy from all pairs of atoms in the repeating box.
//!
//! Every pair i < j is visited once, so the work grows as N^2: 372,816
//! pairs for 864 atoms. The separation of each pair is moved to its nearest
//! copy (periodic.rs). Pairs farther apart than the cutoff r_c are skipped.
//! Each pair inside it adds U(r) - U(r_c) to the energy. That is the cut and
//! shifted potential: it reaches zero at r_c, so the energy does not jump
//! when a pair crosses the cutoff. The force is the plain Lennard-Jones
//! force, which the shift does not change.
//!
//! Label: approximated. Pairs beyond r_c are ignored. At 2.25 sigma the
//! pair energy is -0.031 eps, about 3% of the well depth.

use crate::lj;
use crate::periodic::PeriodicBox;
use crate::vec3::Vec3;

/// Fills `forces` with the total force on each atom, in eps / sigma, and
/// returns the total potential energy, in eps. `positions` are in sigma,
/// and `forces` must have one entry per position. `cutoff` is r_c in
/// sigma and must be at most half the box side.
pub fn compute(positions: &[Vec3], space: PeriodicBox, cutoff: f64, forces: &mut [Vec3]) -> f64 {
    assert_eq!(forces.len(), positions.len());
    assert!(cutoff <= 0.5 * space.side, "cutoff is over half the box");
    let cutoff_sq = cutoff * cutoff;
    let shift = lj::energy(Vec3::new(cutoff, 0.0, 0.0));
    forces.fill(Vec3::new(0.0, 0.0, 0.0));
    let mut potential = 0.0;
    for i in 0..positions.len() {
        for j in (i + 1)..positions.len() {
            let r_ij = space.nearest_image(positions[i] - positions[j]);
            if r_ij.norm_sq() < cutoff_sq {
                // Newton's third law: j feels the opposite of what i feels.
                let f = lj::force(r_ij);
                forces[i] = forces[i] + f;
                forces[j] = forces[j] - f;
                potential += lj::energy(r_ij) - shift;
            }
        }
    }
    potential
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice;
    use crate::rng::Rng;

    const BOX: PeriodicBox = PeriodicBox { side: 8.0 };
    const ZERO: Vec3 = Vec3::new(0.0, 0.0, 0.0);

    #[test]
    fn atoms_interact_across_the_boundary() {
        // x = 0.25 and x = 7.25: the copy of the second atom at -0.75 is
        // exactly 1 sigma away, where U = 0 and the force is 24 (lj.rs).
        let positions = [Vec3::new(0.25, 4.0, 4.0), Vec3::new(7.25, 4.0, 4.0)];
        let mut forces = [ZERO; 2];
        let u = compute(&positions, BOX, 2.5, &mut forces);
        assert_eq!(
            forces,
            [Vec3::new(24.0, 0.0, 0.0), Vec3::new(-24.0, 0.0, 0.0)]
        );
        assert_eq!(u, -lj::energy(Vec3::new(2.5, 0.0, 0.0)));
    }

    #[test]
    fn pairs_beyond_the_cutoff_add_nothing() {
        let positions = [Vec3::new(1.0, 1.0, 1.0), Vec3::new(1.0, 4.0, 1.0)];
        let mut forces = [Vec3::new(9.0, 9.0, 9.0); 2];
        let u = compute(&positions, BOX, 2.5, &mut forces);
        assert_eq!(forces, [ZERO; 2]);
        assert_eq!(u, 0.0);
    }

    #[test]
    #[should_panic(expected = "cutoff is over half the box")]
    fn cutoff_over_half_the_box_is_refused() {
        let positions = [Vec3::new(1.0, 1.0, 1.0), Vec3::new(1.0, 4.0, 1.0)];
        let mut forces = [ZERO; 2];
        compute(&positions, BOX, 4.1, &mut forces);
    }

    #[test]
    fn forces_add_up_to_zero() {
        // 108 lattice atoms, each moved at random by up to 0.1 sigma along
        // each axis. Every pair adds f to one atom and -f to the other, so
        // the sum is zero apart from rounding. Measured: 7.7e-14 eps / sigma,
        // on forces of about 10. Tolerance: 1e-10, our own choice.
        let space = PeriodicBox { side: 10.229 / 2.0 };
        let mut rng = Rng::new(9);
        let mut positions = lattice::fcc(3, space.side);
        for r in &mut positions {
            let x = rng.uniform();
            let y = rng.uniform();
            let z = rng.uniform();
            *r = *r + (Vec3::new(x, y, z) - Vec3::new(0.5, 0.5, 0.5)) * 0.2;
        }
        let mut forces = vec![ZERO; positions.len()];
        compute(&positions, space, 2.25, &mut forces);
        let total = forces.iter().fold(ZERO, |sum, &f| sum + f);
        assert!(total.norm() <= 1e-10, "total force {total:?}");
    }

    #[test]
    fn perfect_lattice_matches_the_shell_sum() {
        // Rahman's density with 108 atoms: L = 10.229 / 2, a = L / 3. With
        // r_c = 2.25 the first three shells are inside the cutoff, at
        // a / sqrt(2), a and a * sqrt(3/2), holding 12, 6 and 24 atoms
        // (lattice.rs). The fourth, at a * sqrt(2) = 2.41, is outside.
        // Each pair is shared by two atoms, so the energy per atom is
        // (12 U(r1) + 6 U(r2) + 24 U(r3)) / 2, with U shifted. By symmetry
        // the force on every atom is zero. Measured: 1.9e-14 relative for
        // the energy and 1.5e-15 for the largest force. Tolerances: 1e-12,
        // our own choice.
        let space = PeriodicBox { side: 10.229 / 2.0 };
        let cutoff = 2.25;
        let positions = lattice::fcc(3, space.side);
        let mut forces = vec![ZERO; positions.len()];
        let total = compute(&positions, space, cutoff, &mut forces);
        let a = space.side / 3.0;
        let shift = lj::energy(Vec3::new(cutoff, 0.0, 0.0));
        let u = |r: f64| lj::energy(Vec3::new(r, 0.0, 0.0)) - shift;
        let r1 = a / 2.0_f64.sqrt();
        let r3 = a * 1.5_f64.sqrt();
        let per_atom = (12.0 * u(r1) + 6.0 * u(a) + 24.0 * u(r3)) / 2.0;
        let want = per_atom * positions.len() as f64;
        let miss = (total - want).abs() / want.abs();
        assert!(miss <= 1e-12, "got {total}, want {want}");
        for f in forces {
            assert!(f.norm() <= 1e-12, "force {f:?}");
        }
    }
}
