//! The Lennard-Jones energy and force between two atoms.
//!
//! Everything here is in reduced units (see units.rs): distances in sigma,
//! energies in eps, forces in eps / sigma. In those units sigma = eps = 1
//! and the potential for two atoms a distance r apart is
//! U(r) = 4 * (1/r^12 - 1/r^6).
//!
//! The 1/r^12 term pushes atoms apart when they overlap. The 1/r^6 term
//! pulls them together at longer range.
//!
//! Label: approximated. A pair potential fitted to argon. It ignores
//! three-body effects and quantum behavior.

use crate::vec3::Vec3;

/// Potential energy of the pair, U = 4 * (1/r^12 - 1/r^6), in eps.
/// `r_ij` is the separation r_i - r_j in sigma. It must not be zero.
pub fn energy(r_ij: Vec3) -> f64 {
    let inv_r2 = 1.0 / r_ij.norm_sq();
    let inv_r6 = inv_r2 * inv_r2 * inv_r2;
    4.0 * (inv_r6 * inv_r6 - inv_r6)
}

/// Force on atom i from atom j, in eps / sigma:
/// F_ij = 24 * (2/r^12 - 1/r^6) * r_ij / r^2.
///
/// `r_ij` is the separation r_i - r_j in sigma. It must not be zero.
/// The force on atom j is the negative of this.
///
/// Label: derived. It is minus the slope of `energy`: F = -dU/dr, along
/// the line that joins the atoms.
pub fn force(r_ij: Vec3) -> Vec3 {
    let inv_r2 = 1.0 / r_ij.norm_sq();
    let inv_r6 = inv_r2 * inv_r2 * inv_r2;
    r_ij * (24.0 * (2.0 * inv_r6 * inv_r6 - inv_r6) * inv_r2)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Tolerance at the minimum. Our own choice. 2^(1/6) cannot be held
    // exactly in an f64, and near the minimum the force changes by about
    // 57 per unit of distance, so rounding can leave a force near 1e-14.
    // 1e-12 gives that a margin of 100.
    const MIN_TOL: f64 = 1e-12;

    // Step and tolerance for the slope test. Our own choices. The miss of
    // a slope measured over a step H each way shrinks as H^2. The worst
    // miss measured with H = 1e-5 was 4.3e-8, on a force of 54. 1e-6
    // gives that a margin of 20.
    const H: f64 = 1e-5;
    const SLOPE_TOL: f64 = 1e-6;

    fn assert_close(got: f64, want: f64, tol: f64) {
        let miss = (got - want).abs();
        assert!(miss <= tol, "got {got:e}, want {want:e}, off by {miss:e}");
    }

    // At separations of 1 and 2 every value below is a whole number times
    // a power of two, which an f64 holds exactly. Those four tests compare
    // with == and the tolerance is zero.

    #[test]
    fn energy_is_zero_at_one_sigma() {
        // 4 * (1 - 1) = 0.
        assert_eq!(energy(Vec3::new(1.0, 0.0, 0.0)), 0.0);
    }

    #[test]
    fn energy_at_two_sigma_is_minus_63_over_1024() {
        // 4 * (1/4096 - 1/64) = -63/1024.
        assert_eq!(energy(Vec3::new(0.0, 2.0, 0.0)), -63.0 / 1024.0);
    }

    #[test]
    fn force_at_one_sigma_is_24_pushing_i_away_from_j() {
        // 24 * (2 - 1) * 1 / 1 = 24, along +x, which points from j to i.
        let f = force(Vec3::new(1.0, 0.0, 0.0));
        assert_eq!(f, Vec3::new(24.0, 0.0, 0.0));
    }

    #[test]
    fn force_at_two_sigma_pulls_i_toward_j() {
        // 24 * (2/4096 - 1/64) * 2 / 4 = -93/512, along -z, toward j.
        let f = force(Vec3::new(0.0, 0.0, 2.0));
        assert_eq!(f, Vec3::new(0.0, 0.0, -93.0 / 512.0));
    }

    #[test]
    fn force_on_j_is_minus_force_on_i() {
        // Newton's third law. Swapping the atoms flips r_ij and leaves r^2
        // alone, so the two forces are exactly opposite.
        let r = Vec3::new(0.9, -0.5, 0.7);
        assert_eq!(force(-r), -force(r));
    }

    #[test]
    fn minimum_is_at_2_to_the_one_sixth_with_depth_minus_one() {
        let r_min = 2.0_f64.powf(1.0 / 6.0);
        let at_min = Vec3::new(r_min, 0.0, 0.0);
        assert_close(energy(at_min), -1.0, MIN_TOL);
        assert_close(force(at_min).norm(), 0.0, MIN_TOL);
        // A minimum, not a maximum: the energy is higher on both sides.
        assert!(energy(Vec3::new(r_min - 0.01, 0.0, 0.0)) > energy(at_min));
        assert!(energy(Vec3::new(r_min + 0.01, 0.0, 0.0)) > energy(at_min));
    }

    #[test]
    fn force_matches_slope_of_energy() {
        // Move atom i a distance H each way along one axis and measure
        //     slope = (U(r + H) - U(r - H)) / (2 * H)
        // The force along that axis must be minus the slope.
        let separations = [
            Vec3::new(0.5, 0.6, -0.55), // r = 0.955, strongly repelling
            Vec3::new(0.9, -0.5, 0.7),  // r = 1.245, attracting
            Vec3::new(-1.5, 1.0, 1.2),  // r = 2.166, weakly attracting
        ];
        let axes = [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];
        for r in separations {
            let f = force(r);
            for (axis, f_part) in axes.into_iter().zip([f.x, f.y, f.z]) {
                let step = axis * H;
                let slope = (energy(r + step) - energy(r - step)) / (2.0 * H);
                assert_close(f_part, -slope, SLOPE_TOL);
            }
        }
    }
}
