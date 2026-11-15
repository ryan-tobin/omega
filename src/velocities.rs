//! Starting velocities for a chosen temperature.
//!
//! In thermal equilibrium each velocity component of an atom of mass m is
//! drawn from a normal distribution with mean 0 and variance kB * T / m,
//! the Maxwell-Boltzmann distribution. In reduced units m = kB = 1, so the
//! variance is T and the spread is sqrt(T).
//!
//! After drawing, the average velocity is subtracted from every atom, so
//! the total momentum is zero and the box as a whole does not drift. That
//! fixes 3 of the 3N velocity components, which is why the temperature of
//! the result is the sum of v^2 divided by 3N - 3, not 3N.
//!
//! Label: assumed. The starting velocities are put in by hand, drawn from
//! the distribution of classical statistical mechanics.

use crate::rng::Rng;
use crate::vec3::Vec3;

/// Velocities for `count` atoms at temperature `temperature` (in eps / kB),
/// with total momentum zero. Unit: sigma / tau.
pub fn maxwell_boltzmann(count: usize, temperature: f64, rng: &mut Rng) -> Vec<Vec3> {
    let spread = temperature.sqrt();
    let mut velocities = Vec::with_capacity(count);
    let mut sum = Vec3::new(0.0, 0.0, 0.0);
    for _ in 0..count {
        let x = rng.normal();
        let y = rng.normal();
        let z = rng.normal();
        let v = Vec3::new(x, y, z) * spread;
        sum = sum + v;
        velocities.push(v);
    }
    let mean = sum * (1.0 / count as f64);
    for v in &mut velocities {
        *v = *v - mean;
    }
    velocities
}

#[cfg(test)]
mod tests {
    use super::*;

    // Rahman's temperature, 94.4 K, in reduced units.
    const T: f64 = 0.7867;

    #[test]
    fn total_momentum_is_zero() {
        // Each atom has mass 1, so the total momentum is the sum of the
        // velocities. Subtracting the mean leaves only rounding: about
        // 1e-16 per atom. The tolerance, 1e-12, is our own choice.
        let mut rng = Rng::new(3);
        let velocities = maxwell_boltzmann(864, T, &mut rng);
        let mut p = Vec3::new(0.0, 0.0, 0.0);
        for v in velocities {
            p = p + v;
        }
        assert!(p.norm() <= 1e-12, "momentum {p:?}");
    }

    #[test]
    fn temperature_of_many_atoms_is_close_to_target() {
        // T = sum of v^2 / (3N - 3). Its relative standard error is
        // sqrt(2 / (3N)), 0.26% for N = 100,000. The tolerance, 1.5%, is
        // about 6 standard errors, our own choice.
        let n = 100_000;
        let mut rng = Rng::new(5);
        let velocities = maxwell_boltzmann(n, T, &mut rng);
        let sum_v2: f64 = velocities.iter().map(|v| v.norm_sq()).sum();
        let t = sum_v2 / (3 * n - 3) as f64;
        assert!((t / T - 1.0).abs() <= 0.015, "temperature {t}");
    }

    #[test]
    fn components_are_independent_with_variance_t() {
        // Each of x, y and z alone must have variance T, and no two of them
        // may be related: the averages of vx * vy, vy * vz and vz * vx are
        // 0. For N = 100,000 the standard errors are 0.0035 for a variance
        // and 0.0025 for a product. Tolerances: 2% of T and 0.012, about
        // 4.5 standard errors, our own choice.
        let n = 100_000;
        let mut rng = Rng::new(17);
        let velocities = maxwell_boltzmann(n, T, &mut rng);
        let mut sq = [0.0; 3];
        let mut cross = [0.0; 3];
        for v in velocities {
            sq = [sq[0] + v.x * v.x, sq[1] + v.y * v.y, sq[2] + v.z * v.z];
            cross = [
                cross[0] + v.x * v.y,
                cross[1] + v.y * v.z,
                cross[2] + v.z * v.x,
            ];
        }
        for s in sq {
            let variance = s / n as f64;
            assert!((variance / T - 1.0).abs() <= 0.02, "variance {variance}");
        }
        for c in cross {
            let mean = c / n as f64;
            assert!(mean.abs() <= 0.012, "average product {mean}");
        }
    }
}
