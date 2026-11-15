//! Temperature and pressure, and the velocity rescaling that brings the
//! atoms to a chosen temperature before anything is measured.
//!
//! Units are reduced (units.rs). Each atom's mass is 1 and kB is 1, so a
//! temperature is in eps / kB (120 K for argon) and a pressure is in
//! eps / sigma^3 (42.15 MPa for argon).

use crate::periodic::PeriodicBox;
use crate::vec3::Vec3;

/// Sum of v^2 over all atoms, twice the kinetic energy, in eps.
fn sum_v2(velocities: &[Vec3]) -> f64 {
    velocities.iter().map(|v| v.norm_sq()).sum()
}

/// Temperature from the kinetic energy, T = sum of v^2 / (3N - 3), in
/// eps / kB.
///
/// On average each velocity component carries kB T / 2 of kinetic energy
/// (equipartition). Of the 3N components, 3 are fixed because the total
/// momentum is zero, which leaves 3N - 3.
///
/// Label: derived, from classical statistical mechanics.
pub fn temperature(velocities: &[Vec3]) -> f64 {
    assert!(velocities.len() >= 2, "temperature needs two or more atoms");
    sum_v2(velocities) / (3 * velocities.len() - 3) as f64
}

/// Pressure from the virial theorem, P = (sum of v^2 + W) / (3V), in
/// eps / sigma^3. V = L^3 is the volume of the box. W is the virial, the
/// sum over pairs of r_ij . F_ij, in eps (pairs.rs).
///
/// The first part is the momentum the atoms carry across a surface, the
/// whole pressure of an ideal gas. The second is the push (W > 0) or pull
/// (W < 0) between atoms on the two sides of it. LAMMPS writes the same
/// sum (compute_pressure.cpp). With T from `temperature`, the first part
/// is (N - 1) T / V. The rho T of PHYSICS_REFERENCE is N T / V, larger
/// by 1 part in N - 1.
///
/// Label: derived. Pairs beyond the cutoff add nothing to W, so this is
/// the pressure of the cut and shifted model, not of the full potential.
pub fn pressure(velocities: &[Vec3], virial: f64, space: PeriodicBox) -> f64 {
    let l = space.side;
    (sum_v2(velocities) + virial) / (3.0 * l * l * l)
}

/// Scales every velocity by sqrt(target / T), so that the temperature
/// becomes `target`, in eps / kB. The total momentum is scaled by the same
/// factor, so a total momentum of zero stays zero.
///
/// Label: assumed. Rescaling is put in by hand and is not part of the
/// motion. It is used only to reach the wanted temperature before
/// measuring, never while measuring.
pub fn rescale(velocities: &mut [Vec3], target: f64) {
    let factor = (target / temperature(velocities)).sqrt();
    for v in velocities {
        *v = *v * factor;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice;
    use crate::pairs;

    // Three atoms with total momentum zero. Sum of v^2 = 5 + 2 + 5 = 12,
    // and 3N - 3 = 6, so T = 2. Every value is exact in an f64.
    const V: [Vec3; 3] = [
        Vec3::new(1.0, 2.0, 0.0),
        Vec3::new(-1.0, 0.0, 1.0),
        Vec3::new(0.0, -2.0, -1.0),
    ];

    #[test]
    fn temperature_of_three_atoms_is_2() {
        assert_eq!(temperature(&V), 2.0);
    }

    #[test]
    fn rescale_scales_every_velocity_by_the_same_factor() {
        // From T = 2 to T = 0.5 the factor is sqrt(0.25) = 0.5, exact.
        let mut v = V;
        rescale(&mut v, 0.5);
        for (a, b) in v.iter().zip(&V) {
            assert_eq!(*a, *b * 0.5);
        }
        assert_eq!(temperature(&v), 0.5);
    }

    #[test]
    fn pressure_of_atoms_too_far_apart_to_interact_is_the_ideal_gas_part() {
        // Box of side 8, no pair within the cutoff of 2.5, so W = 0 and
        // P = 12 / (3 * 512) = 1/128, exact.
        let space = PeriodicBox { side: 8.0 };
        let positions = [
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(5.0, 1.0, 1.0),
            Vec3::new(1.0, 5.0, 1.0),
        ];
        let mut forces = [Vec3::new(0.0, 0.0, 0.0); 3];
        let (_, w) = pairs::compute(&positions, space, 2.5, &mut forces);
        assert_eq!(w, 0.0);
        assert_eq!(pressure(&V, w, space), 1.0 / 128.0);
    }

    #[test]
    fn pressure_of_a_lattice_at_rest_is_minus_the_slope_of_energy_with_volume() {
        // At rest only the virial is left, and P must equal -dU/dV. Grow
        // and shrink the box of 108 lattice atoms at Rahman's density by H
        // in side and measure dU/dV = (U(L + H) - U(L - H)) / (V+ - V-).
        // The nearest shells to the cutoff are 0.16 sigma away, so U
        // changes smoothly with L. Measured miss: 4.6e-8, on a pressure of
        // -6.04. With H = 1e-3 it is 4.6e-6, the H^2 pattern of a two-sided
        // slope. The tolerance, 1e-6, is our own choice.
        const H: f64 = 1e-4;
        let side = 10.229 / 2.0;
        let energy = |l: f64| {
            let positions = lattice::fcc(3, l);
            let mut forces = vec![Vec3::new(0.0, 0.0, 0.0); positions.len()];
            pairs::compute(&positions, PeriodicBox { side: l }, 2.25, &mut forces)
        };
        let (_, w) = energy(side);
        let at_rest = [Vec3::new(0.0, 0.0, 0.0); 108];
        let p = pressure(&at_rest, w, PeriodicBox { side });
        let (big, small) = (side + H, side - H);
        let volume_change = big * big * big - small * small * small;
        let slope = (energy(big).0 - energy(small).0) / volume_change;
        let miss = (p + slope).abs();
        assert!(miss <= 1e-6, "P {p}, -dU/dV {}", -slope);
    }
}
