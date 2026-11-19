//! What the cutoff leaves out: the long-range corrections to the energy
//! and the pressure, and the energy of the cut model without its shift.
//!
//! Pairs beyond the cutoff r_c add nothing in pairs.rs. Far apart, atoms
//! are close to evenly spread, g(r) = 1, so on average each atom has
//! rho 4 pi r^2 dr partners at a distance r to r + dr. Counting each pair
//! once, the energy left out per atom and the pressure left out are
//!
//!   U_tail / N = (rho / 2) integral from r_c to infinity of 4 pi r^2 U(r) dr
//!              = (8/3) pi rho (r_c^-9 / 3 - r_c^-3),
//!   P_tail     = -(2 pi / 3) rho^2 integral of r^3 dU/dr dr
//!              = (16/3) pi rho^2 (2 r_c^-9 / 3 - r_c^-3),
//!
//! for U = 4 (r^-12 - r^-6) (M. P. Allen and D. J. Tildesley, Computer
//! Simulation of Liquids, 2nd ed. (2017), cited from memory; the tests
//! check the closed forms against the integrals). The NIST benchmark
//! tables quote energies and pressures of the potential cut at r_c, not
//! shifted, with these corrections added. So `unshifted` puts the shift
//! back first.
//!
//! Label: approximated. g(r) beyond r_c is not exactly 1. Measured on
//! the same states of a NIST run (nist.rs, five blocks of 100 tau), the
//! energy per atom with its correction is -5.5236, -5.5148, -5.5169 and
//! -5.5163 for r_c = 2.5, 3, 3.5 and 4, where the corrections alone run
//! from -0.4155 to -0.1016.

use crate::lj;
use crate::vec3::Vec3;
use std::f64::consts::PI;

/// The energy per atom that pairs beyond `cutoff` (sigma) would add at
/// number density `density` (sigma^-3), in eps.
pub fn energy_per_atom(density: f64, cutoff: f64) -> f64 {
    let inv3 = 1.0 / (cutoff * cutoff * cutoff);
    8.0 / 3.0 * PI * density * (inv3 * inv3 * inv3 / 3.0 - inv3)
}

/// The pressure that pairs beyond `cutoff` (sigma) would add at number
/// density `density` (sigma^-3), in eps / sigma^3.
pub fn pressure(density: f64, cutoff: f64) -> f64 {
    let inv3 = 1.0 / (cutoff * cutoff * cutoff);
    16.0 / 3.0 * PI * density * density * (2.0 * inv3 * inv3 * inv3 / 3.0 - inv3)
}

/// The potential energy of the potential cut at `cutoff` (sigma) but not
/// shifted, in eps, from the shifted energy `potential` of pairs.rs (eps)
/// and the number of pairs inside the cutoff, `pairs` (flux.rs). Each of
/// those pairs was shifted up by -U(r_c).
pub fn unshifted(potential: f64, pairs: usize, cutoff: f64) -> f64 {
    potential + pairs as f64 * lj::energy(Vec3::new(cutoff, 0.0, 0.0))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats;

    /// The integral from `cutoff` to infinity of `f`, by the trapezoid
    /// rule in s = 1 / r, from 0 to 1 / cutoff in 20,000 steps: dr = -ds / s^2,
    /// so r^2 U(r) dr and r^3 dU/dr dr become polynomials in s.
    fn integral(f: impl Fn(f64) -> f64, cutoff: f64) -> f64 {
        let steps = 20_000;
        let h = 1.0 / cutoff / steps as f64;
        let g = |s: f64| if s == 0.0 { 0.0 } else { f(1.0 / s) / (s * s) };
        let samples: Vec<f64> = (0..=steps).map(|k| g(k as f64 * h)).collect();
        stats::trapezoid(&samples, h)
    }

    #[test]
    fn closed_forms_match_the_integrals() {
        // U(r) = 4 (r^-12 - r^-6), dU/dr = -24 (2 r^-13 - r^-7). In s the
        // integrands are 16 pi (s^8 - s^2) and 24 (2 s^8 - s^2), so the
        // trapezoid rule misses by about h^2 / 12 times the change of
        // slope. Measured: at most 1.25e-9 relative. Tolerance 1e-8, our own
        // choice.
        let u = |r: f64| 4.0 * (r.powi(-12) - r.powi(-6));
        let du = |r: f64| -24.0 * (2.0 * r.powi(-13) - r.powi(-7));
        for (density, cutoff) in [(0.8442, 2.5), (0.776, 3.0), (0.5, 4.0)] {
            let e = density / 2.0 * integral(|r| 4.0 * PI * r * r * u(r), cutoff);
            let p = -2.0 * PI / 3.0 * density * density * integral(|r| r.powi(3) * du(r), cutoff);
            let (got_e, got_p) = (energy_per_atom(density, cutoff), pressure(density, cutoff));
            assert!(((got_e - e) / e).abs() <= 1e-8, "energy {got_e} {e}");
            assert!(((got_p - p) / p).abs() <= 1e-8, "pressure {got_p} {p}");
        }
    }

    #[test]
    fn corrections_at_the_nist_state() {
        // rho 0.776, r_c 3: 3^-3 = 1/27 and 3^-9 = 1/19683, so
        // U_tail / N = (8/3) pi 0.776 (1/59049 - 1/27) = -0.240668 and
        // P_tail = (16/3) pi 0.776^2 (2/59049 - 1/27) = -0.373346, worked
        // out by hand to 6 digits. Tolerance 1e-6.
        assert!((energy_per_atom(0.776, 3.0) + 0.240668).abs() <= 1e-6);
        assert!((pressure(0.776, 3.0) + 0.373346).abs() <= 1e-6);
    }

    #[test]
    fn unshifting_adds_the_shift_back_for_each_pair() {
        // U(2.5) = 4 (2.5^-12 - 2.5^-6) = -0.0163169, so 10 pairs lower
        // the energy by 0.163169. Tolerance 1e-6.
        let got = unshifted(-1.0, 10, 2.5);
        assert!((got + 1.163169).abs() <= 1e-6, "{got}");
    }
}
