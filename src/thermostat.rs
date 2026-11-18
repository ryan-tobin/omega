//! The Nose-Hoover thermostat: holds the atoms at a chosen temperature
//! while they move.
//!
//! A friction xi is added to the motion: dv/dt = F - xi v. It is not fixed.
//! It grows while the atoms are hotter than the target and shrinks while
//! they are colder, at the rate dxi/dt = (T_now / T - 1) / period^2, where
//! T_now is the kinetic temperature (thermo.rs). A negative xi speeds the
//! atoms up. So the temperature swings around T and averages T, and energy
//! flows in and out as it does between a small sample and a large bath
//! around it. Positions and velocities then follow the canonical
//! distribution, the weights exp(-E / kB T) (S. Nose, J. Chem. Phys. 81,
//! 511 (1984); W. G. Hoover, Phys. Rev. A 31, 1695 (1985)).
//!
//! The total energy E now changes, but the motion still conserves
//!
//!   E + g T (period^2 xi^2 / 2 + eta),   deta/dt = xi,
//!
//! with g = 3N - 3 the number of free velocity components. `energy` returns
//! the second part, the energy the bath holds.
//!
//! In time the friction is applied for half a step before and half a step
//! after each velocity Verlet step, the order of G. J. Martyna, M. E.
//! Tuckerman, D. J. Tobias and M. L. Klein, Mol. Phys. 87, 1117 (1996), for
//! one thermostat. Friction held at xi for a time h multiplies every
//! velocity by exp(-xi h). We use (1 - a) / (1 + a) with a = xi h / 2
//! instead, the midpoint rule: it needs no exponential, differs from it by
//! 2 a^3 / 3, and a step with -xi still undoes a step with xi, so the
//! motion can be run backward.
//!
//! Label: approximated. One thermostat gives the canonical distribution
//! only if the motion wanders over every state of a given energy. A liquid
//! does. A few atoms on stiff springs do not.

use crate::thermo;
use crate::vec3::Vec3;

/// A thermostat and its state.
pub struct Thermostat {
    /// The target temperature T, in eps / kB.
    pub temperature: f64,
    /// How fast the friction answers, in tau. LAMMPS, which calls it
    /// Tdamp, suggests about 100 time steps.
    pub period: f64,
    /// The friction xi, in 1 / tau.
    pub friction: f64,
    /// eta, the integral of the friction over time. Unit: none.
    pub eta: f64,
}

impl Thermostat {
    /// A thermostat for `temperature`, in eps / kB, that answers over
    /// `period`, in tau. It starts with no friction.
    pub fn new(temperature: f64, period: f64) -> Self {
        assert!(temperature > 0.0 && period > 0.0);
        Self {
            temperature,
            period,
            friction: 0.0,
            eta: 0.0,
        }
    }

    /// Applies the friction to `velocities`, in sigma / tau, for half of a
    /// time step `dt`, in tau: a quarter-step change of xi, the friction on
    /// the velocities, then another quarter-step change of xi. Call it
    /// before and after each velocity Verlet step.
    pub fn half_step(&mut self, velocities: &mut [Vec3], dt: f64) {
        self.friction += 0.25 * dt * self.rate(velocities);
        let a = 0.25 * dt * self.friction;
        let scale = (1.0 - a) / (1.0 + a);
        for v in velocities.iter_mut() {
            *v = *v * scale;
        }
        self.eta += 0.5 * dt * self.friction;
        self.friction += 0.25 * dt * self.rate(velocities);
    }

    /// How fast the friction changes at these velocities, in 1 / tau^2.
    fn rate(&self, velocities: &[Vec3]) -> f64 {
        let ratio = thermo::temperature(velocities) / self.temperature;
        (ratio - 1.0) / (self.period * self.period)
    }

    /// The energy the bath holds, g T (period^2 xi^2 / 2 + eta), in eps,
    /// for `atoms` atoms. Added to the total energy of the atoms it gives
    /// the quantity the motion conserves.
    pub fn energy(&self, atoms: usize) -> f64 {
        let g = (3 * atoms - 3) as f64;
        let spring = 0.5 * self.period * self.period * self.friction * self.friction;
        g * self.temperature * (spring + self.eta)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Three atoms with total momentum zero, sum of v^2 = 12 and
    // temperature 12 / 6 = 2 (as in thermo.rs).
    const V: [Vec3; 3] = [
        Vec3::new(1.0, 2.0, 0.0),
        Vec3::new(-1.0, 0.0, 1.0),
        Vec3::new(0.0, -2.0, -1.0),
    ];

    #[test]
    fn half_step_matches_a_hand_calculation() {
        // T = 1, period = 1, dt = 1. The atoms are at 2, twice too hot, so
        // xi grows at (2 / 1 - 1) / 1^2 = 1 and after a quarter step is
        // 0.25. Then a = 0.25 * 0.25 = 0.0625, and every velocity is
        // multiplied by 0.9375 / 1.0625 = 15/17. eta = 0.5 * 0.25. The
        // atoms are now at 2 * (15/17)^2 = 450/289, so xi ends at
        // 0.25 + 0.25 * (450/289 - 1). Tolerance 1e-15 for rounding.
        let mut thermostat = Thermostat::new(1.0, 1.0);
        let mut v = V;
        thermostat.half_step(&mut v, 1.0);
        for (got, start) in v.iter().zip(&V) {
            assert_eq!(*got, *start * (15.0 / 17.0));
        }
        assert_eq!(thermostat.eta, 0.125);
        let want = 0.25 + 0.25 * (450.0 / 289.0 - 1.0);
        assert!((thermostat.friction - want).abs() <= 1e-15);
    }

    #[test]
    fn a_half_step_backward_undoes_a_half_step_forward() {
        // The same half step with -dt must bring back the velocities, xi
        // and eta. Tolerance 1e-15 for rounding.
        let mut thermostat = Thermostat::new(1.5, 0.5);
        (thermostat.friction, thermostat.eta) = (0.3, 0.1);
        let mut v = V;
        thermostat.half_step(&mut v, 0.01);
        assert!((v[0] - V[0]).norm() > 1e-4, "the step changed nothing");
        thermostat.half_step(&mut v, -0.01);
        for (got, start) in v.iter().zip(&V) {
            assert!((*got - *start).norm() <= 1e-15);
        }
        assert!((thermostat.friction - 0.3).abs() <= 1e-15);
        assert!((thermostat.eta - 0.1).abs() <= 1e-15);
    }

    #[test]
    fn bath_energy_matches_a_hand_calculation() {
        // 3 atoms: g = 6. T = 2, period = 0.5, xi = 2, eta = 0.25:
        // 6 * 2 * (0.25 * 4 / 2 + 0.25) = 9. Exact.
        let mut thermostat = Thermostat::new(2.0, 0.5);
        (thermostat.friction, thermostat.eta) = (2.0, 0.25);
        assert_eq!(thermostat.energy(3), 9.0);
    }
}
