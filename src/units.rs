//! Reduced units for layer 1, and how big each one is in SI.
//!
//! Inside the engine, lengths are counted in sigma, energies in eps and
//! masses in the mass of one atom. Every other units follows those
//! three. To turn a reduced number into SI, multiple it by the unit:
//! a time of 3.0 reduced is 3.0 * time_s() seconds.

/// Boltzmann constant in J/K.
/// Source: NIST SP 961, CODATA recommended values 2022.
pub const BOLTZMANN_J_PER_K: f64 = 1.380649e-23;

// Atomic mass constant (1 u) in kg. Measured 1.660 539 068 92(52) e-27
/// Source: NIST SP 961, CODATA recommended values 2022.
pub const ATOMIC_MASS_KG: f64 = 1.66053906892e-27;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Units {
    /// Length unit sigma, in meters.
    pub sigma_m: f64,
    /// Energy unit eps, in joules.
    pub eps_j: f64,
    /// Mass unit, the mass of one atom, in kilograms.
    pub mass_kg: f64,
}

/// Argon as modeled by Rahman, Phys. Rev. 136, A405 (1964):
/// sigma = 3.4 angstrom and eps/kB = 120 K. THe mass is 39.95 u, the
/// CIAAW abridged standard atomic weight of argon.
/// Label: assumed. sigma and eps are fitted values put in by hand
pub const ARGON: Units = Units {
    sigma_m: 3.4e-10,
    eps_j: 120.0 * BOLTZMANN_J_PER_K,
    mass_kg: 39.95 * ATOMIC_MASS_KG,
};

impl Units {
    /// Time unit tau = sigma * sqrt(m / eps), in seconds.
    pub fn time_s(&self) -> f64 {
        self.sigma_m * (self.mass_kg / self.eps_j).sqrt()
    }

    /// Temperature unit eps / kB, in kelvin.
    pub fn temperature_k(&self) -> f64 {
        self.eps_j / BOLTZMANN_J_PER_K
    }

    /// Velocity unit sigma / tau = sqrt(eps / m), in meters per second.
    pub fn velocity_m_per_s(&self) -> f64 {
        (self.eps_j / self.mass_kg).sqrt()
    }

    /// Pressure unit eps / sigma^3, in pascals.
    pub fn pressure_pa(&self) -> f64 {
        self.eps_j / self.sigma_m.powi(3)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Every expected value is a row of the reduced-units table in
    // PHYSICS_REFERENCE. The table prints four or five significant figures,
    // so each tolerance is half a unit in the last printed digit. The
    // tolerance is our own choice: the tightest the table supports.
    fn assert_close(got: f64, want: f64, tol: f64) {
        let miss = (got - want).abs();
        assert!(miss <= tol, "got {got:e}, want {want:e}, off by {miss:e}");
    }

    #[test]
    fn argon_energy_unit_in_joules() {
        assert_close(ARGON.eps_j, 1.6568e-21, 0.5e-25);
    }

    #[test]
    fn argon_mass_unit_in_kilograms() {
        assert_close(ARGON.mass_kg, 6.6339e-26, 0.5e-30);
    }

    #[test]
    fn argon_time_unit_is_2_151_picoseconds() {
        assert_close(ARGON.time_s(), 2.151e-12, 0.5e-15);
    }

    #[test]
    fn argon_temperature_unit_is_120_kelvin() {
        assert_close(ARGON.temperature_k(), 120.0, 0.5);
    }

    #[test]
    fn argon_velocity_unit_is_158_0_meters_per_second() {
        assert_close(ARGON.velocity_m_per_s(), 158.0, 0.05);
    }

    #[test]
    fn argon_pressure_unit_is_42_15_megapascals() {
        assert_close(ARGON.pressure_pa(), 42.15e6, 0.005e6);
    }
}
