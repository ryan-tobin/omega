//! N atoms in a repeating box, and the velocity Verlet step that moves them.
//!
//! The state is a struct of arrays: positions, velocities and forces each
//! live in their own flat array, with atom i at index i in all three.
//! Units are reduced (units.rs): positions in sigma, velocities in
//! sigma / tau, forces in eps / sigma, energies in eps, time in tau. Each
//! atom's mass is the mass unit, 1, so its acceleration equals its force.
//!
//! Positions are never wrapped back into the box. The nearest-image rule
//! in periodic.rs works for separations of any length, and unwrapped
//! positions keep the full path of each atom.

use crate::pairs;
use crate::periodic::PeriodicBox;
use crate::vec3::Vec3;

/// The atoms, their box and the force cutoff.
pub struct System {
    /// Position of each atom, in sigma.
    pub positions: Vec<Vec3>,
    /// Velocity of each atom, in sigma / tau.
    pub velocities: Vec<Vec3>,
    /// Total force on each atom at the current positions, in eps / sigma.
    pub forces: Vec<Vec3>,
    /// Potential energy at the current positions, in eps.
    pub potential: f64,
    /// Viral at the current positions, the sum over pairs of
    /// r_ij. F_ij (pairs.rs), in eps.
    pub virial: f64,
    /// The repeating box.
    pub space: PeriodicBox,
    /// Force cutoff r_c, in sigma.
    pub cutoff: f64,
}

impl System {
    /// Builds the system and computes the forces, potential energy and
    /// virial at the starting positions.
    pub fn new(
        positions: Vec<Vec3>,
        velocities: Vec<Vec3>,
        space: PeriodicBox,
        cutoff: f64,
    ) -> Self {
        assert_eq!(positions.len(), velocities.len());
        let mut forces = vec![Vec3::new(0.0, 0.0, 0.0); positions.len()];
        let (potential, virial) = pairs::compute(&positions, space, cutoff, &mut forces);
        Self {
            positions,
            velocities,
            forces,
            potential,
            virial,
            space,
            cutoff,
        }
    }

    /// Moves time forward by `dt`, in tau, with one velocity Verlet step:
    /// half a kick from the old forces, a drift at the new velocity, new
    /// forces, then half a kick from the new forces.
    ///
    /// Label: approximated. Over a fixed stretch of time the step follows
    /// the true motion with an error that shrinks as dt^2. It is
    /// time-reversible. For a smooth force its energy error wobbles without
    /// building up. The force here jumps at the cutoff, which adds a slow
    /// drift (see tests/conservation.rs).
    pub fn step(&mut self, dt: f64) {
        let half_dt = 0.5 * dt;
        for (v, &f) in self.velocities.iter_mut().zip(&self.forces) {
            *v = *v + f * half_dt;
        }
        for (r, &v) in self.positions.iter_mut().zip(&self.velocities) {
            *r = *r + v * dt;
        }
        (self.potential, self.virial) =
            pairs::compute(&self.positions, self.space, self.cutoff, &mut self.forces);
        for (v, &f) in self.velocities.iter_mut().zip(&self.forces) {
            *v = *v + f * half_dt;
        }
    }

    /// Kinetic energy, the sum of v^2 / 2 over all atoms, in eps.
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.velocities.iter().map(|v| v.norm_sq()).sum::<f64>()
    }

    /// Kinetic plus potential energy, in eps.
    pub fn total_energy(&self) -> f64 {
        self.kinetic_energy() + self.potential
    }

    /// Total momentum, the sum of the velocities (each mass is 1), in
    /// m * sigma / tau.
    pub fn momentum(&self) -> Vec3 {
        let mut p = Vec3::new(0.0, 0.0, 0.0);
        for &v in &self.velocities {
            p = p + v;
        }
        p
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two atoms 1.05 sigma apart, at rest, in a box of side 8. They swing
    /// between 1.05 and 1.26 sigma and never reach the cutoff at 2.5, so
    /// the force they feel is smooth. Returns the largest change in total
    /// energy, in eps, over `steps` steps of `dt`.
    fn two_atom_energy_error(dt: f64, steps: usize) -> f64 {
        let positions = vec![Vec3::new(3.475, 4.0, 4.0), Vec3::new(4.525, 4.0, 4.0)];
        let velocities = vec![Vec3::new(0.0, 0.0, 0.0); 2];
        let mut system = System::new(positions, velocities, PeriodicBox { side: 8.0 }, 2.5);
        let e0 = system.total_energy();
        let mut worst: f64 = 0.0;
        for _ in 0..steps {
            system.step(dt);
            worst = worst.max((system.total_energy() - e0).abs());
        }
        worst
    }

    #[test]
    fn energy_error_shrinks_four_times_when_dt_halves() {
        // Velocity Verlet is a second-order method, so its energy error
        // grows as dt^2: halving dt should cut it by a factor of 4. Both
        // runs cover 2 tau, almost three swings of 0.70 tau. Measured: 4.0008.
        // The tolerance, 0.05, is our own choice.
        let coarse = two_atom_energy_error(0.005, 400);
        let fine = two_atom_energy_error(0.0025, 800);
        let ratio = coarse / fine;
        assert!((ratio - 4.0).abs() <= 0.05, "ratio {ratio}");
    }

    #[test]
    fn momentum_is_the_sum_of_the_velocities() {
        // Two atoms 4 sigma apart, beyond the cutoff. Exact sums.
        let positions = vec![Vec3::new(1.0, 1.0, 1.0), Vec3::new(5.0, 1.0, 1.0)];
        let velocities = vec![Vec3::new(1.0, 2.0, 3.0), Vec3::new(0.5, 0.0, -1.0)];
        let system = System::new(positions, velocities, PeriodicBox { side: 8.0 }, 2.5);
        assert_eq!(system.momentum(), Vec3::new(1.5, 2.0, 2.0));
    }
}
