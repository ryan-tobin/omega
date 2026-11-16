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
//!
//! The forces come from all pairs, or from a neighbor list split across
//! threads (forces.rs).

use crate::forces::Method;
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
    /// Virial at the current positions, the sum over pairs of
    /// r_ij . F_ij (pairs.rs), in eps.
    pub virial: f64,
    /// The repeating box.
    pub space: PeriodicBox,
    /// Force cutoff r_c, in sigma.
    pub cutoff: f64,
    /// How the forces are computed.
    pub method: Method,
}

impl System {
    /// Builds the system with forces from all pairs, and computes the
    /// forces, potential energy and virial at the starting positions.
    pub fn new(
        positions: Vec<Vec3>,
        velocities: Vec<Vec3>,
        space: PeriodicBox,
        cutoff: f64,
    ) -> Self {
        Self::with_forces(positions, velocities, space, cutoff, Method::AllPairs)
    }

    /// The same, with the forces computed by `method`.
    pub fn with_forces(
        positions: Vec<Vec3>,
        velocities: Vec<Vec3>,
        space: PeriodicBox,
        cutoff: f64,
        mut method: Method,
    ) -> Self {
        assert_eq!(positions.len(), velocities.len());
        let mut forces = vec![Vec3::new(0.0, 0.0, 0.0); positions.len()];
        let (potential, virial) = method.compute(&positions, space, cutoff, &mut forces);
        Self {
            positions,
            velocities,
            forces,
            potential,
            virial,
            space,
            cutoff,
            method,
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
            self.method
                .compute(&self.positions, self.space, self.cutoff, &mut self.forces);
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