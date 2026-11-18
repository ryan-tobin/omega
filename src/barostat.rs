//! The Monte Carlo barostat: holds the atoms at a chosen pressure by
//! letting the box change size.
//!
//! Every so many steps the barostat proposes a new box side L', drawn
//! evenly within `max_step` of the old side L, with every position scaled
//! by L' / L. It works out the potential energy U' of the scaled atoms and
//!
//!   dW = U' - U + P (L'^3 - L^3) - (3 (N - 1) + 2) T ln(L' / L).
//!
//! The move is kept if dW <= 0, and otherwise with probability
//! exp(-dW / T). If not, the old box stays. Over many moves the box sides
//! then follow the weights L^(3 (N - 1) + 2) exp(-(U + P L^3) / T) of
//! atoms at constant pressure and temperature: a smaller box is paid for
//! by P V, a larger one leaves the atoms more room, and U says what the
//! atoms make of it. The 2 is there because the moves are even in L, not
//! in the volume. N - 1 and not N because the atoms as a whole do not
//! move (thermo.rs counts 3N - 3 velocity components for the same reason),
//! so the pressure of thermo.rs averages P, and atoms that never touch
//! fill a volume of N T / P on average, the ideal gas law (tested).
//!
//! The velocities are left alone: the thermostat holds the temperature.
//! A box smaller than two cutoffs is never accepted (pairs.rs).
//!
//! The method is that of K.-H. Chow and D. M. Ferguson, Comput. Phys.
//! Commun. 91, 283 (1995), and J. Aqvist and others, Chem. Phys. Lett.
//! 384, 288 (2004). OpenMM's MonteCarloBarostat uses it with moves even in
//! the volume and with N.
//!
//! Label: derived, from the statistical mechanics of constant pressure,
//! given exact motion between the moves.

use crate::forces::Method;
use crate::ln::ln;
use crate::periodic::PeriodicBox;
use crate::rng::Rng;
use crate::system::System;
use crate::vec3::Vec3;

/// A barostat and its state.
pub struct Barostat {
    /// The target pressure P, in eps / sigma^3.
    pub pressure: f64,
    /// The temperature T of the thermostat it works beside, in eps / kB.
    pub temperature: f64,
    /// The largest change of the box side in one move, in sigma.
    pub max_step: f64,
    /// Moves tried so far. Unit: none.
    pub tries: u64,
    /// Moves kept so far. Unit: none.
    pub accepted: u64,
    rng: Rng,
}

impl Barostat {
    /// A barostat for `pressure` at `temperature` whose moves change the
    /// box side by up to `max_step`, with random numbers from `seed`.
    /// `tries` and `accepted` are 0 for a new run. For a run that carries
    /// on from a checkpoint they are the saved counts: every try draws two
    /// random numbers, so the generator is moved on by 2 * tries.
    pub fn new(
        pressure: f64,
        temperature: f64,
        max_step: f64,
        seed: u64,
        tries: u64,
        accepted: u64,
    ) -> Self {
        assert!(temperature > 0.0 && max_step > 0.0);
        let mut rng = Rng::new(seed);
        for _ in 0..2 * tries {
            rng.next_u64();
        }
        Self {
            pressure,
            temperature,
            max_step,
            tries,
            accepted,
            rng,
        }
    }

    /// Tries one move of the box of `system`. `trial` computes the forces
    /// in the new box: a method that has not been used yet, since a
    /// neighbor list built for the old box does not fit the new one. If
    /// the move is kept, `system` gets the new box, positions, forces,
    /// potential energy and virial, and `trial` as its method. Returns
    /// true if the move was kept.
    pub fn attempt(&mut self, system: &mut System, mut trial: Method) -> bool {
        self.tries += 1;
        let step = self.max_step * (2.0 * self.rng.uniform() - 1.0);
        // In (0, 1], so its logarithm exists.
        let coin = 1.0 - self.rng.uniform();
        let old = system.space.side;
        let side = old + step;
        if side < 2.0 * system.cutoff {
            return false;
        }
        let scale = side / old;
        let positions: Vec<Vec3> = system.positions.iter().map(|&r| r * scale).collect();
        let space = PeriodicBox { side };
        let mut forces = vec![Vec3::new(0.0, 0.0, 0.0); positions.len()];
        let (potential, virial) = trial.compute(&positions, space, system.cutoff, &mut forces);
        let room = 3.0 * (positions.len() as f64 - 1.0) + 2.0;
        let work = self.pressure * (side * side * side - old * old * old);
        let change = potential - system.potential + work - self.temperature * room * ln(scale);
        if ln(coin) >= -change / self.temperature {
            return false;
        }
        self.accepted += 1;
        system.positions = positions;
        system.forces = forces;
        (system.potential, system.virial) = (potential, virial);
        system.space = space;
        system.method = trial;
        true
    }
}
