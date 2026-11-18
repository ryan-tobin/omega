//! Liquid argon under a thermostat, shared by tests/thermostat.rs,
//! tests/barostat.rs and tests/restart.rs.

use omega::lattice;
use omega::periodic::PeriodicBox;
use omega::rng::Rng;
use omega::system::System;
use omega::thermostat::Thermostat;
use omega::velocities;

/// Rahman's temperature, in eps / kB.
pub const T: f64 = 0.7867;
/// The time step, in tau.
const DT: f64 = 0.005;

/// One time step with the thermostat on.
pub fn step(system: &mut System, thermostat: &mut Thermostat) {
    thermostat.half_step(&mut system.velocities, DT);
    system.step(DT);
    thermostat.half_step(&mut system.velocities, DT);
}

/// 108 atoms at Rahman's density with his cutoff of 2.25 sigma, as in
/// tests/conservation.rs, and a thermostat set to his temperature with a
/// period of 0.5 tau, 100 time steps. The atoms start on the lattice at
/// twice his temperature, which melts it, and settle for 4000 steps.
pub fn liquid(seed: u64) -> (System, Thermostat) {
    let space = PeriodicBox { side: 10.229 / 2.0 };
    let positions = lattice::fcc(3, space.side);
    let mut rng = Rng::new(seed);
    let start = velocities::maxwell_boltzmann(positions.len(), 2.0 * T, &mut rng);
    let mut system = System::new(positions, start, space, 2.25);
    let mut thermostat = Thermostat::new(T, 0.5);
    for _ in 0..4000 {
        step(&mut system, &mut thermostat);
    }
    (system, thermostat)
}
