//! The transport run: a Lennard-Jones liquid made ready at a chosen
//! temperature, then left alone while its currents are recorded.
//!
//! Getting there. Atoms start on an fcc lattice at the density of
//! viscardy.rs, with velocities for twice the temperature, and run under
//! the Nose-Hoover thermostat (thermostat.rs) until the lattice is gone.
//! The thermostat then stays on while the total energy is averaged. Last,
//! every velocity is scaled so that the total energy equals that average,
//! and the thermostat comes off. A box at the average energy of the
//! canonical states at T moves, on its own, at an average temperature
//! close to T, and the motion from there conserves energy, which the
//! Green-Kubo integrals ask for: a thermostat would damp the currents.
//!
//! Measuring. The run is cut into blocks. In each the currents of flux.rs
//! are sampled every SAMPLE_EVERY steps into the correlations of
//! greenkubo.rs, integrated out to MAX_LAG samples. Each block gives one
//! viscosity, read at VISCOSITY_LAG, and one conductivity, read at
//! CONDUCTIVITY_LAG, and the spread of the blocks gives the error bars
//! (stats.rs). Where to read is our choice, from a trial of 24 blocks of
//! 864 atoms: the viscosity rose from 3.233 +- 0.049 at 2 tau to 3.275 +-
//! 0.057 at 3 tau and then stayed within its error bar, 3.29 +- 0.08 at
//! 5 tau and 3.33 +- 0.11 at 10. The conductivity reached 6.97 +- 0.09 by
//! 1 tau and then stayed, 6.96 +- 0.15 at 2 tau and 6.82 +- 0.35 at 10.
//! Reading later only adds noise.
//!
//! Labels: the thermostat and the scaling are assumed, used only to get
//! there; the time step is approximated (system.rs).

use crate::flux;
use crate::forces::Method;
use crate::greenkubo::GreenKubo;
use crate::lattice;
use crate::periodic::PeriodicBox;
use crate::rng::Rng;
use crate::system::System;
use crate::thermo;
use crate::thermostat::Thermostat;
use crate::velocities;
use crate::viscardy::{CUTOFF, DENSITY, TEMPERATURE};

/// The time step, in tau. Ours, as in M2 to M8; Viscardy et al. used 0.003.
pub const DT: f64 = 0.005;
/// Steps under the thermostat to melt, then to average the energy.
pub const MELT_STEPS: usize = 20_000;
pub const AVERAGE_STEPS: usize = 20_000;
/// The thermostat's period, in tau: 100 steps, as LAMMPS suggests.
pub const PERIOD: f64 = 0.5;
/// The neighbor list's skin, in sigma, as since M5.
pub const SKIN: f64 = 0.3;
/// Steps between samples of the currents: 0.025 tau.
pub const SAMPLE_EVERY: usize = 5;
/// The last lag of the integrals, 5 tau, and where the viscosity (3 tau)
/// and the conductivity (2 tau) are read, all in samples.
pub const MAX_LAG: usize = 200;
pub const VISCOSITY_LAG: usize = 120;
pub const CONDUCTIVITY_LAG: usize = 80;
/// Steps in a block: 250 tau.
pub const BLOCK_STEPS: usize = 50_000;

/// What one block measured.
pub struct Block {
    /// Average temperature, in eps / kB.
    pub temperature: f64,
    /// Viscosity and conductivity with the integrals stopped at each lag,
    /// 0 to MAX_LAG (units in greenkubo.rs).
    pub viscosity: Vec<f64>,
    pub conductivity: Vec<f64>,
    /// The viscosity at VISCOSITY_LAG from the off-diagonal and from the
    /// diagonal parts of the stress alone.
    pub parts: (f64, f64),
    /// Largest change of the total energy from the start of the block, as
    /// a fraction of it. Unit: none.
    pub energy_change: f64,
}

/// 4 `cells`^3 atoms at the state of viscardy.rs, ready to measure, with
/// forces from a neighbor list split across `threads` threads. `seed`
/// draws the starting velocities.
pub fn liquid(cells: usize, seed: u64, threads: usize) -> System {
    let atoms = 4 * cells * cells * cells;
    let space = PeriodicBox {
        side: (atoms as f64 / DENSITY).cbrt(),
    };
    let positions = lattice::fcc(cells, space.side);
    let mut rng = Rng::new(seed);
    let hot = velocities::maxwell_boltzmann(atoms, 2.0 * TEMPERATURE, &mut rng);
    let method = Method::listed(SKIN, threads);
    let mut system = System::with_forces(positions, hot, space, CUTOFF, method);
    let mut thermostat = Thermostat::new(TEMPERATURE, PERIOD);
    let mut energy = 0.0;
    for k in 0..MELT_STEPS + AVERAGE_STEPS {
        thermostat.half_step(&mut system.velocities, DT);
        system.step(DT);
        thermostat.half_step(&mut system.velocities, DT);
        if k >= MELT_STEPS {
            energy += system.total_energy();
        }
    }
    let target = energy / AVERAGE_STEPS as f64 - system.potential;
    assert!(target > 0.0, "no kinetic energy left to give");
    let scale = (target / system.kinetic_energy()).sqrt();
    for v in &mut system.velocities {
        *v = *v * scale;
    }
    system
}

/// Runs one block of BLOCK_STEPS steps on `system`, with the currents
/// measured on `threads` threads. The system must use a neighbor list.
pub fn block(system: &mut System, threads: usize) -> Block {
    let space = system.space;
    let volume = space.side * space.side * space.side;
    let h = DT * SAMPLE_EVERY as f64;
    let mut gk = GreenKubo::new(MAX_LAG, h, volume);
    let start = system.total_energy();
    let mut energy_change: f64 = 0.0;
    for k in 1..=BLOCK_STEPS {
        system.step(DT);
        energy_change = energy_change.max((system.total_energy() - start).abs());
        if k % SAMPLE_EVERY == 0 {
            let Method::Listed { list, .. } = &system.method else {
                panic!("the transport run needs a neighbor list");
            };
            let (r, v) = (&system.positions, &system.velocities);
            let currents = flux::measure(r, v, space, system.cutoff, list, threads);
            gk.add(&currents, thermo::temperature(v));
        }
    }
    let (off, diagonal) = gk.viscosity_parts();
    Block {
        temperature: gk.temperature(),
        viscosity: gk.viscosity(),
        conductivity: gk.conductivity(),
        parts: (off[VISCOSITY_LAG], diagonal[VISCOSITY_LAG]),
        energy_change: energy_change / start.abs(),
    }
}
