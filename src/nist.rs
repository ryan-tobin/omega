//! The NIST benchmark for the Lennard-Jones liquid: 500 atoms at T = 0.85
//! and density 0.776, the potential cut at 3 sigma, with the long-range
//! corrections of tail.rs added to the energy and the pressure.
//!
//! NIST Standard Reference Simulation Website, Lennard-Jones fluid,
//! mmlapps.nist.gov/srs/LJ_PURE, read 2026-11-18:
//! - mc.htm: canonical Monte Carlo, N = 500, "3 sigma + standard long
//!   range corrections", 2.5e8 trials, pressure from the virial, error
//!   bars one standard deviation from blocks of 5e7 trials. At T = 0.85,
//!   density 0.776: energy per atom -5.5121 +- 0.000455, pressure
//!   0.0067714 +- 0.00177.
//! - md.htm: molecular dynamics at constant energy, velocity Verlet with
//!   dt = 0.005, N = 500, potential and force cut (not shifted) at 3
//!   sigma. Its row at density 0.776 has T = 0.851(1), energy -5.517(1)
//!   and pressure 0.030(6). It ran hotter than 0.85, and a check of the
//!   draft roadmap could not reproduce it, so it is shown, not tested.
//!
//! Our run holds T = 0.85 with the Nose-Hoover thermostat, which gives
//! canonical weights exp(-U / T), but for the shifted potential: our
//! forces are minus the slope of the shifted energy. The Monte Carlo
//! weighs states by the cut potential, which is lower by n U(r_c) when n
//! pairs are inside the cutoff, U(3) = -0.0055 eps each. So each of our
//! samples is also counted with the weight w = exp(-n U(r_c) / T), and
//! the weighted average is the Monte Carlo's average from our states:
//! sum of w A / sum of w, for any A (R. W. Zwanzig, J. Chem. Phys. 22,
//! 1420 (1954), cited from memory; the identity is one line of algebra).
//! The weights are taken relative to the block's first sample, which
//! cancels in the ratio and keeps the exponential in range.
//!
//! The pressure is rho T + W / (3V) + P_tail, with T the thermostat's: in
//! the canonical states the velocities carry exactly rho T, which is what
//! the Monte Carlo, which has no velocities, adds to its virial.

use crate::forces::Method;
use crate::lattice;
use crate::lj;
use crate::periodic::PeriodicBox;
use crate::rng::Rng;
use crate::system::System;
use crate::tail;
use crate::thermostat::Thermostat;
use crate::vec3::Vec3;
use crate::velocities;
use crate::{flux, thermo};

/// The state: temperature in eps / kB, number density in sigma^-3,
/// cutoff in sigma, and 4 * 5^3 = 500 atoms.
pub const TEMPERATURE: f64 = 0.85;
pub const DENSITY: f64 = 0.776;
pub const CUTOFF: f64 = 3.0;
pub const CELLS: usize = 5;
/// The Monte Carlo energy per atom (eps) and pressure (eps / sigma^3),
/// each with its error.
pub const MC_ENERGY: (f64, f64) = (-5.5121, 0.000455);
pub const MC_PRESSURE: (f64, f64) = (0.0067714, 0.00177);
/// The molecular dynamics row: temperature, energy and pressure.
pub const MD_ROW: [(f64, f64); 3] = [(0.851, 0.001), (-5.517, 0.001), (0.030, 0.006)];

/// Time step (tau), thermostat period (tau), steps to melt and settle,
/// steps per block and steps between samples. Our own choices.
pub const DT: f64 = 0.005;
pub const PERIOD: f64 = 0.5;
pub const SETTLE_STEPS: usize = 20_000;
pub const BLOCK_STEPS: usize = 20_000;
pub const SAMPLE_EVERY: usize = 10;

/// One block's averages: energy per atom (eps) and pressure
/// (eps / sigma^3), both with the corrections, as our run weighs its
/// states and as the cut potential weighs them, and the temperature.
pub struct Averages {
    pub energy: f64,
    pub pressure: f64,
    pub cut_energy: f64,
    pub cut_pressure: f64,
    pub temperature: f64,
}

/// Melts 500 atoms at the NIST state and runs `blocks` blocks of
/// BLOCK_STEPS steps under the thermostat, on `threads` threads.
pub fn run(seed: u64, blocks: usize, threads: usize) -> Vec<Averages> {
    let atoms = 4 * CELLS * CELLS * CELLS;
    let space = PeriodicBox {
        side: (atoms as f64 / DENSITY).cbrt(),
    };
    let volume = space.side * space.side * space.side;
    let step_height = lj::energy(Vec3::new(CUTOFF, 0.0, 0.0));
    let mut rng = Rng::new(seed);
    let hot = velocities::maxwell_boltzmann(atoms, 2.0 * TEMPERATURE, &mut rng);
    let positions = lattice::fcc(CELLS, space.side);
    let method = Method::listed(0.3, threads);
    let mut system = System::with_forces(positions, hot, space, CUTOFF, method);
    let mut thermostat = Thermostat::new(TEMPERATURE, PERIOD);
    let mut step = |system: &mut System| {
        thermostat.half_step(&mut system.velocities, DT);
        system.step(DT);
        thermostat.half_step(&mut system.velocities, DT);
    };
    for _ in 0..SETTLE_STEPS {
        step(&mut system);
    }
    let mut out = Vec::new();
    for _ in 0..blocks {
        let (mut u, mut p, mut t) = (0.0, 0.0, 0.0);
        let (mut w, mut wu, mut wp, mut first) = (0.0, 0.0, 0.0, None);
        for k in 1..=BLOCK_STEPS {
            step(&mut system);
            if k % SAMPLE_EVERY == 0 {
                let Method::Listed { list, .. } = &system.method else {
                    panic!("built with a list");
                };
                let (r, v) = (&system.positions, &system.velocities);
                let pairs = flux::measure(r, v, space, CUTOFF, list, threads).pairs;
                let u_k = tail::unshifted(system.potential, pairs, CUTOFF) / atoms as f64;
                let p_k = DENSITY * TEMPERATURE + system.virial / (3.0 * volume);
                let n0 = *first.get_or_insert(pairs) as f64;
                let w_k = (-(pairs as f64 - n0) * step_height / TEMPERATURE).exp();
                (u, p, t) = (u + u_k, p + p_k, t + thermo::temperature(v));
                (w, wu, wp) = (w + w_k, wu + w_k * u_k, wp + w_k * p_k);
            }
        }
        let samples = (BLOCK_STEPS / SAMPLE_EVERY) as f64;
        let (u_tail, p_tail) = (
            tail::energy_per_atom(DENSITY, CUTOFF),
            tail::pressure(DENSITY, CUTOFF),
        );
        out.push(Averages {
            energy: u / samples + u_tail,
            pressure: p / samples + p_tail,
            cut_energy: wu / w + u_tail,
            cut_pressure: wp / w + p_tail,
            temperature: t / samples,
        });
    }
    out
}
