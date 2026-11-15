//! Prints the size in SI of each reduced unit for argon, the Lennard-Jones
//! energy and force between two atoms, the energy minimum in SI, a run of
//! 864 atoms in a repeating box the size of Rahman's, and how closely that
//! run retraces its path when every velocity is flipped.

use omega::lattice;
use omega::lj;
use omega::periodic::PeriodicBox;
use omega::rng::Rng;
use omega::system::System;
use omega::units::ARGON;
use omega::vec3::Vec3;
use omega::velocities;

// From Rahman (PHYSICS_REFERENCE), in reduced units: the box side, the
// cutoff, and his temperature, used only to draw the starting velocities.
// About half the starting kinetic energy then goes into the vibrations of
// the lattice, so the run ends cooler. The time step and the seed are ours.
const SIDE: f64 = 10.229;
const TEMPERATURE: f64 = 0.7867;
const CUTOFF: f64 = 2.25;
const DT: f64 = 0.005;
const SEED: u64 = 1964;

fn main() {
    println!("Argon reduced units in SI");
    println!("length       {:.4e} m", ARGON.sigma_m);
    println!("energy       {:.4e} J", ARGON.eps_j);
    println!("mass         {:.4e} kg", ARGON.mass_kg);
    println!("time         {:.4e} s", ARGON.time_s());
    println!("temperature  {:.4e} K", ARGON.temperature_k());
    println!("velocity     {:.4e} m/s", ARGON.velocity_m_per_s());
    println!("pressure     {:.4e} Pa", ARGON.pressure_pa());

    // Two atoms on the x axis, a distance r apart, in reduced units: r in
    // sigma, U in eps, F in eps / sigma. A positive force pushes the atoms
    // apart. A negative force pulls them together.
    println!();
    println!("Two atoms, Lennard-Jones, reduced units");
    println!("     r         U         F");
    for r in [0.95, 1.0, 1.1, 1.2, 1.5, 2.0, 2.5] {
        let r_ij = Vec3::new(r, 0.0, 0.0);
        let u = lj::energy(r_ij);
        let f = lj::force(r_ij).x;
        println!("{r:6.2} {u:9.4} {f:9.4}");
    }

    // The bottom of the energy well at r = 2^(1/6), converted to SI for
    // argon by multiplying by the length unit and the energy unit.
    let r_min = 2.0_f64.powf(1.0 / 6.0);
    let u_min = lj::energy(Vec3::new(r_min, 0.0, 0.0));
    println!();
    println!("Lennard-Jones energy minimum, in SI for argon");
    println!("separation   {:.4e} m", r_min * ARGON.sigma_m);
    println!("energy       {:.4e} J", u_min * ARGON.eps_j);

    run_many_atoms();
    run_forward_and_back();
}

/// 864 atoms (4 * 6^3) on an fcc lattice filling Rahman's box, with
/// velocities drawn for his temperature.
fn rahman_box() -> System {
    let space = PeriodicBox { side: SIDE };
    let positions = lattice::fcc(6, SIDE);
    let mut rng = Rng::new(SEED);
    let velocities = velocities::maxwell_boltzmann(positions.len(), TEMPERATURE, &mut rng);
    System::new(positions, velocities, space, CUTOFF)
}

/// Steps the box forward 10,000 times and prints the energies, in eps,
/// every 1000 steps. Time is in tau.
fn run_many_atoms() {
    let mut system = rahman_box();
    let e0 = system.total_energy();
    let mut worst: f64 = 0.0;
    println!();
    println!("864 atoms, box side {SIDE}, cutoff {CUTOFF}, dt {DT}, seed {SEED}");
    println!(" step   time    kinetic  potential      total");
    for step in 0..=10_000 {
        if step > 0 {
            system.step(DT);
        }
        let e = system.total_energy();
        worst = worst.max((e - e0).abs() / e0.abs());
        if step % 1000 == 0 {
            let time = step as f64 * DT;
            let kinetic = system.kinetic_energy();
            let potential = system.potential;
            println!("{step:5} {time:6.1} {kinetic:10.4} {potential:10.4} {e:10.4}");
        }
    }
    println!("largest change in total energy, relative  {worst:.3e}");
    println!(
        "total momentum at the end                 {:.3e}",
        system.momentum().norm()
    );
}

/// Runs the box forward, flips every velocity, runs it the same number of
/// steps again and prints how far the farthest atom ends from its start.
fn run_forward_and_back() {
    println!();
    println!("Forward, then backward: farthest atom from its start, in sigma");
    for steps in [1000, 2000] {
        let mut system = rahman_box();
        let start = system.positions.clone();
        for _ in 0..steps {
            system.step(DT);
        }
        for v in &mut system.velocities {
            *v = -*v;
        }
        for _ in 0..steps {
            system.step(DT);
        }
        let mut farthest: f64 = 0.0;
        for (r, &r0) in system.positions.iter().zip(&start) {
            farthest = farthest.max((*r - r0).norm());
        }
        println!("{steps:5} steps each way  {farthest:.3e}");
    }
}
