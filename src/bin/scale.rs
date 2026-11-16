//! Times the force methods of forces.rs from 864 atoms to about a million:
//! all pairs, and a neighbor list on 1 thread and on every thread the
//! machine offers.
//!
//! Each size is an fcc lattice at Rahman's density with velocities for
//! his temperature, his cutoff and our time step (rahman.rs), stepped
//! STEPS times by velocity Verlet. The list is built again about every
//! 10 steps (the builds column counts every build, the first included).
//! A time is the wall-clock time of those steps divided by STEPS, after
//! the system is built. All pairs runs only up to ALL_PAIRS_UP_TO atoms,
//! since its work grows as N^2.
//!
//! After the steps, the program prints how far apart the runs ended up:
//! the list runs on 1 and on many threads must agree exactly, and all
//! pairs and the list to rounding.

use omega::forces::Method;
use omega::lattice;
use omega::periodic::PeriodicBox;
use omega::rahman::{CUTOFF, DT, SIDE, TEMPERATURE};
use omega::rng::Rng;
use omega::system::System;
use omega::velocities;
use std::thread;
use std::time::Instant;

// Our own choices: the seed as in M2 to M4, LAMMPS's default skin in
// reduced units (neighbors.rs), 100 steps, and fcc lattices of 6 to 63
// cells a side, 864 to 1,000,188 atoms, about 4 times more at each size.
const SEED: u64 = 1964;
const SKIN: f64 = 0.3;
const STEPS: usize = 100;
const CELLS: [usize; 6] = [6, 10, 16, 25, 40, 63];
const ALL_PAIRS_UP_TO: usize = 4000;

/// 4 * cells^3 atoms on an fcc lattice at Rahman's density, with
/// velocities for his temperature, and forces computed by `method`.
fn argon(cells: usize, method: Method) -> System {
    let space = PeriodicBox {
        side: SIDE / 6.0 * cells as f64,
    };
    let positions = lattice::fcc(cells, space.side);
    let mut rng = Rng::new(SEED);
    let velocities = velocities::maxwell_boltzmann(positions.len(), TEMPERATURE, &mut rng);
    System::with_forces(positions, velocities, space, CUTOFF, method)
}

/// Steps `system` STEPS times. Returns the wall-clock time per step, in
/// milliseconds.
fn time(system: &mut System) -> f64 {
    let start = Instant::now();
    for _ in 0..STEPS {
        system.step(DT);
    }
    start.elapsed().as_secs_f64() * 1e3 / STEPS as f64
}

/// The farthest any atom of `a` is from the same atom of `b`, in sigma.
fn apart(a: &System, b: &System) -> f64 {
    let pairs = a.positions.iter().zip(&b.positions);
    pairs.map(|(&p, &q)| (p - q).norm()).fold(0.0, f64::max)
}

fn main() {
    let threads = thread::available_parallelism().map_or(1, |n| n.get());
    println!("fcc lattice at Rahman's density, box side {SIDE} per 864 atoms");
    println!("cutoff {CUTOFF}, skin {SKIN}, dt {DT}, seed {SEED}, {STEPS} steps");
    println!("{threads} threads available");
    println!();
    println!("milliseconds per step");
    println!("    atoms  all pairs   list, 1 thread   list, {threads} threads  builds");
    let mut checks = Vec::new();
    for cells in CELLS {
        let mut one = argon(cells, Method::listed(SKIN, 1));
        let n = one.positions.len();
        let t_one = time(&mut one);
        let mut many = argon(cells, Method::listed(SKIN, threads));
        let t_many = time(&mut many);
        let builds = many.method.builds();
        let mut all = String::from("         -");
        let mut gap = None;
        if n <= ALL_PAIRS_UP_TO {
            let mut system = argon(cells, Method::AllPairs);
            all = format!("{:10.3}", time(&mut system));
            gap = Some(apart(&system, &many));
        }
        println!("{n:9} {all} {t_one:16.3} {t_many:16.3} {builds:7}");
        checks.push((n, apart(&one, &many), gap));
    }
    println!();
    println!("farthest apart after {STEPS} steps, in sigma");
    println!("    atoms  list, 1 and {threads} threads  all pairs and list");
    for (n, threads_gap, all_gap) in checks {
        let all_gap = all_gap.map_or(String::from("-"), |g| format!("{g:.1e}"));
        println!("{n:9} {threads_gap:24.1e} {all_gap:>19}");
    }
}