//! Viscosity and thermal conductivity of the Lennard-Jones liquid at
//! T = 0.722, density 0.8442, cutoff 2.5 sigma, from the Green-Kubo
//! integrals (transport.rs), next to Viscardy, Servantie and Gaspard's
//! numbers for the same box (viscardy.rs).
//!
//! cargo run --release --bin transport -- [cells] [blocks]
//!
//! runs 4 cells^3 atoms (default 6: 864 atoms) for `blocks` blocks of
//! 250 tau (default 50, 2.5 million steps), seed 1, on every thread the
//! machine offers.

use omega::stats;
use omega::transport::{self, BLOCK_STEPS, DT, SAMPLE_EVERY};
use omega::transport::{CONDUCTIVITY_LAG as KAPPA_AT, MAX_LAG, VISCOSITY_LAG as ETA_AT};
use omega::units::{ARGON, BOLTZMANN_J_PER_K};
use omega::viscardy;
use std::time::Instant;

/// Where the integrals are shown stopped, in samples of 0.025 tau.
const SHOWN: [usize; 6] = [20, 40, 80, 120, 160, MAX_LAG];

fn main() {
    let args: Vec<usize> = std::env::args()
        .skip(1)
        .map(|a| a.parse().expect("a number"))
        .collect();
    let cells = args.first().copied().unwrap_or(6);
    let blocks = args.get(1).copied().unwrap_or(50);
    let threads = std::thread::available_parallelism().map_or(1, |n| n.get());
    let clock = Instant::now();
    let mut system = transport::liquid(cells, 1, threads);
    let atoms = system.positions.len();
    println!("{atoms} atoms, {blocks} blocks of {BLOCK_STEPS} steps, dt {DT}, {threads} threads");
    println!("block  temperature  viscosity  conductivity  energy change");
    let mut runs = Vec::new();
    for b in 0..blocks {
        let block = transport::block(&mut system, threads);
        let (eta, lambda) = (block.viscosity[ETA_AT], block.conductivity[KAPPA_AT]);
        let (t, de) = (block.temperature, block.energy_change);
        println!("{b:5}  {t:11.5}  {eta:9.4}  {lambda:12.4}  {de:13.2e}");
        runs.push(block);
    }
    let mean = |f: &dyn Fn(&transport::Block) -> f64| {
        stats::mean_and_error(&runs.iter().map(f).collect::<Vec<f64>>())
    };
    let (t, t_err) = mean(&|b| b.temperature);
    println!(
        "temperature {t:.5} +- {t_err:.5}, aimed at {}",
        viscardy::TEMPERATURE
    );
    println!("stopped at (tau)  viscosity         conductivity");
    for lag in SHOWN {
        let (eta, eta_err) = mean(&|b| b.viscosity[lag]);
        let (lambda, lambda_err) = mean(&|b| b.conductivity[lag]);
        let at = (lag * SAMPLE_EVERY) as f64 * DT;
        println!("{at:16.2}  {eta:.4} +- {eta_err:.4}  {lambda:.4} +- {lambda_err:.4}");
    }
    let (eta, eta_err) = mean(&|b| b.viscosity[ETA_AT]);
    let (lambda, lambda_err) = mean(&|b| b.conductivity[KAPPA_AT]);
    let (off, off_err) = mean(&|b| b.parts.0);
    let (diagonal, diagonal_err) = mean(&|b| b.parts.1);
    println!(
        "viscosity at 3 tau {eta:.4} +- {eta_err:.4}: off-diagonal parts {off:.4} +- {off_err:.4}, diagonal {diagonal:.4} +- {diagonal_err:.4}"
    );
    println!("conductivity at 2 tau {lambda:.4} +- {lambda_err:.4}");
    if let Some(row) = viscardy::row_for(atoms) {
        let off = |ours: f64, theirs: f64| 100.0 * (ours - theirs) / theirs;
        let (v, v_err) = row.viscosity;
        let (c, c_err) = row.conductivity;
        println!(
            "Viscardy et al. viscosity {v} +- {v_err}: ours {:+.1}%",
            off(eta, v)
        );
        println!(
            "Viscardy et al. conductivity {c} +- {c_err}: ours {:+.1}%",
            off(lambda, c)
        );
    }
    let own = |x: f64, e: f64| 100.0 * e / x;
    println!(
        "our error: viscosity {:.1}%, conductivity {:.1}%",
        own(eta, eta_err),
        own(lambda, lambda_err),
    );

    let (m, eps, sigma) = (ARGON.mass_kg, ARGON.eps_j, ARGON.sigma_m);
    let eta_unit = (m * eps).sqrt() / (sigma * sigma);
    let lambda_unit = BOLTZMANN_J_PER_K / (sigma * sigma) * (eps / m).sqrt();
    let (eta_si, lambda_si) = (1e3 * eta * eta_unit, lambda * lambda_unit);
    println!(
        "as argon at 86.6 K: viscosity {eta_si:.4} mPa s, conductivity {lambda_si:.4} W/(m K)"
    );
    println!("real argon there (roadmap, NIST WebBook): 0.295 mPa s, 0.136 W/(m K)");
    println!("took {:.0} s", clock.elapsed().as_secs_f64());
}
