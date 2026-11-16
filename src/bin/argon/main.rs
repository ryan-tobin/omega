//! Rahman's liquid argon, measured and compared with his paper: A. Rahman,
//! Phys. Rev. 136, A405 (1964).
//!
//! Melts 864 atoms in his box at his temperature (rahman.rs), lets them
//! settle, then measures 10 blocks of 20 tau of free motion (run.rs) and
//! prints our numbers next to his (compare.rs). Writes the measured curves
//! as CSV in the folder `output/argon` (files.rs), and the plots of
//! docs/liquid-argon.md as SVG in the folder `docs` (plots.rs).

mod compare;
mod files;
mod plots;

use omega::rahman::{self, CUTOFF, DT, RESCALE_STEPS, SIDE, TEMPERATURE};
use omega::run::{self, Run};
use omega::structure;
use omega::units::ARGON;
use std::io;

// Our own choices: the seed, as in M2 and M3, 2000 free steps (10 tau) for
// the liquid to settle after the last rescaling, and 10 blocks.
const SEED: u64 = 1964;
const SETTLE_STEPS: usize = 2000;
const BLOCKS: usize = 10;

/// The measured curves, as points (x, y) in reduced units: g against r
/// (sigma), S against k (1 / sigma), and the MSD (sigma^2), the VACF
/// ((sigma / tau)^2) and VACF / VACF(0) against t (tau), averaged over
/// the blocks.
pub struct Curves {
    pub g: Vec<(f64, f64)>,
    pub s: Vec<(f64, f64)>,
    pub msd: Vec<(f64, f64)>,
    pub vacf: Vec<(f64, f64)>,
    pub normal: Vec<(f64, f64)>,
}

impl Curves {
    /// The curves of `run`. S(k) is computed at k = 0.05, 0.10, ... 30.
    fn of(run: &Run) -> Self {
        let h = run::SAMPLE_EVERY as f64 * DT;
        let at = |values: Vec<f64>| -> Vec<(f64, f64)> {
            let points = values.into_iter().enumerate();
            points.map(|(k, v)| (k as f64 * h, v)).collect()
        };
        let g = run.rdf.g();
        let density = 864.0 / (SIDE * SIDE * SIDE);
        let s = (1..=600).map(|i| 0.05 * f64::from(i));
        let s = s
            .map(|k| (k, structure::structure_factor(&g, density, k)))
            .collect();
        let vacf = at(run.curve(|b| &b.vacf));
        let normal = vacf.iter().map(|&(t, v)| (t, v / vacf[0].1)).collect();
        let msd = at(run.curve(|b| &b.msd));
        Self {
            g,
            s,
            msd,
            vacf,
            normal,
        }
    }
}

/// Angstrom per sigma, picoseconds per tau, and 1e-5 cm^2 / s per
/// sigma^2 / tau, for argon: the units Rahman wrote in.
pub fn si() -> (f64, f64, f64) {
    let cm = ARGON.sigma_m * 100.0;
    let d = 1e5 * cm * cm / ARGON.time_s();
    (ARGON.sigma_m * 1e10, ARGON.time_s() * 1e12, d)
}

fn main() -> io::Result<()> {
    let mut system = rahman::liquid(SEED);
    for _ in 0..SETTLE_STEPS {
        system.step(DT);
    }
    let (steps, every) = (run::BLOCK_STEPS, run::SAMPLE_EVERY);
    println!("864 atoms, box side {SIDE}, cutoff {CUTOFF}, dt {DT}, seed {SEED}");
    println!("melted at T = {TEMPERATURE} for {RESCALE_STEPS} steps, then {SETTLE_STEPS} free");
    println!("measured in {BLOCKS} blocks of {steps} free steps, sampled every {every}");
    println!();
    println!("block       T    D (MSD)   D (VACF)");
    let mut run = Run::new(&system);
    for i in 1..=BLOCKS {
        let b = run.add_block(&mut system, DT);
        let (t, dm, dv) = (b.temperature, b.d_msd, b.d_vacf);
        println!("{i:5} {t:7.4} {dm:10.5} {dv:10.5}");
    }
    let change = run.energy_change;
    println!("largest relative change in total energy  {change:.3e}");
    let curves = Curves::of(&run);
    compare::print(&run, &curves);
    files::write(&run, &curves)?;
    plots::write(&run, &curves)
}
