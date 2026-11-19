//! Viscosity and thermal conductivity of the Lennard-Jones liquid at
//! T = 0.722, density 0.8442, cutoff 2.5 (transport.rs), against
//! Viscardy, Servantie and Gaspard (viscardy.rs). Two runs of 50 blocks
//! of 250 tau, 2.5 million steps each, 864 and 108 atoms, shared by the
//! tests. Every pass line is our own choice. They took 48 minutes with
//! --release on 4 cores, so a plain cargo test skips them:
//!
//! cargo test --release --test transport -- --ignored

use omega::stats;
use omega::transport::{self, Block, CONDUCTIVITY_LAG, VISCOSITY_LAG};
use omega::viscardy::{self, TEMPERATURE};
use std::sync::OnceLock;

/// Means and standard errors over the blocks of one run.
struct Measured {
    temperature: (f64, f64),
    viscosity: (f64, f64),
    conductivity: (f64, f64),
    off: (f64, f64),
    diagonal: (f64, f64),
}

fn summarize(cells: usize) -> Measured {
    let threads = std::thread::available_parallelism().map_or(1, |n| n.get());
    let mut system = transport::liquid(cells, 1, threads);
    let blocks: Vec<Block> = (0..50)
        .map(|_| transport::block(&mut system, threads))
        .collect();
    let mean = |f: &dyn Fn(&Block) -> f64| {
        stats::mean_and_error(&blocks.iter().map(f).collect::<Vec<f64>>())
    };
    Measured {
        temperature: mean(&|b| b.temperature),
        viscosity: mean(&|b| b.viscosity[VISCOSITY_LAG]),
        conductivity: mean(&|b| b.conductivity[CONDUCTIVITY_LAG]),
        off: mean(&|b| b.parts.0),
        diagonal: mean(&|b| b.parts.1),
    }
}

fn big() -> &'static Measured {
    static RUN: OnceLock<Measured> = OnceLock::new();
    RUN.get_or_init(|| summarize(6))
}

fn small() -> &'static Measured {
    static RUN: OnceLock<Measured> = OnceLock::new();
    RUN.get_or_init(|| summarize(3))
}

/// Checks that `ours` is within 5% of `theirs` and that its own error is
/// at most 2% of it: the roadmap's lines for M9.
fn close(name: &str, ours: (f64, f64), theirs: f64) {
    let off = (ours.0 - theirs).abs() / theirs;
    assert!(off <= 0.05, "{name} {ours:?} is {off:.3} from {theirs}");
    assert!(ours.1 <= 0.02 * ours.0, "{name} {ours:?}: error over 2%");
}

#[test]
#[ignore]
fn free_runs_stay_at_the_published_temperature() {
    // The energy set from the thermostat's average gives a free run near
    // 0.722. Measured: 864 atoms 0.72274 +- 0.00004, 108 atoms 0.72297 +-
    // 0.00016. Pass line 0.5% (the roadmap's).
    for run in [big(), small()] {
        let (t, _) = run.temperature;
        assert!((t - TEMPERATURE).abs() <= 0.005 * TEMPERATURE, "T = {t}");
    }
}

#[test]
#[ignore]
fn viscosity_of_864_atoms_matches_viscardy() {
    // Theirs: 3.247 +- 0.033. Measured: 3.258 +- 0.040, +0.3%.
    let theirs = viscardy::row_for(864).expect("in the table").viscosity.0;
    close("viscosity", big().viscosity, theirs);
}

#[test]
#[ignore]
fn conductivity_of_864_atoms_matches_viscardy() {
    // Theirs: 6.948 +- 0.07. Measured: 6.990 +- 0.089, +0.6%.
    let theirs = viscardy::row_for(864).expect("in the table").conductivity.0;
    close("conductivity", big().conductivity, theirs);
}

#[test]
#[ignore]
fn off_diagonal_and_diagonal_stress_give_one_viscosity() {
    // The two halves of greenkubo.rs's average must agree if the box has
    // no favored direction. They come from the same run, so their errors
    // are not independent; 3 of the larger error is used. Measured:
    // 3.215 +- 0.053 and 3.323 +- 0.069, 0.108 apart against 0.206.
    let s = big();
    let gap = (s.off.0 - s.diagonal.0).abs();
    assert!(
        gap <= 3.0 * s.off.1.max(s.diagonal.1),
        "{:?} {:?}",
        s.off,
        s.diagonal
    );
}

#[test]
#[ignore]
fn the_small_box_is_less_viscous_as_published() {
    // Theirs for 108 atoms: 3.057 +- 0.045 and 6.663 +- 0.10, below 864.
    // Measured: 3.147 +- 0.041 (+2.9%) and 6.638 +- 0.108 (-0.4%).
    let row = viscardy::row_for(108).expect("in the table");
    let (s, b) = (small(), big());
    close("viscosity", s.viscosity, row.viscosity.0);
    close("conductivity", s.conductivity, row.conductivity.0);
    assert!(
        s.viscosity.0 < b.viscosity.0,
        "{:?} {:?}",
        s.viscosity,
        b.viscosity
    );
}
