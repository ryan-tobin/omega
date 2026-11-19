//! The NIST benchmark (nist.rs): 500 atoms at T = 0.85 and density 0.776,
//! cut at 3 sigma, with long-range corrections, against NIST's Monte
//! Carlo. One run of 20 blocks of 100 tau is shared by the tests. About 5
//! minutes with --release on 4 cores, so a plain cargo test skips them:
//!
//! cargo test --release --test nist -- --ignored

use omega::nist::{self, Averages, MC_ENERGY, MC_PRESSURE};
use omega::stats;
use std::sync::OnceLock;

fn run() -> &'static [Averages] {
    static RUN: OnceLock<Vec<Averages>> = OnceLock::new();
    let threads = std::thread::available_parallelism().map_or(1, |n| n.get());
    RUN.get_or_init(|| nist::run(1, 20, threads))
}

/// The mean of `f` over the blocks and its standard error.
fn mean(f: fn(&Averages) -> f64) -> (f64, f64) {
    stats::mean_and_error(&run().iter().map(f).collect::<Vec<f64>>())
}

/// How many combined standard errors `ours` is from `theirs`.
fn apart(ours: (f64, f64), theirs: (f64, f64)) -> f64 {
    (ours.0 - theirs.0).abs() / (ours.1 * ours.1 + theirs.1 * theirs.1).sqrt()
}

#[test]
#[ignore]
fn the_thermostat_holds_the_temperature() {
    // Measured: 0.85000 +- 0.00004 over 20 blocks. Pass line 0.1%, our
    // own choice: the energy changes by about 0.002 per 0.001 of T.
    let (t, _) = mean(|a| a.temperature);
    assert!((t - 0.85).abs() <= 0.00085, "T = {t}");
}

#[test]
#[ignore]
fn energy_weighed_as_cut_matches_the_monte_carlo() {
    // Measured over 20 blocks: -5.51124 +- 0.00026, 1.6 combined standard
    // errors from -5.5121 +- 0.00046. Unweighted it is -5.51197 +- 0.00024,
    // 0.3 apart. Pass line 3, our own choice.
    let ours = mean(|a| a.cut_energy);
    let n = apart(ours, MC_ENERGY);
    assert!(
        n <= 3.0,
        "energy {ours:?}, {n:.1} standard errors from NIST"
    );
}

#[test]
#[ignore]
fn pressure_weighed_as_cut_matches_the_monte_carlo() {
    // Measured over 20 blocks: 0.0097 +- 0.0016, 1.2 combined standard
    // errors from 0.0068 +- 0.0018. Our own weighting moves it by +0.009:
    // unweighted it is 0.0009 +- 0.0015, 2.5 apart. Pass line 3, our own
    // choice.
    let ours = mean(|a| a.cut_pressure);
    let n = apart(ours, MC_PRESSURE);
    assert!(
        n <= 3.0,
        "pressure {ours:?}, {n:.1} standard errors from NIST"
    );
}
