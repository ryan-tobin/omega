//! Boiling and melting of argon, against published results for the same
//! model: the Lennard-Jones potential cut and shifted at 2.5 sigma, which
//! is what pairs.rs computes with a cutoff of 2.5. Each test runs one of
//! the files in scenarios/ and takes minutes, so a plain cargo test skips
//! them: cargo test --release --test phases -- --ignored
//!
//! In reduced units the published numbers are:
//!
//! Liquid beside its vapor at T = 0.8: densities 0.730 and 0.0198 (J.
//! Vrabec, G. K. Kedia, G. Fuchs and H. Hasse, Mol. Phys. 104, 1509 (2006),
//! as quoted by M. Horsch and others, arXiv:1110.4466, Table 1), and
//! 0.730(1) and 0.020(1) from a run of the same kind as ours (S. W. Sides,
//! G. S. Grest and M.-D. Lacasse, arXiv:cond-mat/9907177, Table II).
//!
//! Crystal beside its melt at T = 0.7: pressure 0.928, with a volume per
//! atom of 1.0452 in the crystal and 1.1714 in the liquid, so densities
//! 0.9568 and 0.8537, and at T = 0.8: pressure 2.185 (U. R. Pedersen,
//! arXiv:1303.5165, Table I).
//!
//! Every tolerance is our own choice and says why.

use omega::drive::Run;
use omega::profile;
use omega::record::{self, SLICES};
use omega::scenario;
use std::path::Path;

/// Runs the scenario in `text` in a folder of its own and returns the rows
/// of `file` written by its last stage.
fn last_stage(name: &str, text: &str, file: &str) -> Vec<Vec<f64>> {
    let s = scenario::parse(text).unwrap();
    let folder = Path::new(env!("CARGO_TARGET_TMPDIR")).join(name);
    let threads = std::thread::available_parallelism().map_or(1, |n| n.get());
    let mut run = Run::new(&s, &folder, threads, false).unwrap();
    for index in 0..s.stages.len() {
        run.stage(index, None).unwrap();
    }
    let last = run.state.step;
    let first = last - s.stages.last().unwrap().steps as u64;
    record::rows(&folder.join(file), first, last).unwrap()
}

#[test]
#[ignore = "about 3 minutes with --release"]
fn liquid_and_vapor_densities_at_0_8_match_the_published_ones() {
    let rows = last_stage("boil", include_str!("../scenarios/boil.txt"), "profile.csv");
    let mean = |k: usize| rows.iter().map(|r| r[k]).sum::<f64>() / rows.len() as f64;
    let slices: Vec<f64> = (1..=SLICES).map(mean).collect();
    let (liquid, vapor) = profile::plateaus(&slices);
    println!("liquid {liquid:.4} (published 0.730), vapor {vapor:.4} (published 0.0198)");
    assert!(
        (liquid / 0.730 - 1.0).abs() <= 0.01,
        "liquid density {liquid}"
    );
    assert!(
        (vapor / 0.0198 - 1.0).abs() <= 0.15,
        "vapor density {vapor}"
    );
}

#[test]
#[ignore = "about 20 seconds with --release"]
fn crystal_at_the_melting_point_has_the_published_density() {
    let rows = last_stage(
        "solid",
        include_str!("../scenarios/solid.txt"),
        "thermo.csv",
    );
    let (density, _) = record::average(&rows, 5);
    println!("crystal {density:.4} (published 0.9568)");
    assert!(
        (density / 0.9568 - 1.0).abs() <= 0.005,
        "crystal density {density}"
    );
}

#[test]
#[ignore = "about 20 seconds with --release"]
fn liquid_at_the_melting_point_has_the_published_density() {
    let rows = last_stage(
        "liquid",
        include_str!("../scenarios/liquid.txt"),
        "thermo.csv",
    );
    let (density, _) = record::average(&rows, 5);
    println!("liquid {density:.4} (published 0.8537)");
    assert!(
        (density / 0.8537 - 1.0).abs() <= 0.01,
        "liquid density {density}"
    );
}

#[test]
#[ignore = "about 7 minutes with --release"]
fn crystal_beside_its_melt_at_0_7_settles_at_the_published_pressure() {
    // In a box that does not change, a pressure above the melting pressure
    // makes the crystal grow, which lowers the pressure, and the other way
    // round: the box settles at the melting pressure. Measured 0.939 with
    // this seed and 0.930 with seed 1, each with a standard error of 0.02
    // from five blocks. Tolerance: 0.07, which along the published melting
    // line (2.185 at T = 0.8) is under 1% of the melting temperature. Both
    // phases must still be there: the energy per atom must put the
    // crystal's share between 0.3 and 0.8, taking -6.20 for the crystal
    // and -5.28 for the liquid, as scenarios/solid.txt and liquid.txt
    // measure them. Measured: 0.55.
    let rows = last_stage("melt", include_str!("../scenarios/melt.txt"), "thermo.csv");
    let ((pressure, error), (energy, _)) = (record::average(&rows, 3), record::average(&rows, 4));
    let crystal = (-5.28 - energy) / (6.20 - 5.28);
    println!("pressure {pressure:.4} +- {error:.4} (published 0.928), {crystal:.2} crystal");
    assert!((pressure - 0.928).abs() <= 0.07, "pressure {pressure}");
    assert!(crystal > 0.3 && crystal < 0.8, "crystal share {crystal}");
}
