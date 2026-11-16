//! Validation against a published result: liquid argon as Rahman measured
//! it, Phys. Rev. 136, A405 (1964) (rahman.rs has his numbers and where
//! they come from).
//!
//! The run is the argon program's, cut to its first 5 blocks (100 tau):
//! melt and hold 864 atoms at his temperature, settle for 10 tau, then
//! measure. It is made once and shared by the tests below. It takes about
//! a minute and a half with --release and far longer without, so these
//! tests are ignored by a plain `cargo test`. Run them with
//! `cargo test --release -- --ignored`.
//!
//! Every pass line is our own choice and says why.

use omega::peaks;
use omega::rahman::{self, DT, SIDE};
use omega::run::{self, Run};
use omega::structure;
use omega::units::ARGON;
use std::sync::OnceLock;

const BLOCKS: usize = 5;

/// The shared run, made by the first test that needs it.
fn rahman_run() -> &'static Run {
    static RUN: OnceLock<Run> = OnceLock::new();
    RUN.get_or_init(|| {
        let mut system = rahman::liquid(1964);
        for _ in 0..2000 {
            system.step(DT);
        }
        let mut run = Run::new(&system);
        for _ in 0..BLOCKS {
            run.add_block(&mut system, DT);
        }
        run
    })
}

/// sigma^2 / tau in cm^2 / s.
fn cm2_per_s() -> f64 {
    let cm = ARGON.sigma_m * 100.0;
    cm * cm / ARGON.time_s()
}

#[test]
#[ignore = "about 90 s with --release: cargo test --release -- --ignored"]
fn diffusion_is_within_10_percent_of_rahman() {
    // PHYSICS_REFERENCE's pass line. Measured: 2.54e-5 cm^2/s, 4.5% above
    // his 2.43e-5, with a standard error of 0.8% over the 5 blocks.
    let (d, _) = rahman_run().mean(|b| b.d_msd);
    let off = d * cm2_per_s() / rahman::DIFFUSION_CM2_PER_S - 1.0;
    assert!(
        off.abs() <= 0.10,
        "D {:.3e} cm^2/s, {off:+.3}",
        d * cm2_per_s()
    );
}

#[test]
#[ignore = "about 90 s with --release: cargo test --release -- --ignored"]
fn diffusion_from_the_vacf_agrees_with_the_msd() {
    // Two ways to the same D: the slope of the MSD and the integral of the
    // VACF. They use the same motion, so they should agree to well within
    // the error of either. Measured: 0.2% apart. Pass line 3%.
    let run = rahman_run();
    let (from_msd, _) = run.mean(|b| b.d_msd);
    let (from_vacf, _) = run.mean(|b| b.d_vacf);
    let off = from_vacf / from_msd - 1.0;
    assert!(off.abs() <= 0.03, "MSD {from_msd}, VACF {from_vacf}");
}

#[test]
#[ignore = "about 90 s with --release: cargo test --release -- --ignored"]
fn first_peak_of_g_is_where_rahman_found_it() {
    // His first peak is quoted to 0.01 angstrom at best. Ours moves by
    // about 0.03 angstrom with the bins and the stretch of run. Measured:
    // 3.697, 0.027 angstrom over his. Pass line 0.1 angstrom, 3% of sigma.
    let (r, _) = peaks::highest(&rahman_run().rdf.g(), 0.9, 1.4);
    let r = r * ARGON.sigma_m * 1e10;
    assert!(
        (r - rahman::FIRST_PEAK_ANGSTROM).abs() <= 0.1,
        "first peak {r} angstrom"
    );
}

#[test]
#[ignore = "about 90 s with --release: cargo test --release -- --ignored"]
fn peaks_of_s_are_where_rahman_found_them() {
    // His four peaks are given to 0.1 in k sigma. Ours carry ripples from
    // cutting g(r) off at 5.1 sigma. Measured: within 1.2%. Pass line 2%.
    let g = rahman_run().rdf.g();
    let density = 864.0 / (SIDE * SIDE * SIDE);
    let s: Vec<(f64, f64)> = (1..=600)
        .map(|i| 0.05 * f64::from(i))
        .map(|k| (k, structure::structure_factor(&g, density, k)))
        .collect();
    for (&(lo, hi), his) in rahman::S_PEAK_WINDOWS.iter().zip(rahman::S_PEAKS) {
        let (k, _) = peaks::highest(&s, lo, hi);
        assert!((k / his - 1.0).abs() <= 0.02, "peak at {k}, his {his}");
    }
}

#[test]
#[ignore = "about 90 s with --release: cargo test --release -- --ignored"]
fn vacf_turns_negative_when_rahman_found() {
    // His time is given to 0.01 ps, and his time step was 0.01 ps. Ours is
    // read from samples 0.054 ps apart. Measured: 0.318 ps, 3.7% under his.
    // Pass line 10%.
    let run = rahman_run();
    let vacf = run.curve(|b| &b.vacf);
    let h = run::SAMPLE_EVERY as f64 * DT;
    let points: Vec<(f64, f64)> = vacf
        .iter()
        .enumerate()
        .map(|(k, &v)| (k as f64 * h, v))
        .collect();
    let zero = peaks::first_below(&points, 0.0).expect("the VACF never turns negative");
    let ps = zero * ARGON.time_s() * 1e12;
    assert!(
        (ps / rahman::VACF_ZERO_PS - 1.0).abs() <= 0.10,
        "first zero {ps} ps"
    );
}
