//! Rahman's liquid again, with forces from a neighbor list split across
//! threads (forces.rs), against the same run with forces from all pairs.
//!
//! M4's argon program (src/bin/argon) melts 864 atoms in Rahman's box with
//! forces from all pairs, lets them settle for 2000 steps and measures 10
//! blocks. This test melts them the same way, then switches to the list
//! for the settling and the blocks. The two runs start from the same atoms
//! with the same energy and part only as rounding differences grow through
//! chaos (tests/scale.rs), so they are two runs of the same liquid, and
//! their averages must agree within their error bars. It takes about a
//! minute with --release, so a plain `cargo test` skips it. Run it with
//! `cargo test --release -- --ignored`.

use omega::forces::Method;
use omega::pairs;
use omega::rahman::{self, CUTOFF, DT};
use omega::run::{Block, Run};
use omega::system::System;
use omega::vec3::Vec3;

/// M4's argon run with forces from all pairs, in reduced units (the argon
/// program prints them in K and cm^2 / s): the mean over its 10 blocks and
/// its standard error, for the temperature (eps / kB), and D from the MSD
/// and from the VACF (sigma^2 / tau).
const ALL_PAIRS: [(f64, f64); 3] = [(0.78994, 0.00027), (0.04785, 0.00046), (0.04762, 0.00040)];

#[test]
#[ignore = "about a minute with --release: cargo test --release -- --ignored"]
fn listed_liquid_matches_the_all_pairs_liquid() {
    // While the liquid settles, the forces must match all pairs every 10
    // steps. Measured: within 7.2e-14 eps / sigma. Tolerance 1e-12, as in
    // tests/forces.rs. Then T and D must agree with all pairs within 4
    // standard errors of their difference, our own choice: about 7% for D
    // and 0.2% for T. Measured: T 0.78976 +- 0.00034, D 0.04839 +-
    // 0.00071 from the MSD and 0.04815 +- 0.00073 from the VACF, at most
    // 0.64 standard errors of the difference from all pairs.
    let melted = rahman::liquid(1964);
    let (r, v, space) = (melted.positions, melted.velocities, melted.space);
    let mut system = System::with_forces(r, v, space, CUTOFF, Method::listed(0.3, 4));
    let mut forces = vec![Vec3::new(0.0, 0.0, 0.0); 864];
    for step in 1..=2000 {
        system.step(DT);
        if step % 10 == 0 {
            pairs::compute(&system.positions, space, CUTOFF, &mut forces);
            let pairs = system.forces.iter().zip(&forces);
            let miss = pairs.map(|(&a, &b)| (a - b).norm()).fold(0.0, f64::max);
            assert!(miss <= 1e-12, "step {step}: forces {miss:e} apart");
        }
    }
    let mut run = Run::new(&system);
    for _ in 0..10 {
        run.add_block(&mut system, DT);
    }
    let picks: [fn(&Block) -> f64; 3] = [|b| b.temperature, |b| b.d_msd, |b| b.d_vacf];
    for (pick, (all, all_error)) in picks.into_iter().zip(ALL_PAIRS) {
        let (listed, error) = run.mean(pick);
        let line = 4.0 * (error * error + all_error * all_error).sqrt();
        let message = format!("{listed} +- {error} against {all} +- {all_error}");
        assert!((listed - all).abs() <= line, "{message}");
    }
}