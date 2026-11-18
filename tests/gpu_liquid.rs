//! Rahman's liquid with forces from the GPU (gpu/forces.rs), against the
//! same run with forces from all pairs on the CPU.
//!
//! As in tests/listed_liquid.rs: the 864 atoms are melted with forces from
//! all pairs (rahman.rs), then the GPU takes over for the 2000 settling
//! steps and the 10 measured blocks of M4's argon program. The two runs
//! start from the same atoms and part within a few hundred steps
//! (tests/gpu_runs.rs), so they are two runs of the same liquid, and their
//! averages must agree within their error bars. It needs a GPU and takes
//! a minute or two with --release:
//! cargo test --release --features gpu -- --ignored

#![cfg(feature = "gpu")]

use omega::forces::Method;
use omega::gpu::Gpu;
use omega::gpu::forces::GpuForces;
use omega::rahman::{self, CUTOFF, DT};
use omega::run::{Block, Run};
use omega::system::System;

/// M4's argon run with forces from all pairs, in reduced units, as in
/// tests/listed_liquid.rs: the mean over its 10 blocks and its standard
/// error, for the temperature (eps / kB), and D from the MSD and from the
/// VACF (sigma^2 / tau).
const ALL_PAIRS: [(f64, f64); 3] = [(0.78994, 0.00027), (0.04785, 0.00046), (0.04762, 0.00040)];

#[test]
#[ignore = "needs a GPU and a minute or two: cargo test --release --features gpu -- --ignored"]
fn gpu_liquid_matches_the_all_pairs_liquid() {
    // T and D must agree with all pairs within 4 standard errors of their
    // difference, our own choice, as in tests/listed_liquid.rs.
    let gpu = Gpu::new().expect("this test needs a GPU");
    let melted = rahman::liquid(1964);
    let (r, v, space) = (melted.positions, melted.velocities, melted.space);
    let method = Method::Gpu(Box::new(GpuForces::new(&gpu)));
    let mut system = System::with_forces(r, v, space, CUTOFF, method);
    for _ in 0..2000 {
        system.step(DT);
    }
    let mut run = Run::new(&system);
    for _ in 0..10 {
        run.add_block(&mut system, DT);
    }
    let picks: [fn(&Block) -> f64; 3] = [|b| b.temperature, |b| b.d_msd, |b| b.d_vacf];
    for (pick, (all, all_error)) in picks.into_iter().zip(ALL_PAIRS) {
        let (on_gpu, error) = run.mean(pick);
        let line = 4.0 * (error * error + all_error * all_error).sqrt();
        let message = format!("{on_gpu} +- {error} against {all} +- {all_error}");
        assert!((on_gpu - all).abs() <= line, "{message}");
    }
}
