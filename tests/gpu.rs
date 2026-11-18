//! Forces from the GPU (gpu/forces.rs) against the CPU at one moment: a
//! pair with a known answer, and all pairs (pairs.rs) on shaken lattices.
//!
//! These tests need a GPU and the `gpu` feature: cargo test --features gpu.
//! Every tolerance is our own choice and says why. One rounding step of a
//! 32-bit float is 6e-8 of the number, so the lines here are 1e-5 where
//! the CPU tests (tests/forces.rs) have 1e-12.
#![cfg(feature = "gpu")]

use omega::forces::Method;
use omega::gpu::Gpu;
use omega::gpu::forces::GpuForces;
use omega::lattice;
use omega::pairs;
use omega::periodic::PeriodicBox;
use omega::rahman::{CUTOFF, SIDE};
use omega::rng::Rng;
use omega::vec3::Vec3;
use std::sync::OnceLock;

const ZERO: Vec3 = Vec3::new(0.0, 0.0, 0.0);
// A pair this close to the cutoff, in sigma, can count for the GPU and not
// for the CPU, or the other way round: the GPU rounds each separation to
// about 1e-7 sigma. Such a pair moves two forces by 0.081 eps / sigma, the
// force at the cutoff (PHYSICS_REFERENCE). So each atom's force is compared
// with all pairs at a cutoff this much shorter and this much longer, and
// the nearer of the two counts.
const EDGE: f64 = 1e-6;

/// Forces computed on this machine's GPU, shared by all tests.
fn on_gpu() -> Method {
    static GPU: OnceLock<Gpu> = OnceLock::new();
    let gpu = GPU.get_or_init(|| Gpu::new().expect("these tests need a GPU"));
    Method::Gpu(Box::new(GpuForces::new(gpu)))
}

/// The farthest any entry of `a` is from the same entry of `b`.
fn apart(a: &[Vec3], b: &[Vec3]) -> f64 {
    let pairs = a.iter().zip(b);
    pairs.map(|(&p, &q)| (p - q).norm()).fold(0.0, f64::max)
}

/// Forces from `gpu` against all pairs at `positions`: the largest miss of
/// a force as a fraction of the largest force, and the relative misses of
/// the potential energy and of the virial.
fn misses(gpu: &mut Method, positions: &[Vec3], space: PeriodicBox) -> [f64; 3] {
    let zeros = || vec![ZERO; positions.len()];
    let (mut got, mut exact, mut inner, mut outer) = (zeros(), zeros(), zeros(), zeros());
    let (u, w) = gpu.compute(positions, space, CUTOFF, &mut got);
    let (u0, _) = pairs::compute(positions, space, CUTOFF, &mut exact);
    let (_, w_in) = pairs::compute(positions, space, CUTOFF - EDGE, &mut inner);
    let (_, w_out) = pairs::compute(positions, space, CUTOFF + EDGE, &mut outer);
    let miss = |i: usize| (got[i] - inner[i]).norm().min((got[i] - outer[i]).norm());
    let force = (0..got.len()).map(miss).fold(0.0, f64::max) / apart(&exact, &zeros());
    let virial = (w / w_in - 1.0).abs().min((w / w_out - 1.0).abs());
    [force, (u / u0 - 1.0).abs(), virial]
}

#[test]
fn two_atoms_one_sigma_apart_across_the_boundary_feel_24() {
    // The case of pairs.rs, in a box of side 8 with cutoff 2: x = 0.25 and
    // x = 7.25 are 1 sigma apart through the face of the box. There
    // F = 24, U = 0 and r . F = 24 (lj.rs), and the shifted energy is
    // -U(2) = 63/1024. Every step is exact in an f32, and the software GPU
    // here gives these values exactly. The tolerance, 1e-5, leaves a real
    // GPU room to round its division.
    let positions = [Vec3::new(0.25, 4.0, 4.0), Vec3::new(7.25, 4.0, 4.0)];
    let mut forces = [ZERO; 2];
    let (u, w) = on_gpu().compute(&positions, PeriodicBox { side: 8.0 }, 2.0, &mut forces);
    let want = [Vec3::new(24.0, 0.0, 0.0), Vec3::new(-24.0, 0.0, 0.0)];
    assert!(apart(&forces, &want) <= 1e-5, "{forces:?}");
    assert!((u - 63.0 / 1024.0).abs() <= 1e-5, "{u}");
    assert!((w - 24.0).abs() <= 1e-5, "{w}");
}

#[test]
fn the_two_atoms_of_a_pair_get_exactly_opposite_forces() {
    // Two atoms 1.2 sigma apart through a face of the box, at positions
    // that are no whole number of ticks. Each force is a single term, so
    // no sum is rounded: the two must be opposite to the last bit
    // (forces.wgsl). By lj.rs the force is about 2.3 eps / sigma.
    let positions = [Vec3::new(0.3, 4.1, 3.9), Vec3::new(7.4, 4.5, 3.2)];
    let mut forces = [ZERO; 2];
    on_gpu().compute(&positions, PeriodicBox { side: 8.0 }, 2.0, &mut forces);
    assert!((forces[0].norm() - 2.3).abs() <= 0.05, "{forces:?}");
    assert_eq!(forces[0], -forces[1]);
}

#[test]
fn gpu_forces_agree_with_all_pairs() {
    // Lattices of 108, 500 and 864 atoms (3, 5 and 6 cells a side), each
    // atom moved at random by up to 0.1 sigma along each axis. Then the
    // same atoms with two in three carried whole box lengths away, which
    // must change nothing. One GpuForces serves every size, so its buffers
    // are reused and made again. Measured on a software GPU over 72 such
    // arrangements: forces within 1.6e-6 of the largest force, energy
    // within 1.5e-7 and virial within 6.2e-7. Lines 1e-5: at least 6 times
    // the worst measured, and far below the 1e-3 that one missed pair
    // would leave.
    let mut gpu = on_gpu();
    for cells in [3, 5, 6] {
        let space = PeriodicBox {
            side: SIDE / 6.0 * cells as f64,
        };
        let mut rng = Rng::new(9);
        let mut positions = lattice::fcc(cells, space.side);
        for r in &mut positions {
            let (x, y, z) = (rng.uniform(), rng.uniform(), rng.uniform());
            *r = *r + (Vec3::new(x, y, z) - Vec3::new(0.5, 0.5, 0.5)) * 0.2;
        }
        let found = misses(&mut gpu, &positions, space);
        assert!(found.iter().all(|&m| m <= 1e-5), "{cells}: {found:?}");
        let away = Vec3::new(space.side, -2.0 * space.side, 3.0 * space.side);
        for (i, r) in positions.iter_mut().enumerate() {
            *r = *r + away * (i % 3) as f64;
        }
        let found = misses(&mut gpu, &positions, space);
        assert!(found.iter().all(|&m| m <= 1e-5), "{cells} moved: {found:?}");
    }
}
