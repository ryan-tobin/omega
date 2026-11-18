//! The ways to compute the forces: from all pairs (pairs.rs), from a
//! neighbor list (neighbors.rs) with the atoms split across threads
//! (parallel.rs), or, with the `gpu` feature, on a GPU (gpu/forces.rs).
//!
//! With a list, the force on each atom is summed over its own neighbors.
//! So every pair is worked out twice, once for each of its atoms, where
//! pairs.rs works it out once and gives the other atom the opposite force.
//! That doubles the work per pair, but lets each thread write the forces
//! of its own atoms only. The force of each pair is the same number either
//! way: the nearest copy of r_j - r_i is exactly minus that of r_i - r_j,
//! and lj::force(-r) is exactly -lj::force(r). Only the order of the sums
//! differs, so the two ways agree to rounding (tested). Each atom's sum
//! runs in the order of its list, and the energy and virial are added up
//! atom by atom once every thread has finished, so the result is the same
//! to the last bit for any number of threads (tested).
//!
//! Label: approximated, as in pairs.rs: cut and shifted at r_c.

use crate::lj;
use crate::neighbors::NeighborList;
use crate::pairs;
use crate::parallel;
use crate::periodic::PeriodicBox;
use crate::vec3::Vec3;

const ZERO: Vec3 = Vec3::new(0.0, 0.0, 0.0);

/// How a system's forces are computed.
pub enum Method {
    /// Every pair, on one thread (pairs.rs).
    AllPairs,
    /// A neighbor list, built again whenever it goes stale, with the work
    /// split across `threads` threads.
    Listed { list: NeighborList, threads: usize },
    /// On the GPU, in 32-bit floats (gpu/forces.rs)
    #[cfg(feature = "gpu")]
    Gpu(Box<crate::gpu::forces::GpuForces>),
}

impl Method {
    /// A neighbor list with a skin of `skin` sigma, not yet built, and the
    /// work split across `threads` threads, at least 1.
    pub fn listed(skin: f64, threads: usize) -> Self {
        assert!(threads >= 1, "need at least one thread");
        Self::Listed {
            list: NeighborList::new(skin),
            threads,
        }
    }

    /// Fills `forces` with the total force on each atom, in eps / sigma,
    /// and returns the total potential energy and the virial W, both in
    /// eps, as pairs::compute does. `positions` are in sigma, and `cutoff`
    /// is r_c in sigma, at most half the box side and the same at every
    /// call. With a list, it is first built again if it has gone stale.
    pub fn compute(
        &mut self,
        positions: &[Vec3],
        space: PeriodicBox,
        cutoff: f64,
        forces: &mut [Vec3],
    ) -> (f64, f64) {
        match self {
            Self::AllPairs => pairs::compute(positions, space, cutoff, forces),
            Self::Listed { list, threads } => {
                if list.is_stale(positions) {
                    list.build(positions, space, cutoff, *threads);
                }
                listed(positions, space, cutoff, list, *threads, forces)
            }
            #[cfg(feature = "gpu")]
            Self::Gpu(gpu) => gpu.compute(positions, space, cutoff, forces),
        }
    }

    /// How many times the neighbor list has been built, 0 for all pairs.
    pub fn builds(&self) -> usize {
        match self {
            Self::Listed { list, .. } => list.builds,
            _ => 0,
        }
    }
}

/// The forces, potential energy and virial of `compute`, from the pairs
/// on `list`, which must be up to date.
fn listed(
    positions: &[Vec3],
    space: PeriodicBox,
    cutoff: f64,
    list: &NeighborList,
    threads: usize,
    forces: &mut [Vec3],
) -> (f64, f64) {
    assert_eq!(forces.len(), positions.len());
    assert!(cutoff <= 0.5 * space.side, "cutoff is over half the box");
    let cutoff_sq = cutoff * cutoff;
    let shift = lj::energy(Vec3::new(cutoff, 0.0, 0.0));
    // Each atom's force, and its shares of the energy and the virial.
    let mut each = vec![(ZERO, 0.0, 0.0); positions.len()];
    parallel::for_each(&mut each, threads, |i, out| {
        let (mut f, mut u, mut w) = (ZERO, 0.0, 0.0);
        for &j in list.of(i) {
            let r_ij = space.nearest_image(positions[i] - positions[j]);
            if r_ij.norm_sq() < cutoff_sq {
                let f_ij = lj::force(r_ij);
                f = f + f_ij;
                u += lj::energy(r_ij) - shift;
                w += r_ij.dot(f_ij);
            }
        }
        *out = (f, u, w);
    });
    let (mut potential, mut virial) = (0.0, 0.0);
    for (force, (f, u, w)) in forces.iter_mut().zip(each) {
        *force = f;
        potential += u;
        virial += w;
    }
    // Each pair was counted once under each of its two atoms.
    (0.5 * potential, 0.5 * virial)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice;

    #[test]
    fn list_is_built_again_only_once_it_goes_stale() {
        // 500 lattice atoms at Rahman's density. Moving one atom 0.16
        // sigma, over half the skin of 0.3, must bring a second build.
        let space = PeriodicBox {
            side: 10.229 / 6.0 * 5.0,
        };
        let mut positions = lattice::fcc(5, space.side);
        let mut forces = vec![ZERO; positions.len()];
        let mut method = Method::listed(0.3, 2);
        assert_eq!(method.builds(), 0);
        method.compute(&positions, space, 2.25, &mut forces);
        method.compute(&positions, space, 2.25, &mut forces);
        assert_eq!(method.builds(), 1);
        positions[0] = positions[0] + Vec3::new(0.16, 0.0, 0.0);
        method.compute(&positions, space, 2.25, &mut forces);
        assert_eq!(method.builds(), 2);
        assert_eq!(Method::AllPairs.builds(), 0);
    }
}
