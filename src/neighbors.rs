//! The Verlet neighbor list: for each atom, the atoms within r_c + s of
//! it, where r_c is the cutoff and s a margin called the skin.
//!
//! The forces need only the pairs within r_c. A list that also holds the
//! pairs out to r_c + s serves for many steps. If no atom has moved more
//! than s / 2 since the list was built, no pair has come closer by more
//! than s (triangle inequality), so a pair that was not on the list is
//! still at least r_c apart. The list is built again the first time some
//! atom has moved farther than that. L. Verlet, Phys. Rev. 159, 98 (1967),
//! is credited with the list (M. P. Allen, Introduction to Molecular
//! Dynamics Simulation (2004)); his paper was not read. LAMMPS rebuilds by
//! the same rule by default, and its default skin in reduced units is 0.3
//! sigma (its neigh_modify and neighbor commands).
//!
//! The list is built from cells at least r_c + s wide (cells.rs), so a
//! build costs a fixed amount per atom. Each pair is listed twice, once
//! under each of its atoms, so that each thread of forces.rs can work out
//! the forces on its own atoms alone. Positions are never wrapped, so an
//! atom has moved |r - r_built| since the build.
//!
//! Label: derived. While no atom has moved more than s / 2, the list holds
//! every pair within r_c (tests/forces.rs checks the forces from it).

use crate::cells::Cells;
use crate::parallel;
use crate::periodic::PeriodicBox;
use crate::vec3::Vec3;

/// For each atom, the atoms near it.
pub struct NeighborList {
    /// The skin s, in sigma.
    pub skin: f64,
    /// How many times the list has been built. Unit: none.
    pub builds: usize,
    /// The positions when the list was last built, in sigma.
    pub built_at: Vec<Vec3>,
    /// The neighbors of each atom, as atom indices.
    near: Vec<Vec<usize>>,
}

impl NeighborList {
    /// An empty list with a skin of `skin` sigma, not yet built.
    pub fn new(skin: f64) -> Self {
        assert!(skin >= 0.0, "the skin cannot be negative");
        Self {
            skin,
            builds: 0,
            built_at: Vec::new(),
            near: Vec::new(),
        }
    }

    /// True if the list was not built for this many atoms, or if some atom
    /// has moved more than half the skin since it was. `positions` in sigma.
    pub fn is_stale(&self, positions: &[Vec3]) -> bool {
        let limit = 0.25 * self.skin * self.skin;
        positions.len() != self.built_at.len()
            || positions
                .iter()
                .zip(&self.built_at)
                .any(|(&r, &r0)| (r - r0).norm_sq() > limit)
    }

    /// Lists, for each atom of `positions` (in sigma), every other atom
    /// whose nearest copy is closer than `cutoff` + skin, in sigma,
    /// searching only the cells around it. The atoms are split across
    /// `threads` threads. Each atom's neighbors are in the order the cells
    /// give them, the same for any number of threads.
    pub fn build(&mut self, positions: &[Vec3], space: PeriodicBox, cutoff: f64, threads: usize) {
        let reach = cutoff + self.skin;
        let reach_sq = reach * reach;
        let cells = Cells::new(positions, space, reach);
        let m = cells.per_side;
        let around: Vec<Vec<usize>> = (0..m * m * m).map(|c| cells.around(c)).collect();
        self.near.resize_with(positions.len(), Vec::new);
        parallel::for_each(&mut self.near, threads, |i, near| {
            near.clear();
            for &c in &around[cells.of_atom[i]] {
                for &j in cells.atoms_in(c) {
                    let r = space.nearest_image(positions[i] - positions[j]);
                    if j != i && r.norm_sq() < reach_sq {
                        near.push(j);
                    }
                }
            }
        });
        self.built_at.clear();
        self.built_at.extend_from_slice(positions);
        self.builds += 1;
    }

    /// The neighbors of atom `i`, as atom indices.
    pub fn of(&self, i: usize) -> &[usize] {
        &self.near[i]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Rng;

    /// `n` atoms at random from -L to 2L along each axis of a box of side
    /// `side`, so most are outside it.
    fn scattered(n: usize, side: f64, seed: u64) -> Vec<Vec3> {
        let mut rng = Rng::new(seed);
        let mut x = || (3.0 * rng.uniform() - 1.0) * side;
        (0..n).map(|_| Vec3::new(x(), x(), x())).collect()
    }

    #[test]
    fn list_holds_the_same_pairs_as_checking_all_of_them() {
        // Rahman's density in boxes of 5, 3 and 2 fcc cells: 500, 108 and
        // 32 atoms. The first gives 3 cells a side, the second 2, and the
        // third, with r_c + s = 2.0 over half its side, 1. Exact.
        for (n, side, cutoff, seed) in [
            (500, 8.524, 2.25, 1),
            (108, 5.1145, 2.25, 2),
            (32, 3.41, 1.7, 3),
        ] {
            let space = PeriodicBox { side };
            let positions = scattered(n, side, seed);
            let mut list = NeighborList::new(0.3);
            list.build(&positions, space, cutoff, 3);
            let reach_sq = (cutoff + 0.3) * (cutoff + 0.3);
            for (i, &ri) in positions.iter().enumerate() {
                let near = |&j: &usize| space.nearest_image(ri - positions[j]).norm_sq() < reach_sq;
                let want: Vec<usize> = (0..n).filter(|&j| j != i).filter(near).collect();
                let mut got = list.of(i).to_vec();
                got.sort_unstable();
                assert_eq!(got, want, "atom {i} of {n}");
            }
        }
    }

    #[test]
    fn list_goes_stale_when_an_atom_moves_over_half_the_skin() {
        let space = PeriodicBox { side: 8.524 };
        let mut positions = scattered(500, space.side, 4);
        let mut list = NeighborList::new(0.3);
        assert!(list.is_stale(&positions), "never built");
        list.build(&positions, space, 2.25, 2);
        assert!(!list.is_stale(&positions));
        positions[7] = positions[7] + Vec3::new(0.0, 0.149, 0.0);
        assert!(!list.is_stale(&positions), "moved 0.149");
        positions[7] = positions[7] + Vec3::new(0.0, 0.002, 0.0);
        assert!(list.is_stale(&positions), "moved 0.151");
    }
}
