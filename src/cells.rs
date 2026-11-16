//! Cells: the box cut into m x m x m smaller cubes, with each atom filed
//! under the cube it sits in.
//!
//! If every cell is at least r wide, the atoms within r of an atom can
//! only be in its own cell or in the 26 cells around it. A cell holds
//! about rho w^3 atoms however large the box is, so filing the atoms and
//! searching the cells around each one costs a fixed amount per atom. The
//! work grows in step with N, where checking all pairs (pairs.rs) grows as
//! N^2. This is the cell method. M. P. Allen (2004) credits it to R. W.
//! Hockney and J. W. Eastwood, Computer Simulation Using Particles (1981),
//! who chain each cell's atoms in a linked list, and others to B. Quentrec
//! and C. Brot, J. Comput. Phys. 13, 430 (1973). Neither was read here.
//! We sort the atoms by cell instead (a counting sort).
//!
//! An atom is filed by its position wrapped into the box, so atoms that
//! have left the box (system.rs never wraps them) land in the right cell.
//! With fewer than 3 cells along a side, the 27 cells around a cell
//! include some cells twice, through the repeating box. Each is kept once.
//!
//! Nothing here changes the physics: the search finds exactly the pairs
//! that checking every pair finds (tested in neighbors.rs, with 3, 2 and
//! 1 cells a side).

use crate::periodic::PeriodicBox;
use crate::vec3::Vec3;

/// The atoms of a box, filed by cell.
pub struct Cells {
    /// Cells along each side of the box, m.
    pub per_side: usize,
    /// The cell of each atom, by atom index.
    pub of_atom: Vec<usize>,
    /// The atoms of cell c are atoms[start[c]..start[c + 1]].
    start: Vec<usize>,
    atoms: Vec<usize>,
}

impl Cells {
    /// Files every atom of `positions`, in sigma, in cells at least
    /// `width` sigma wide: as many cells along each side of `space` as fit,
    /// and at least 1.
    pub fn new(positions: &[Vec3], space: PeriodicBox, width: f64) -> Self {
        let l = space.side;
        let m = ((l / width) as usize).max(1);
        let w = l / m as f64;
        // rem_euclid can round a tiny negative x up to exactly L, which
        // would be cell m, so the index is held to m - 1.
        let along = |x: f64| ((x.rem_euclid(l) / w) as usize).min(m - 1);
        let cell = |r: &Vec3| (along(r.x) * m + along(r.y)) * m + along(r.z);
        let of_atom: Vec<usize> = positions.iter().map(cell).collect();
        // Count the atoms of each cell, turn the counts into where each
        // cell starts, then place the atoms in index order.
        let mut start = vec![0; m * m * m + 1];
        for &c in &of_atom {
            start[c] += 1;
        }
        let mut sum = 0;
        for s in &mut start {
            let count = *s;
            *s = sum;
            sum += count;
        }
        let mut next = start.clone();
        let mut atoms = vec![0; positions.len()];
        for (i, &c) in of_atom.iter().enumerate() {
            atoms[next[c]] = i;
            next[c] += 1;
        }
        Self {
            per_side: m,
            of_atom,
            start,
            atoms,
        }
    }

    /// The atoms in cell `c`, lowest index first.
    pub fn atoms_in(&self, c: usize) -> &[usize] {
        &self.atoms[self.start[c]..self.start[c + 1]]
    }

    /// Cell `c` and the 26 cells around it, through the repeating box,
    /// each once, in increasing order.
    pub fn around(&self, c: usize) -> Vec<usize> {
        let m = self.per_side;
        let (x, y, z) = (c / (m * m), c / m % m, c % m);
        let mut cells = Vec::with_capacity(27);
        for dx in [m - 1, 0, 1] {
            for dy in [m - 1, 0, 1] {
                for dz in [m - 1, 0, 1] {
                    cells.push(((x + dx) % m * m + (y + dy) % m) * m + (z + dz) % m);
                }
            }
        }
        cells.sort_unstable();
        cells.dedup();
        cells
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Rng;

    const RAHMAN: PeriodicBox = PeriodicBox { side: 10.229 };

    #[test]
    fn every_atom_is_filed_once_in_the_cell_it_sits_in() {
        // Rahman's box with cells at least 2.55 wide: 10.229 / 2.55 = 4.01,
        // so 4 a side, 2.557 wide. 1000 atoms at random from -L to 2L
        // along each axis, so most are outside the box. Wrapped, each must
        // lie inside its cell, and each cell must list its atoms once,
        // lowest index first. 1e-12 allows for rounding on the faces
        // between cells, our own choice.
        let mut rng = Rng::new(5);
        let l = RAHMAN.side;
        let mut random = || (3.0 * rng.uniform() - 1.0) * l;
        let positions: Vec<Vec3> = (0..1000)
            .map(|_| Vec3::new(random(), random(), random()))
            .collect();
        let cells = Cells::new(&positions, RAHMAN, 2.55);
        assert_eq!(cells.per_side, 4);
        let mut listed = Vec::new();
        for c in 0..64 {
            let atoms = cells.atoms_in(c);
            assert!(atoms.windows(2).all(|p| p[0] < p[1]), "cell {c}");
            listed.extend_from_slice(atoms);
            for &i in atoms {
                assert_eq!(cells.of_atom[i], c);
                let r = positions[i];
                for (x, k) in [(r.x, c / 16), (r.y, c / 4 % 4), (r.z, c % 4)] {
                    let (x, low) = (x.rem_euclid(l), k as f64 * l / 4.0);
                    assert!(low - 1e-12 <= x && x <= low + l / 4.0 + 1e-12);
                }
            }
        }
        listed.sort_unstable();
        assert_eq!(listed, (0..1000).collect::<Vec<_>>());
    }

    #[test]
    fn an_atom_just_below_zero_lands_in_the_last_cell() {
        // -1e-17 wrapped is L - 1e-17, which rounds to exactly L. It must
        // be filed in the last cell along each axis, 63 of 64, not past it.
        let cells = Cells::new(&[Vec3::new(-1e-17, -1e-17, -1e-17)], RAHMAN, 2.55);
        assert_eq!(cells.of_atom, [63]);
    }
}