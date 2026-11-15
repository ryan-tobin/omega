//! Starting positions: the face-centered cubic (fcc) lattice.
//!
//! A cube of side a has an atom at one corner and at the center of each of
//! the three faces that meet there: at (0, 0, 0), (a/2, a/2, 0),
//! (a/2, 0, a/2) and (0, a/2, a/2). Repeating that cube n times along each
//! axis fills a box of side L = n * a with 4 * n^3 atoms. For n = 6 that is
//! 864 atoms, the number Rahman used.
//!
//! In this lattice each atom has 12 nearest neighbors at a / sqrt(2), then
//! 6 at a, then 24 at a * sqrt(3/2). Those counts are the first terms of
//! OEIS A004015, the theta series of the fcc lattice.
//!
//! Label: assumed. Solid argon takes this lattice, which makes it a natural
//! starting arrangement, but the starting positions are put in by hand.

use crate::vec3::Vec3;

/// The four atoms of one cube, as fractions of its side a.
const BASIS: [(f64, f64, f64); 4] = [
    (0.0, 0.0, 0.0),
    (0.5, 0.5, 0.0),
    (0.5, 0.0, 0.5),
    (0.0, 0.5, 0.5),
];

/// Positions of 4 * cells^3 atoms on an fcc lattice that fills a box of
/// side `side`, in sigma. Every position lies in [0, side).
pub fn fcc(cells: usize, side: f64) -> Vec<Vec3> {
    let a = side / cells as f64;
    let mut positions = Vec::with_capacity(4 * cells * cells * cells);
    for ix in 0..cells {
        for iy in 0..cells {
            for iz in 0..cells {
                let corner = Vec3::new(ix as f64, iy as f64, iz as f64);
                for (bx, by, bz) in BASIS {
                    positions.push((corner + Vec3::new(bx, by, bz)) * a);
                }
            }
        }
    }
    positions
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::periodic::PeriodicBox;

    #[test]
    fn six_cells_give_864_atoms_inside_the_box() {
        let side = 10.229;
        let positions = fcc(6, side);
        assert_eq!(positions.len(), 864);
        for r in positions {
            for c in [r.x, r.y, r.z] {
                assert!((0.0..side).contains(&c));
            }
        }
    }

    #[test]
    fn neighbor_shells_hold_12_then_6_then_24_atoms() {
        // Cube side a = 1 in a box of side 3. Every position is a whole
        // number or a half, so every squared distance is exact: 0.5 for
        // the first shell, 1 for the second and 1.5 for the third. The
        // third shell, at 1.22, is inside L/2 = 1.5, so the nearest image
        // gives each neighbor once.
        let space = PeriodicBox { side: 3.0 };
        let positions = fcc(3, space.side);
        for (i, &ri) in positions.iter().enumerate() {
            let mut shells = [0; 3];
            for (j, &rj) in positions.iter().enumerate() {
                if i == j {
                    continue;
                }
                let r2 = space.nearest_image(ri - rj).norm_sq();
                assert!(r2 >= 0.5, "atoms {i} and {j} are too close");
                for (k, shell) in [0.5, 1.0, 1.5].into_iter().enumerate() {
                    if r2 == shell {
                        shells[k] += 1;
                    }
                }
            }
            assert_eq!(shells, [12, 6, 24], "atom {i}");
        }
    }
}
