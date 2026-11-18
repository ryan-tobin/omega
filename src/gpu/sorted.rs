//! The atoms sorted by cell, written out as the bytes the GPU reads, and the
//! GPU's answer read back into atom order.
//!
//! The GPU works in 32-bit floats, which hold about 7 digits. A position in
//! a box 100 sigma wide would keep only 5 digits after the point. So each
//! coordinate is sent as a whole number of ticks instead, a tick being
//! 2^-32 of the box side. A 32-bit whole number wraps around at 2^32 just
//! as the box wraps around at its side. The difference of two of them, read
//! as a signed number, is then the separation to the nearest copy
//! (periodic.rs), with no rounding and no box length in the sum. Only that
//! small difference becomes a 32-bit float (forces.wgsl).
//!
//! Label: approximated. Each coordinate is rounded to the nearest tick: by
//! up to 1.2e-9 sigma in Rahman's box, and 1.3e-8 sigma in the box of a
//! million atoms.

use crate::cells::Cells;
use crate::periodic::PeriodicBox;
use crate::vec3::Vec3;

/// Bytes sent for each atom: x, y, z and its cell, each a u32.
pub const ATOM_BYTES: usize = 16;
/// Bytes read back for each atom: five f32 (forces.wgsl).
pub const RESULT_BYTES: usize = 20;
/// Ticks along each side of the box: 2^32.
pub const TICKS: f64 = 4_294_967_296.0;

/// The atoms in cell order. A slot is a place in that order.
pub struct Sorted {
    /// Cells along each side of the box, m.
    pub per_side: usize,
    /// The atom in each slot, as an atom index.
    pub order: Vec<usize>,
    /// ATOM_BYTES for each slot: the atom's position in ticks, and its
    /// cell.
    pub atoms: Vec<u8>,
    /// One u32 for each cell and one more: the slots of cell c run from
    /// entry c up to, not including, entry c + 1.
    pub starts: Vec<u8>,
}

impl Sorted {
    /// Sorts the atoms of `positions`, in sigma, into cells at least `width`
    /// sigma wide (cells.rs), keeping atom order inside each cell.
    pub fn new(positions: &[Vec3], space: PeriodicBox, width: f64) -> Self {
        let cells = Cells::new(positions, space, width);
        let m = cells.per_side;
        // Rounded to a whole number of ticks, then cut to its lowest 32
        // bits. The cut wraps a position outside the box back into it.
        let ticks = |x: f64| ((x / space.side * TICKS).round() as i64 as u32).to_le_bytes();
        let mut order = Vec::with_capacity(positions.len());
        let mut atoms = Vec::with_capacity(positions.len() * ATOM_BYTES);
        let mut starts = Vec::with_capacity((m * m * m + 1) * 4);
        for c in 0..m * m * m {
            starts.extend_from_slice(&(order.len() as u32).to_le_bytes());
            for &i in cells.atoms_in(c) {
                let r = positions[i];
                for word in [ticks(r.x), ticks(r.y), ticks(r.z), (c as u32).to_le_bytes()] {
                    atoms.extend_from_slice(&word);
                }
                order.push(i);
            }
        }
        starts.extend_from_slice(&(order.len() as u32).to_le_bytes());
        Self {
            per_side: m,
            order,
            atoms,
            starts,
        }
    }

    /// Reads the GPU's `results`, RESULT_BYTES for each slot, into `forces`
    /// in atom order, in eps / sigma. Returns the potential energy and the
    /// virial, in eps: half the sums over the atoms, since each pair was
    /// counted under both of its atoms. The sums are made in 64-bit floats.
    pub fn unpack(&self, results: &[u8], forces: &mut [Vec3]) -> (f64, f64) {
        let (mut potential, mut virial) = (0.0, 0.0);
        for (record, &i) in results.chunks_exact(RESULT_BYTES).zip(&self.order) {
            let v: [f64; 5] = std::array::from_fn(|k| {
                let b = &record[4 * k..4 * k + 4];
                f32::from_le_bytes([b[0], b[1], b[2], b[3]]) as f64
            });
            forces[i] = Vec3::new(v[0], v[1], v[2]);
            potential += v[3];
            virial += v[4];
        }
        (0.5 * potential, 0.5 * virial)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const BOX: PeriodicBox = PeriodicBox { side: 8.0 };
    const ATOMS: [Vec3; 3] = [
        Vec3::new(0.500_000_001_4, 0.5, 0.5),
        Vec3::new(-0.5, 9.25, 3.0),
        Vec3::new(2.0, 2.0, 2.0),
    ];

    #[test]
    fn atoms_are_sent_by_cell_in_ticks_that_wrap_like_the_box() {
        // A box of side 8 in cells 2 wide: 4 a side, 64 in all. One sigma
        // is 2^32 / 8 = 2^29 ticks. Atom 0 is in cell 0, and its x is three
        // quarters of a tick past 0.5 sigma, which rounds up. Atom 1 at
        // (-0.5, 9.25, 3) wraps to (7.5, 1.25, 3), which is cell (3, 0, 1),
        // number 49. Atom 2 sits on the corner of cell (1, 1, 1), number
        // 21. So the slots hold atoms 0, 2 and 1.
        let sorted = Sorted::new(&ATOMS, BOX, 2.0);
        assert_eq!(sorted.per_side, 4);
        assert_eq!(sorted.order, [0, 2, 1]);
        // Each slot's x, y, z and cell. A quarter of a sigma is 2^27 ticks.
        let q = 1u32 << 27;
        let rows = [
            [2 * q + 1, 2 * q, 2 * q, 0],
            [8 * q, 8 * q, 8 * q, 21],
            [30 * q, 5 * q, 12 * q, 49],
        ];
        let want: Vec<u8> = rows.concat().iter().flat_map(|w| w.to_le_bytes()).collect();
        assert_eq!(sorted.atoms, want);
        // Cell c starts at the number of atoms in the cells before it, and
        // one more entry closes the last cell.
        let before = |c: u32| [0, 21, 49].iter().filter(|&&cell| cell < c).count() as u32;
        let want: Vec<u8> = (0..65).flat_map(|c| before(c).to_le_bytes()).collect();
        assert_eq!(sorted.starts, want);
    }

    #[test]
    fn results_come_back_in_atom_order_with_the_sums_halved() {
        // Slot 0 holds atom 0, slot 1 atom 2 and slot 2 atom 1. Exact: the
        // energies add to 1 + 2 + 4 = 7 and the virials to 8 + 16 + 32 =
        // 56, and each is halved.
        let sorted = Sorted::new(&ATOMS, BOX, 2.0);
        let rows: [[f32; 5]; 3] = [
            [1.0, 2.0, 3.0, 1.0, 8.0],
            [4.0, 5.0, 6.0, 2.0, 16.0],
            [7.0, 8.0, 9.0, 4.0, 32.0],
        ];
        let results: Vec<u8> = rows.concat().iter().flat_map(|v| v.to_le_bytes()).collect();
        let mut forces = [Vec3::new(0.0, 0.0, 0.0); 3];
        let sums = sorted.unpack(&results, &mut forces);
        assert_eq!(forces[0], Vec3::new(1.0, 2.0, 3.0));
        assert_eq!(forces[1], Vec3::new(7.0, 8.0, 9.0));
        assert_eq!(forces[2], Vec3::new(4.0, 5.0, 6.0));
        assert_eq!(sums, (3.5, 28.0));
    }
}
