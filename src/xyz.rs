//! Atom positions written as extended XYZ text, a format other tools read
//! (ASE, and Open Babel as plain XYZ), so they can check what we compute.
//!
//! A frame is the number of atoms on one line, then one line of key=value
//! pairs, then one line per atom: its element and its x, y and z. The
//! pairs give the box (Lattice: its three edges, one after another), the
//! columns of the atom lines (Properties) and the time (Time). This is
//! the layout of the example in ASE's extended XYZ reader, ase/io/extxyz.py.
//! A reader of plain XYZ keeps that line as a comment. A file can hold
//! many frames, one after another.
//!
//! XYZ coordinates are read as angstroms (Open Babel's XYZ format notes,
//! and ASE), so positions and the box are written in angstroms. Positions
//! are wrapped back into the box, where a viewer expects them. The time
//! is written in picoseconds.

use crate::periodic::PeriodicBox;
use crate::units::Units;
use crate::vec3::Vec3;
use std::io::{self, Write};

/// Writes one frame: atoms of element `symbol` at `positions`, in sigma,
/// in the box `space`, at time `time`, in tau. `units` turns lengths into
/// angstroms and the time into picoseconds, written with 6 decimals.
///
/// Each coordinate is first wrapped into the box: x - L floor(x / L),
/// computed by `f64::rem_euclid`. Its documentation says that for a
/// negative x very close to 0 the result can be L itself, not 0. A viewer
/// draws both at the same face of the box.
pub fn write_frame(
    out: &mut impl Write,
    symbol: &str,
    positions: &[Vec3],
    space: PeriodicBox,
    time: f64,
    units: Units,
) -> io::Result<()> {
    let angstrom = units.sigma_m * 1e10;
    let side = space.side * angstrom;
    let ps = time * units.time_s() * 1e12;
    writeln!(out, "{}", positions.len())?;
    writeln!(
        out,
        "Lattice=\"{side:.6} 0.0 0.0 0.0 {side:.6} 0.0 0.0 0.0 {side:.6}\" \
         Properties=species:S:1:pos:R:3 Time={ps:.6}"
    )?;
    let wrap = |c: f64| c.rem_euclid(space.side) * angstrom;
    for r in positions {
        let (x, y, z) = (wrap(r.x), wrap(r.y), wrap(r.z));
        writeln!(out, "{symbol} {x:.6} {y:.6} {z:.6}")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::ARGON;

    #[test]
    fn frame_of_two_argon_atoms() {
        // A box of side 8 sigma is 27.2 angstrom. The second atom starts
        // outside the box and is wrapped to (7.5, 0.5, 4) sigma. One sigma
        // is 3.4 angstrom and one tau is 2.151442 ps, so t = 2 tau is
        // 4.302884 ps.
        let positions = [Vec3::new(1.0, 2.0, 3.0), Vec3::new(-0.5, 8.5, 4.0)];
        let space = PeriodicBox { side: 8.0 };
        let mut out = Vec::new();
        write_frame(&mut out, "Ar", &positions, space, 2.0, ARGON).unwrap();
        let text = String::from_utf8(out).unwrap();
        let want = "2\n\
            Lattice=\"27.200000 0.0 0.0 0.0 27.200000 0.0 0.0 0.0 27.200000\" \
            Properties=species:S:1:pos:R:3 Time=4.302884\n\
            Ar 3.400000 6.800000 10.200000\n\
            Ar 25.500000 1.700000 13.600000\n";
        assert_eq!(text, want);
    }
}
