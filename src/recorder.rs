//! Saving a run to a folder, for replay.rs to read back.
//!
//! A saved run is three CSV files (csv.rs), in reduced units:
//!
//! - setup.csv, one row: the number of atoms, the box side and the force
//!   cutoff, both in sigma.
//! - thermo.csv, one row per saved moment: the time in tau, then the
//!   temperature (eps / kB), the potential energy per atom (eps), the
//!   pressure (eps / sigma^3) and the total energy (eps). These are the
//!   columns of the measure program.
//! - atoms.csv, one row per atom for each saved moment, in the same
//!   order: its position in sigma, never wrapped, and its velocity in
//!   sigma / tau.
//!
//! csv.rs writes every number exactly, so the files hold the state of
//! the run at each saved moment, not a rounded copy of it.

use crate::csv;
use crate::system::System;
use crate::thermo;
use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::Path;

/// The column names of setup.csv, thermo.csv, and atoms.csv
pub const SETUP: [&str; 3] = ["atoms", "side", "cutoff"];
pub const THERMO: [&str; 5] = ["time", "temperature", "potential", "pressure", "energy"];
pub const ATOMS: [&str; 6] = ["x", "y", "z", "vx", "vy", "vz"];

/// The row of thermo.csv for `system` at `time`, in tau: the time, then
/// the temperature, potential energy per atom, pressure and total energy
/// of the system as it is now, in the units given at the top.
pub fn measure(time: f64, system: &System) -> [f64; 5] {
    let velocities = &system.velocities;
    [
        time,
        thermo::temperature(velocities),
        system.potential / system.positions.len() as f64,
        thermo::pressure(velocities, system.virial, system.space),
        system.total_energy(),
    ]
}

/// Saves a run to a folder, one moment at a time.
pub struct Recorder {
    thermo: BufWriter<File>,
    atoms: BufWriter<File>,
}

impl Recorder {
    /// Makes the folder `dir` with its three files, and writes setup.csv
    /// for `system`: its number of atoms, box side and cutoff.
    pub fn create(dir: &Path, system: &System) -> io::Result<Self> {
        fs::create_dir_all(dir)?;
        let n = system.positions.len() as f64;
        let mut setup = File::create(dir.join("setup.csv"))?;
        csv::write_header(&mut setup, &SETUP)?;
        csv::write_row(&mut setup, &[n, system.space.side, system.cutoff])?;
        let mut thermo = BufWriter::new(File::create(dir.join("thermo.csv"))?);
        let mut atoms = BufWriter::new(File::create(dir.join("atoms.csv"))?);
        csv::write_header(&mut thermo, &THERMO)?;
        csv::write_header(&mut atoms, &ATOMS)?;
        Ok(Self { thermo, atoms })
    }

    /// Saves `system` as it is now, at `time` in tau: one row of
    /// thermo.csv, which is returned, and one row of atoms.csv per atom.
    pub fn add(&mut self, time: f64, system: &System) -> io::Result<[f64; 5]> {
        let row = measure(time, system);
        csv::write_row(&mut self.thermo, &row)?;
        for (r, v) in system.positions.iter().zip(&system.velocities) {
            csv::write_row(&mut self.atoms, &[r.x, r.y, r.z, v.x, v.y, v.z])?;
        }
        Ok(row)
    }

    /// Writes out whatever is still held in memory. Call it once, last.
    pub fn finish(mut self) -> io::Result<()> {
        self.thermo.flush()?;
        self.atoms.flush()
    }
}