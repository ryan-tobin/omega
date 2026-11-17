//! A saved run read back: what the viewer replays.
//!
//! recorder.rs saves a run as three CSV files and says what they hold.
//! Reading them back makes a frame of each saved moment. Measuring a
//! frame again builds a System from its atoms alone and takes the same
//! four measurements the run took. csv.rs writes and reads every number
//! exactly, so the atoms are where the run had them, to the last bit.
//! The forces come from all pairs, in the same order as in a run that
//! took all pairs too, so the measurements come out equal to the saved
//! row, also to the last bit (tested in tests/replay.rs). That is how the
//! viewer shows that what it draws is the run that was saved.
//!
//! Taking all pairs suits runs of a few thousand atoms, not a million.
//! A run on the neighbor list (neighbors.rs) adds the same pairs in
//! another order, so its saved rows could differ from these in the last
//! digits.
//!
//! Label: derived. Nothing is modeled here: a frame is the run's own state.

use crate::csv;
use crate::periodic::PeriodicBox;
use crate::recorder;
use crate::system::System;
use crate::vec3::Vec3;
use std::fs::File;
use std::io::{self, BufReader};
use std::path::Path;

/// One frame, a snapshot of the run at a saved moment: the position of
/// every atom, in sigma and never wrapped, and its velocity, in
/// sigma / tau.
pub struct Snapshot {
    pub positions: Vec<Vec3>,
    pub velocities: Vec<Vec3>,
}

/// A saved run read back.
pub struct Replay {
    /// The repeating box of the run.
    pub space: PeriodicBox,
    /// its force cutoff, in sigma.
    pub cutoff: f64,
    /// The rows of thermo.csv as the run saved them, one per frame.
    pub saved: Vec<[f64; 5]>,
    /// The frames, in the order they were saved.
    pub frames: Vec<Snapshot>,
}

/// Columns `from`, `from + 1` and `from + 2` of rows of atoms.csv, as the
/// x, y and z of one vector per row.
fn vectors(rows: &[f64], from: usize) -> Vec<Vec3> {
    let vector = |a: &[f64]| Vec3::new(a[from], a[from + 1], a[from + 2]);
    rows.chunks(6).map(vector).collect()
}

impl Replay {
    /// Reads the run saved in the folder `dir`. Files with other columns
    /// than recorder.rs writes, a run with no atoms or no frames, and an
    /// atoms.csv without one row per atom for every row of thermo.csv
    /// are errors.
    pub fn load(dir: &Path) -> io::Result<Self> {
        let read = |name: &str| csv::read(BufReader::new(File::open(dir.join(name))?));
        let (setup, thermo, atoms) = (read("setup.csv")?, read("thermo.csv")?, read("atoms.csv")?);
        let named = setup.names == recorder::SETUP
            && thermo.names == recorder::THERMO 
            && atoms.names == recorder::ATOMS;
        if !named || setup.values.len() != 3 {
            return Err(io::Error::other("not the files of a saved run"));
        }
        let (n, side, cutoff) = (setup.values[0] as usize, setup.values[1], setup.values[2]);
        let saved: Vec<[f64; 5]> = thermo.rows().map(|row| row.try_into().unwrap()).collect();
        if n == 0 || saved.is_empty() {
            return Err(io::Error::other("a saved run with no atoms or no frames"));
        }
        if atoms.values.len() != 6 * n * saved.len() {
            return Err(io::Error::other("atoms.csv does not fit thermo.csv"));
        }
        let frames = atoms.values.chunks(6 * n).map(|rows| Snapshot {
            positions: vectors(rows, 0),
            velocities: vectors(rows, 3),
        });
        Ok(Self {
            space: PeriodicBox { side },
            cutoff,
            saved,
            frames: frames.collect(),
        })
    }

    /// Frame `k` measured again from its atoms alone: the time in tau,
    /// copied from thermo.csv, then the temperature, potential energy
    /// per atom, pressure and total energy, as recorder::measure gives
    /// them.
    pub fn measure(&self, k: usize) -> [f64; 5] {
        let frame = &self.frames[k];
        let (r, v) = (frame.positions.clone(), frame.velocities.clone());
        let system = System::new(r, v, self.space, self.cutoff);
        recorder::measure(self.saved[k][0], &system)
    }

    /// Every frame measured again, in order: the numbers the view shows.
    pub fn measure_all(&self) -> Vec<[f64; 5]> {
        (0..self.frames.len()).map(|k| self.measure(k)).collect()
    }
}

