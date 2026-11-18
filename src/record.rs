//! The files a run keeps in its folder.
//!
//! thermo.csv has a row for every sample: steps done, time in tau,
//! temperature in eps / kB, pressure in eps / sigma^3, potential energy
//! per atom in eps, and number density in 1 / sigma^3. profile.csv has the
//! density along x (profile.rs): steps done, then one column per slice, in
//! 1 / sigma^3. checkpoint.bin is the last saved state (checkpoint.rs).
//!
//! A run that is picked up again keeps the rows up to its checkpoint and
//! drops the rows after it, which the steps to come write again. So the
//! files of a run that was stopped and picked up are the same, byte for
//! byte, as those of a run that never stopped (tested in tests/drive.rs).

use crate::checkpoint::Checkpoint;
use crate::csv;
use crate::profile;
use crate::stats;
use crate::system::System;
use crate::thermo;
use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

/// The columns of thermo.csv.
pub const THERMO: [&str; 6] = [
    "step",
    "time",
    "temperature",
    "pressure",
    "potential",
    "density",
];

/// The slices of the density profile. Our own choice.
pub const SLICES: usize = 60;

/// The open files of a run.
pub struct Record {
    folder: PathBuf,
    thermo: BufWriter<File>,
    profile: BufWriter<File>,
}

impl Record {
    /// Opens the files in `folder` and keeps the rows they hold up to step
    /// `kept`. 0 keeps none.
    pub fn open(folder: &Path, kept: u64) -> io::Result<Self> {
        fs::create_dir_all(folder)?;
        let open = |name: &str, header: &[&str]| -> io::Result<BufWriter<File>> {
            let path = folder.join(name);
            let old = fs::read_to_string(&path).unwrap_or_default();
            let mut file = BufWriter::new(File::create(&path)?);
            csv::write_header(&mut file, header)?;
            for line in old.lines().skip(1) {
                let step = line.split(',').next().and_then(|w| w.parse::<f64>().ok());
                if step.is_some_and(|step| step <= kept as f64) {
                    writeln!(file, "{line}")?;
                }
            }
            Ok(file)
        };
        let names: Vec<String> = (0..SLICES).map(|k| format!("slice{k}")).collect();
        let mut columns = vec!["step"];
        columns.extend(names.iter().map(String::as_str));
        let (thermo, profile) = (open("thermo.csv", &THERMO)?, open("profile.csv", &columns)?);
        let folder = folder.to_path_buf();
        Ok(Self {
            folder,
            thermo,
            profile,
        })
    }

    /// Adds a row to thermo.csv: `system` after `step` steps of `dt`.
    pub fn sample(&mut self, step: u64, dt: f64, system: &System) -> io::Result<()> {
        let (n, side) = (system.positions.len() as f64, system.space.side);
        let t = thermo::temperature(&system.velocities);
        let p = thermo::pressure(&system.velocities, system.virial, system.space);
        let (u, density) = (system.potential / n, n / (side * side * side));
        csv::write_row(
            &mut self.thermo,
            &[step as f64, step as f64 * dt, t, p, u, density],
        )
    }

    /// Adds a row to profile.csv: `system` after `step` steps.
    pub fn slices(&mut self, step: u64, system: &System) -> io::Result<()> {
        let mut row = vec![step as f64];
        row.extend(profile::along_x(&system.positions, system.space, SLICES));
        csv::write_row(&mut self.profile, &row)
    }

    /// Writes out the rows so far and saves `state` as checkpoint.bin. The
    /// bytes go to a second file first, which then takes the name, so a
    /// run that is killed half way through leaves the old checkpoint whole.
    pub fn save(&mut self, state: &Checkpoint) -> io::Result<()> {
        self.thermo.flush()?;
        self.profile.flush()?;
        let new = self.folder.join("checkpoint.new");
        fs::write(&new, state.to_bytes())?;
        fs::rename(new, self.folder.join("checkpoint.bin"))
    }
}

/// The checkpoint saved in `folder`.
pub fn load(folder: &Path) -> io::Result<Checkpoint> {
    let bytes = fs::read(folder.join("checkpoint.bin"))?;
    Checkpoint::from_bytes(&bytes).map_err(io::Error::other)
}

/// The rows of the CSV file at `path` whose first column, the step, is
/// above `first` and at most `last`.
pub fn rows(path: &Path, first: u64, last: u64) -> io::Result<Vec<Vec<f64>>> {
    let mut rows = Vec::new();
    for line in fs::read_to_string(path)?.lines().skip(1) {
        let row: Result<Vec<f64>, _> = line.split(',').map(str::parse).collect();
        let row = row.map_err(io::Error::other)?;
        if row[0] > first as f64 && row[0] <= last as f64 {
            rows.push(row);
        }
    }
    Ok(rows)
}

/// The mean of column `column` of `rows`, and its standard error from the
/// means of five blocks of rows in a row (stats.rs). Needs 10 rows or more.
pub fn average(rows: &[Vec<f64>], column: usize) -> (f64, f64) {
    let mean =
        |block: &[Vec<f64>]| block.iter().map(|r| r[column]).sum::<f64>() / block.len() as f64;
    let blocks: Vec<f64> = rows
        .chunks_exact(rows.len() / 5)
        .take(5)
        .map(mean)
        .collect();
    stats::mean_and_error(&blocks)
}
