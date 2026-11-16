//! The measured curves and block values written as CSV files in the
//! folder `output/argon`, in reduced units. The measure program's files
//! are in `output`, and these are kept apart from them.
//!
//! Columns: rdf.csv r (sigma) and g; sk.csv k (1 / sigma) and S; msd.csv
//! t (tau) and MSD (sigma^2); vacf.csv t (tau), VACF ((sigma / tau)^2) and
//! VACF / VACF(0); blocks.csv the block number, T (eps / kB), and D from
//! the MSD and from the VACF (sigma^2 / tau).

use super::Curves;
use omega::csv;
use omega::run::Run;
use std::fs::{self, File};
use std::io::{self, BufWriter};

/// Writes `rows` under the column `names` to output/argon/`name`.
fn table(name: &str, names: &[&str], rows: &[Vec<f64>]) -> io::Result<()> {
    let mut out = BufWriter::new(File::create(format!("output/argon/{name}"))?);
    csv::write_header(&mut out, names)?;
    rows.iter()
        .try_for_each(|row| csv::write_row(&mut out, row))
}

/// The points as rows of two numbers.
fn rows(points: &[(f64, f64)]) -> Vec<Vec<f64>> {
    points.iter().map(|&(x, y)| vec![x, y]).collect()
}

/// Writes the five CSV files.
pub fn write(run: &Run, curves: &Curves) -> io::Result<()> {
    fs::create_dir_all("output/argon")?;
    table("rdf.csv", &["r", "g"], &rows(&curves.g))?;
    table("sk.csv", &["k", "s"], &rows(&curves.s))?;
    table("msd.csv", &["t", "msd"], &rows(&curves.msd))?;
    let both = curves.vacf.iter().zip(&curves.normal);
    let both: Vec<Vec<f64>> = both.map(|(&(t, v), &(_, n))| vec![t, v, n]).collect();
    table("vacf.csv", &["t", "vacf", "normalized"], &both)?;
    let blocks = run.blocks.iter().enumerate();
    let blocks = blocks.map(|(i, b)| vec![(i + 1) as f64, b.temperature, b.d_msd, b.d_vacf]);
    let names = ["block", "temperature", "d_msd", "d_vacf"];
    table("blocks.csv", &names, &blocks.collect::<Vec<_>>())
}
