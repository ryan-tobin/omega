//! Measurements written as CSV text: one line of column names, then one
//! line of numbers per row, separated by commas. And the same text read
//! back.
//!
//! Each number is written by Rust's plain `{}` formatting for f64. It
//! prints the shortest decimal that reads back as exactly the same f64
//! (Rust 1.100.0, library/core/src/fmt/float.rs and the notes at the top
//! of library/core/src/num/imp/flt2dec/mod.rs). So a CSV file holds the
//! measured values exactly, not rounded copies. It never uses an exponent,
//! so a number far from 1, such as 1e-20, prints with many zeros.
//!
//! Reading uses Rust's f64 parser, which returns the f64 closest to the
//! decimal (the notes at the top of library/core/src/num/imp/dec2flt/mod.rs).
//! So a number written here and read back is the same f64, bit for bit.

use std::io::{self, BufRead, Write};

/// A CSV file read back: its column names and all of its numbers.
pub struct Table {
    /// The column names, from the first line.
    pub names: Vec<String>,
    /// Every number, row after row, one per column in each row.
    pub values: Vec<f64>,
}

impl Table {
    /// The rows in file order, each one number per column.
    pub fn rows(&self) -> std::slice::ChunksExact<'_, f64> {
        self.values.chunks_exact(self.names.len())
    }
}

/// Writes the column names, separated by commas, as one line.
pub fn write_header(out: &mut impl Write, names: &[&str]) -> io::Result<()> {
    writeln!(out, "{}", names.join(","))
}

/// Writes one row of numbers, separated by commas, as one line. The
/// program that writes a file says what each column holds, in which unit.
pub fn write_row(out: &mut impl Write, values: &[f64]) -> io::Result<()> {
    let text: Vec<String> = values.iter().map(|v| v.to_string()).collect();
    writeln!(out, "{}", text.join(","))
}

/// Reads CSV text as `write_header` and `write_row` write it. The numbers
/// carry the units they were written in. A line that does not hold one
/// number per column is an error that names the line.
pub fn read(input: impl BufRead) -> io::Result<Table> {
    let mut lines = input.lines();
    let header = lines.next().transpose()?.unwrap_or_default();
    let names: Vec<String> = header.split(',').map(String::from).collect();
    let mut values = Vec::new();
    for (k, line) in lines.enumerate() {
        let bad = || io::Error::other(format!("line {}: not one number per column", k + 2));
        let before = values.len();
        for field in line?.split(',') {
            values.push(field.parse().map_err(|_| bad())?);
        }
        if values.len() - before != names.len() {
            return Err(bad());
        }
    }
    Ok(Table { names, values })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn header_and_rows_are_comma_separated_lines() {
        let mut out = Vec::new();
        write_header(&mut out, &["time", "temperature"]).unwrap();
        write_row(&mut out, &[0.0, 0.7867]).unwrap();
        write_row(&mut out, &[0.05, -1.5]).unwrap();
        let text = String::from_utf8(out).unwrap();
        assert_eq!(text, "time,temperature\n0,0.7867\n0.05,-1.5\n");
    }

    #[test]
    fn numbers_read_back_exactly() {
        // Values with no short decimal form: 0.1 + 0.2, which needs 17
        // significant digits, 1/3, -2/3 and pi, then the smallest and the
        // largest normal f64.
        let values = [
            0.1 + 0.2,
            1.0 / 3.0,
            -2.0 / 3.0,
            std::f64::consts::PI,
            f64::MIN_POSITIVE,
            f64::MAX,
        ];
        let mut out = Vec::new();
        write_row(&mut out, &values).unwrap();
        let text = String::from_utf8(out).unwrap();
        let back: Vec<f64> = text
            .trim_end()
            .split(',')
            .map(|s| s.parse().unwrap())
            .collect();
        assert_eq!(back, values);
    }

    #[test]
    fn a_table_reads_back_as_it_was_written() {
        // Exact: the names, and the same f64s as in the test above.
        let rows = [[0.1 + 0.2, -2.0 / 3.0], [f64::MIN_POSITIVE, f64::MAX]];
        let mut out = Vec::new();
        write_header(&mut out, &["a", "b"]).unwrap();
        for row in &rows {
            write_row(&mut out, row).unwrap();
        }
        let table = read(out.as_slice()).unwrap();
        assert_eq!(table.names, ["a", "b"]);
        assert_eq!(table.rows().collect::<Vec<_>>(), rows);
    }

    #[test]
    fn a_line_with_a_missing_or_broken_number_is_refused() {
        for text in ["a,b\n1,2\n3\n", "a,b\n1,2\n3,x\n", "a,b\n1,2\n3,4,5\n"] {
            let error = read(text.as_bytes()).err().unwrap().to_string();
            assert_eq!(error, "line 3: not one number per column", "{text:?}");
        }
    }
}