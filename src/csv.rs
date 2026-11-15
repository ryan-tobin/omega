//! Measurements written as CSV text: one line of column names, then one
//! line of numbers per row, separated by commas.
//!
//! Each number is written by Rust's plain `{}` formatting for f64. It
//! prints the shortest decimal that reads back as exactly the same f64
//! (Rust 1.100.0, library/core/src/fmt/float.rs and the notes at the top
//! of library/core/src/num/imp/flt2dec/mod.rs). So a CSV file holds the
//! measured values exactly, not rounded copies. It never uses an exponent,
//! so a number far from 1, such as 1e-20, prints with many zeros.

use std::io::{self, Write};

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
}
