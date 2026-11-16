//! Our numbers next to Rahman's, in his units (rahman.rs has his).
//!
//! The windows the peaks are looked for in are our own choice: in sigma
//! for g(r), first peak 0.9 to 1.4, first trough 1.3 to 1.9, second peak
//! 1.8 to 2.5, and in tau for the lowest point of the VACF, 0 to 1.

use super::{Curves, si};
use omega::peaks;
use omega::rahman;
use omega::run::Run;
use omega::units::ARGON;

/// One row of the table: what is measured, our value, Rahman's. A row
/// without his value ends after ours, with no spaces.
fn row(what: &str, ours: String, his: String) {
    let line = format!("{what:28}{ours:24}{his}");
    println!("{}", line.trim_end());
}

/// Prints the table.
pub fn print(run: &Run, curves: &Curves) {
    let (angstrom, ps, d) = si();
    let g = &curves.g;
    let (r1, g1) = peaks::highest(g, 0.9, 1.4);
    let (r2, g2) = peaks::lowest(g, 1.3, 1.9);
    let (r3, g3) = peaks::highest(g, 1.8, 2.5);
    let zero = peaks::first_below(&curves.normal, 0.0).expect("the VACF never turns negative");
    let (t_low, c_low) = peaks::lowest(&curves.normal, 0.0, 1.0);
    let (t, t_error) = run.mean(|b| b.temperature);
    let (dm, dm_error) = run.mean(|b| b.d_msd);
    let (dv, dv_error) = run.mean(|b| b.d_vacf);
    let kelvin = ARGON.temperature_k();
    let his_d = rahman::DIFFUSION_CM2_PER_S * 1e5;
    let none = String::new;
    println!();
    row("", "ours".to_string(), "Rahman 1964".to_string());
    let ours = format!("{:.2} +- {:.2}", t * kelvin, t_error * kelvin);
    row("temperature, K", ours, "94.4".to_string());
    let ours = format!("{:.2} (g = {g1:.2})", r1 * angstrom);
    row(
        "g(r) first peak, angstrom",
        ours,
        rahman::FIRST_PEAK_ANGSTROM.to_string(),
    );
    let ours = format!("{:.2} (g = {g2:.2})", r2 * angstrom);
    row("g(r) first trough, angstrom", ours, none());
    let ours = format!("{:.2} (g = {g3:.2})", r3 * angstrom);
    row("g(r) second peak, angstrom", ours, none());
    for (i, &(lo, hi)) in rahman::S_PEAK_WINDOWS.iter().enumerate() {
        let (k, height) = peaks::highest(&curves.s, lo, hi);
        let what = format!("S(k) peak {}, k sigma", i + 1);
        let his = format!("{}, X-ray {}", rahman::S_PEAKS[i], rahman::S_PEAKS_XRAY[i]);
        row(&what, format!("{k:.2} (S = {height:.2})"), his);
    }
    let ours = format!("{:.3} +- {:.3}", dm * d, dm_error * d);
    row("D from MSD, 1e-5 cm^2/s", ours, his_d.to_string());
    row(
        "D from VACF, 1e-5 cm^2/s",
        format!("{:.3} +- {:.3}", dv * d, dv_error * d),
        none(),
    );
    let ours = format!("{:.3}", zero * ps);
    row(
        "VACF first below zero, ps",
        ours,
        rahman::VACF_ZERO_PS.to_string(),
    );
    let ours = format!("{:.3} (VACF = {c_low:.3})", t_low * ps);
    row("VACF lowest point, ps", ours, none());
    println!();
    println!(
        "D from MSD is {:+.1}% from Rahman's",
        100.0 * (dm * d / his_d - 1.0)
    );
}
