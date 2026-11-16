//! The four plots of docs/liquid-argon.md, written as SVG files in the
//! folder `docs`, in Rahman's units: angstrom, picoseconds and
//! 1e-5 cm^2 / s. Each shows our curve in black with Rahman's numbers
//! (rahman.rs) in red.

use super::{Curves, si};
use omega::plot::{self, Axis, Mark, Series};
use omega::rahman::{self, DT};
use omega::run::{self, Run};
use omega::stats;
use std::fs;
use std::io;

/// Writes one plot to docs/`name`.
fn save(name: &str, title: &str, x: Axis, y: Axis, series: &[Series]) -> io::Result<()> {
    fs::write(format!("docs/{name}"), plot::svg(title, &x, &y, series))
}

/// Writes the four plots.
pub fn write(run: &Run, curves: &Curves) -> io::Result<()> {
    fs::create_dir_all("docs")?;
    let (angstrom, ps, d) = si();

    // g(r), with Rahman's first peak as a dashed line.
    let ours: Vec<(f64, f64)> = curves.g.iter().map(|&(r, g)| (r * angstrom, g)).collect();
    let peak = rahman::FIRST_PEAK_ANGSTROM;
    let his = [(peak, 0.0), (peak, 3.5)];
    let series = [
        Series::new("ours", "black", Mark::Line, &ours),
        Series::new("Rahman's first peak", "firebrick", Mark::Dashed, &his),
    ];
    let x = Axis::new("r (angstrom)", (0.0, 16.0), 2.0, 0);
    let y = Axis::new("g(r)", (0.0, 3.5), 0.5, 1);
    save(
        "argon-gr.svg",
        "Radial distribution function g(r)",
        x,
        y,
        &series,
    )?;

    // S(k), with Rahman's peaks as dots on our curve, at the k nearest
    // each of his.
    let nearest = |k: f64| {
        curves
            .s
            .iter()
            .min_by(|a, b| (a.0 - k).abs().total_cmp(&(b.0 - k).abs()))
    };
    let his: Vec<(f64, f64)> = rahman::S_PEAKS
        .iter()
        .filter_map(|&k| nearest(k))
        .copied()
        .collect();
    let series = [
        Series::new("ours", "black", Mark::Line, &curves.s),
        Series::new(
            "Rahman's peaks, on our curve",
            "firebrick",
            Mark::Dots,
            &his,
        ),
    ];
    let x = Axis::new("k sigma", (0.0, 30.0), 5.0, 0);
    let y = Axis::new("S(k)", (0.0, 3.5), 0.5, 1);
    save("argon-sk.svg", "Structure factor S(k)", x, y, &series)?;

    // The VACF, with Rahman's first zero as a dot.
    let ours: Vec<(f64, f64)> = curves.normal.iter().map(|&(t, c)| (t * ps, c)).collect();
    let his = [(rahman::VACF_ZERO_PS, 0.0)];
    let series = [
        Series::new("ours", "black", Mark::Line, &ours),
        Series::new("Rahman's first zero", "firebrick", Mark::Dots, &his),
    ];
    let x = Axis::new("t (ps)", (0.0, 3.0), 0.5, 1);
    let y = Axis::new("VACF / VACF(0)", (-0.2, 1.0), 0.2, 1);
    save(
        "argon-vacf.svg",
        "Velocity autocorrelation function",
        x,
        y,
        &series,
    )?;

    // D: the VACF integrated up to t, which levels off at D, against D
    // from the slope of the MSD and Rahman's D.
    let h = run::SAMPLE_EVERY as f64 * DT;
    let vacf: Vec<f64> = curves.vacf.iter().map(|&(_, v)| v).collect();
    let up_to = |k: usize| {
        (
            k as f64 * h * ps,
            stats::trapezoid(&vacf[..=k], h) / 3.0 * d,
        )
    };
    let mut running = vec![(0.0, 0.0)];
    running.extend((1..vacf.len()).map(up_to));
    let level = |value: f64| [(0.0, value), (11.0, value)];
    let (dm, _) = run.mean(|b| b.d_msd);
    let (msd, his) = (level(dm * d), level(rahman::DIFFUSION_CM2_PER_S * 1e5));
    let series = [
        Series::new(
            "ours, VACF integrated up to t",
            "black",
            Mark::Line,
            &running,
        ),
        Series::new(
            "ours, from the slope of the MSD",
            "steelblue",
            Mark::Dashed,
            &msd,
        ),
        Series::new("Rahman", "firebrick", Mark::Dashed, &his),
    ];
    let x = Axis::new("t (ps)", (0.0, 10.0), 2.0, 0);
    let y = Axis::new("D (1e-5 cm^2/s)", (0.0, 4.0), 0.5, 1);
    save(
        "argon-diffusion.svg",
        "Self-diffusion coefficient D",
        x,
        y,
        &series,
    )
}
