//! The NIST Lennard-Jones benchmark (nist.rs): energy and pressure of 500
//! atoms at T = 0.85, density 0.776, cut at 3 sigma with long-range
//! corrections, next to NIST's Monte Carlo and molecular dynamics rows.
//!
//! cargo run --release --bin nist -- [blocks]
//!
//! runs `blocks` blocks of 100 tau (default 20), seed 1, on every thread
//! the machine offers.

use omega::nist::{self, Averages};
use omega::stats;

fn main() {
    let blocks = std::env::args()
        .nth(1)
        .map_or(20, |a| a.parse().expect("a number"));
    let threads = std::thread::available_parallelism().map_or(1, |n| n.get());
    let runs = nist::run(1, blocks, threads);
    println!("500 atoms, {blocks} blocks of 100 tau, {threads} threads");
    println!("                         our states        weighed as the cut potential");
    println!("block  temperature  energy    pressure  energy    pressure");
    for (b, a) in runs.iter().enumerate() {
        let (t, u, p) = (a.temperature, a.energy, a.pressure);
        let (cu, cp) = (a.cut_energy, a.cut_pressure);
        println!("{b:5}  {t:11.5}  {u:.5}  {p:8.5}  {cu:.5}  {cp:8.5}");
    }
    let mean =
        |f: fn(&Averages) -> f64| stats::mean_and_error(&runs.iter().map(f).collect::<Vec<f64>>());
    let (t, t_err) = mean(|a| a.temperature);
    let (u, u_err) = mean(|a| a.energy);
    let (p, p_err) = mean(|a| a.pressure);
    let (cu, cu_err) = mean(|a| a.cut_energy);
    let (cp, cp_err) = mean(|a| a.cut_pressure);
    println!("temperature {t:.5} +- {t_err:.5}");
    println!("our states:          energy {u:.5} +- {u_err:.5}, pressure {p:.5} +- {p_err:.5}");
    println!("weighed as cut:      energy {cu:.5} +- {cu_err:.5}, pressure {cp:.5} +- {cp_err:.5}");
    let ((mu, mu_err), (mp, mp_err)) = (nist::MC_ENERGY, nist::MC_PRESSURE);
    println!("NIST Monte Carlo:    energy {mu:.5} +- {mu_err:.5}, pressure {mp:.5} +- {mp_err:.5}");
    let [(dt, _), (du, _), (dp, _)] = nist::MD_ROW;
    println!("NIST MD (at T {dt}): energy {du}, pressure {dp}");
}
