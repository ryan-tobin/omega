//! Runs a scenario file (scenario.rs) and prints what each stage measured.
//!
//!   cargo run --release --bin scenario -- scenarios/boil.txt
//!
//! The files of the run go to output/<name of the scenario file>
//! (record.rs). With `--resume` after the file name the run carries on from
//! the checkpoint in that folder. With `--stop <step>` it stops after that
//! step of the run.

use omega::drive::Run;
use omega::profile;
use omega::record::{self, SLICES, THERMO};
use omega::scenario;
use std::path::Path;
use std::time::Instant;

fn main() {
    if let Err(message) = run() {
        eprintln!("scenario: {message}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), String> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let file = args.first().ok_or("give a scenario file")?;
    let resume = args.iter().any(|a| a == "--resume");
    let stop = match args.iter().position(|a| a == "--stop") {
        Some(i) => Some(
            args.get(i + 1)
                .and_then(|n| n.parse::<u64>().ok())
                .ok_or("--stop needs a step")?,
        ),
        None => None,
    };
    let text = std::fs::read_to_string(file).map_err(|e| format!("{file}: {e}"))?;
    let s = scenario::parse(&text).map_err(|e| format!("{file}: {e}"))?;
    let name = Path::new(file).file_stem().ok_or("the file has no name")?;
    let folder = Path::new("output").join(name);
    let failed = |e: std::io::Error| format!("{}: {e}", folder.display());
    let threads = std::thread::available_parallelism().map_or(1, |n| n.get());
    let mut run = Run::new(&s, &folder, threads, resume).map_err(failed)?;
    let (atoms, from) = (run.system.positions.len(), run.state.step);
    println!(
        "{atoms} atoms on {threads} threads, from step {from}, files in {}",
        folder.display()
    );
    let mut first = 0;
    for (index, stage) in s.stages.iter().enumerate() {
        let clock = Instant::now();
        run.stage(index, stop).map_err(failed)?;
        let last = first + stage.steps as u64;
        let seconds = clock.elapsed().as_secs_f64();
        println!(
            "stage {}, steps {} to {last}: {seconds:.1} s",
            stage.name,
            first + 1
        );
        let rows = record::rows(&folder.join("thermo.csv"), first, last).map_err(failed)?;
        if rows.len() >= 10 {
            for (column, name) in THERMO.iter().enumerate().skip(2) {
                let (mean, error) = record::average(&rows, column);
                println!("  {name:<11} {mean:8.4} +- {error:.4}");
            }
        }
        let rows = record::rows(&folder.join("profile.csv"), first, last).map_err(failed)?;
        if !rows.is_empty() {
            let mean = |k: usize| rows.iter().map(|r| r[k]).sum::<f64>() / rows.len() as f64;
            let slices: Vec<f64> = (1..=SLICES).map(mean).collect();
            let (middle, ends) = profile::plateaus(&slices);
            println!("  density in the middle {middle:.4}, at the ends {ends:.4}");
        }
        let (tries, kept) = run.state.barostat;
        if stage.barostat.is_some() && run.state.step > first && run.state.step <= last {
            println!("  barostat kept {kept} of {tries} moves");
        }
        first = last;
    }
    if run.state.step < first {
        println!(
            "stopped after step {}: carry on with --resume",
            run.state.step
        );
    }
    Ok(())
}
