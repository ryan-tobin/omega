//! Saves a short run of Rahman's liquid for the viewer to replay, then
//! reads it back and checks it the way the viewer does.
//!
//! The run is the start of the measure program's: 864 atoms melted and
//! brought to Rahman's temperature (rahman.rs), then left to move freely.
//! Every SAVE_EVERY steps the atoms and their measurements are saved to
//! the folder output/replay (recorder.rs). The rows printed here are rows
//! of the measure program's table.
//!
//! Reading the run back measures every frame again from its atoms alone
//! (replay.rs). Each must give exactly the numbers that were saved.

use omega::rahman::{self, CUTOFF, DT, RESCALE_STEPS, SIDE, TEMPERATURE};
use omega::recorder::Recorder;
use omega::replay::Replay;
use std::io;
use std::path::Path;

// The seed of M2 to M5. Our own choices: 2000 free steps (10 tau), saved
// every 10 steps (0.05 tau), which is 201 frames. Between two frames an
// atom at the root-mean-square speed, 1.54 sigma / tau, moves 0.08 sigma.
const SEED: u64 = 1964;
const STEPS: usize = 2000;
const SAVE_EVERY: usize = 10;
const FOLDER: &str = "output/replay";

fn main() -> io::Result<()> {
    let mut system = rahman::liquid(SEED);
    let atoms = system.positions.len();
    let mut recorder = Recorder::create(Path::new(FOLDER), &system)?;
    println!("{atoms} atoms, box side {SIDE}, cutoff {CUTOFF}, dt {DT}, seed {SEED}");
    println!("rescaled to T = {TEMPERATURE} for {RESCALE_STEPS} steps, then free");
    println!("  time       T     U/N        P          E");
    for step in 0..=STEPS {
        if step > 0 {
            system.step(DT);
        }
        if step % SAVE_EVERY != 0 {
            continue;
        }
        let [time, t, u, p, e] = recorder.add(step as f64 * DT, &system)?;
        if step % 1000 == 0 {
            println!("{time:6.1} {t:7.4} {u:7.4} {p:8.4} {e:10.4}");
        }
    }
    recorder.finish()?;

    let replay = Replay::load(Path::new(FOLDER))?;
    let frames = replay.frames.len();
    let shown = replay.measure_all();
    let pairs = shown.iter().zip(&replay.saved);
    let equal = pairs.filter(|(a, b)| a == b).count();
    println!();
    println!("saved {frames} frames of {atoms} atoms to {FOLDER}");
    println!("read back: {equal} of {frames} frames, measured again, equal thermo.csv exactly");
    Ok(())
}