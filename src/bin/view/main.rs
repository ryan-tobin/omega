//! The viewer: a window that replays a saved run (replay.rs).
//!
//! It shows the atoms of one frame at a time, colored by speed
//! (scene.rs), with the frame's time, temperature, potential energy per
//! atom, pressure and total energy over them (hud.rs). Those numbers are
//! not copied from the run's thermo.csv: each frame is measured again
//! from its atoms. On starting, the program prints how many frames give
//! exactly the numbers the run saved.
//!
//! Every pixel of the picture is computed by our own code. Two crates
//! from outside reach the screen (TECH_STACK): winit opens the window
//! and passes on the keys and the mouse, and softbuffer puts our pixels
//! in the window. They are built only with `--features viewer`.
//!
//! Space plays or pauses. The left and right arrows step one frame back
//! or on. Dragging with the mouse turns the box. Escape quits. The first
//! argument names the folder of the run, output/replay if left out.

mod events;
mod viewer;

use omega::replay::Replay;
use softbuffer::Context;
use std::error::Error;
use std::path::Path;
use viewer::Viewer;
use winit::event_loop::EventLoop;

fn main() -> Result<(), Box<dyn Error>> {
    let folder = std::env::args().nth(1);
    let folder = folder.unwrap_or(String::from("output/replay"));
    let replay = Replay::load(Path::new(&folder))
        .map_err(|e| format!("no saved run in {folder} ({e}): run the record program first"))?;
    let (frames, atoms) = (replay.frames.len(), replay.frames[0].positions.len());
    let shown = replay.measure_all();
    let pairs = shown.iter().zip(&replay.saved);
    let equal = pairs.filter(|(a, b)| a == b).count();
    println!("replay of {folder}: {frames} frames of {atoms} atoms");
    println!("{equal} of {frames} frames, measured again, equal thermo.csv exactly");
    println!("space: play or pause   arrows: one frame back or on   drag: turn   esc: quit");

    let event_loop = EventLoop::new()?;
    let context = Context::new(event_loop.owned_display_handle())?;
    let mut viewer = Viewer::new(replay, shown, context);
    event_loop.run_app(&mut viewer)?;
    Ok(())
}