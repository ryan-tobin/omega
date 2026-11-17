//! The proof of M6: what the viewer shows is the run that was saved.
//!
//! Each test saves a short run of 108 atoms, started as in
//! tests/conservation.rs, to a folder of its own, reads it back and
//! removes the folder. The folders are made in target/tmp, the place
//! Cargo names in CARGO_TARGET_TMPDIR for the files of tests like these.

use omega::hud;
use omega::lattice;
use omega::periodic::PeriodicBox;
use omega::picture::Picture;
use omega::recorder::Recorder;
use omega::replay::Replay;
use omega::rng::Rng;
use omega::scene::{self, Camera};
use omega::system::System;
use omega::thermo;
use omega::velocities;
use std::fs;
use std::path::{Path, PathBuf};

const DT: f64 = 0.005;
const BACKGROUND: u32 = 0x0010_1418;

/// Saves 5 moments, 20 steps apart, of a run of 108 atoms to the folder
/// `name`. Returns the folder, the rows as the run measured them, and
/// the system as it was at the last saved moment.
fn saved_run(name: &str) -> (PathBuf, Vec<[f64; 5]>, System) {
    let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join(name);
    let space = PeriodicBox { side: 10.229 / 2.0 };
    let positions = lattice::fcc(3, space.side);
    let velocities = velocities::maxwell_boltzmann(positions.len(), 0.7867, &mut Rng::new(6));
    let mut system = System::new(positions, velocities, space, 2.25);
    let mut recorder = Recorder::create(&dir, &system).unwrap();
    let mut rows = Vec::new();
    for moment in 0..5 {
        (0..20 * moment.min(1)).for_each(|_| system.step(DT));
        rows.push(recorder.add(moment as f64 * 20.0 * DT, &system).unwrap());
    }
    recorder.finish().unwrap();
    (dir, rows, system)
}

/// What `Replay::load` says against the folder `dir`.
fn refusal(dir: &Path) -> String {
    let error = Replay::load(dir).err();
    error.expect("the run was read").to_string()
}

#[test]
fn every_frame_measured_again_equals_the_row_that_was_saved() {
    // Exact, with no tolerance: the files hold every number exactly, and
    // measuring a frame again does the same sums in the same order
    // (replay.rs). The atoms have left the box by the last frame, so
    // positions that were never wrapped are part of what is tested.
    let (dir, rows, system) = saved_run("replay-equal");
    let replay = Replay::load(&dir).unwrap();
    let last = &replay.frames[4];
    assert_eq!(replay.frames.len(), 5);
    assert!(last.positions == system.positions && last.velocities == system.velocities);
    assert!(last.positions.iter().any(|r| r.x < 0.0));
    assert_eq!(replay.saved, rows);
    assert_eq!(replay.measure_all(), rows);
    assert_ne!(rows[3], rows[4]);
    // The columns of a row: time, T, U / N, P and E.
    let v = &system.velocities;
    let t = thermo::temperature(v);
    let p = thermo::pressure(v, system.virial, system.space);
    let (u, e) = (system.potential / 108.0, system.total_energy());
    assert_eq!(rows[4], [4.0 * 20.0 * DT, t, u, p, e]);
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn the_numbers_read_off_the_picture_are_the_numbers_in_the_file() {
    // Each frame is drawn as the viewer draws it: the atoms, then the
    // text over them. In this small picture the text lies over more than
    // 1000 pixels of atoms (measured: 1093 to 2338). The text is read
    // back from the pixels and must be the row of thermo.csv to the
    // digits shown. Tolerance: half of the last digit shown, 0.005 for
    // the time and 0.00005 for the rest, plus 1e-9 for the rounding of
    // the comparison itself.
    let (dir, _, _) = saved_run("replay-picture");
    let replay = Replay::load(&dir).unwrap();
    let full = scene::full_speed(&replay.frames[0].velocities);
    let (yaw, pitch) = (0.6, 0.4);
    let camera = Camera { yaw, pitch };
    let labels = ["TIME ", "T    ", "U/N  ", "P    ", "E    "];
    for (k, frame) in replay.frames.iter().enumerate() {
        let mut picture = Picture::new(400, 300, BACKGROUND);
        scene::draw(&mut picture, frame, replay.space, &camera, full);
        hud::draw(&mut picture, &hud::lines(replay.measure(k), false), 2);
        for (line, label) in labels.iter().enumerate() {
            let text = hud::read(&picture, line, 15, 2);
            let shown: f64 = text[5..].trim().parse().expect("a number");
            let half = if line == 0 { 0.005 } else { 0.00005 };
            let miss = (shown - replay.saved[k][line]).abs();
            assert!(text.starts_with(label), "frame {k}: {text}");
            assert!(miss <= half + 1e-9, "frame {k}: {text}");
        }
    }
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn a_frame_with_one_atom_changed_no_longer_gives_the_saved_row() {
    // The numbers shown come from the atoms, not from thermo.csv: moving
    // one atom by a millionth of sigma changes the potential energy, and
    // making one a millionth faster changes the temperature. Only the
    // time is copied.
    let (dir, rows, _) = saved_run("replay-changed");
    let mut replay = Replay::load(&dir).unwrap();
    replay.frames[2].positions[17].x += 1e-6;
    replay.frames[3].velocities[17] = replay.frames[3].velocities[17] * (1.0 + 1e-6);
    let (moved, faster) = (replay.measure(2), replay.measure(3));
    assert_eq!((moved[0], moved[1]), (rows[2][0], rows[2][1]));
    assert_ne!(moved[2], rows[2][2]);
    assert_ne!(faster[1], rows[3][1]);
    assert_eq!(faster[2], rows[3][2]);
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn files_that_are_not_one_saved_run_are_refused() {
    let (dir, _, _) = saved_run("replay-refused");
    let atoms = fs::read_to_string(dir.join("atoms.csv")).unwrap();
    let thermo = fs::read_to_string(dir.join("thermo.csv")).unwrap();

    // One atom of the last frame is missing, then there is one too many.
    let short = atoms.trim_end().rsplit_once('\n').unwrap().0;
    for wrong in [format!("{short}\n"), format!("{atoms}0,0,0,0,0,0\n")] {
        fs::write(dir.join("atoms.csv"), wrong).unwrap();
        assert_eq!(refusal(&dir), "atoms.csv does not fit thermo.csv");
    }
    fs::write(dir.join("atoms.csv"), &atoms).unwrap();

    // No frames at all.
    let header = thermo.lines().next().unwrap();
    fs::write(dir.join("thermo.csv"), format!("{header}\n")).unwrap();
    assert_eq!(refusal(&dir), "a saved run with no atoms or no frames");

    // Another table of five columns under the name thermo.csv.
    fs::write(dir.join("thermo.csv"), "a,b,c,d,e\n0,1,2,3,4\n").unwrap();
    assert_eq!(refusal(&dir), "not the files of a saved run");

    // No folder.
    fs::remove_dir_all(&dir).unwrap();
    assert!(Replay::load(&dir).is_err());
}