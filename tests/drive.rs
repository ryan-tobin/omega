//! The scenario runner (drive.rs) and its files (record.rs): a run that
//! is stopped, killed and picked up again leaves the same files, byte for
//! byte, as the same run left alone.

use omega::barostat::Barostat;
use omega::checkpoint::Checkpoint;
use omega::drive::{Run, SKIN};
use omega::forces::Method;
use omega::periodic::PeriodicBox;
use omega::scenario::{self, Scenario};
use omega::system::System;
use omega::thermostat::Thermostat;
use omega::{record, start, thermo};
use std::fs;
use std::path::{Path, PathBuf};

/// 500 atoms of liquid argon: 120 steps under the thermostat, then 240
/// with the barostat as well, a move tried every 5 steps.
const TEXT: &str = "cells 5
side 8.524
temperature 1.6
seed 11
cutoff 2.25
dt 0.005
checkpoint 50

stage warm
steps 120
thermostat 0.7867 0.5
sample 10

stage press
steps 240
thermostat 0.7867 0.5
barostat 1.0 0.05 5
sample 10
profile 20
";

/// An empty folder for the files of one run of one test.
fn folder(name: &str) -> PathBuf {
    let folder = Path::new(env!("CARGO_TARGET_TMPDIR")).join(name);
    let _ = fs::remove_dir_all(&folder);
    folder
}

/// Runs the stages of `s` in `folder` up to step `stop`.
fn run(s: &Scenario, folder: &Path, resume: bool, stop: Option<u64>) -> u64 {
    let mut run = Run::new(s, folder, 2, resume).unwrap();
    for index in 0..s.stages.len() {
        run.stage(index, stop).unwrap();
    }
    run.state.step
}

#[test]
fn a_run_stopped_killed_and_picked_up_leaves_the_same_files() {
    let s = scenario::parse(TEXT).unwrap();
    let (whole, parts) = (folder("drive-whole"), folder("drive-parts"));
    assert_eq!(run(&s, &whole, false, None), 360);
    // Stopped in stage 0, picked up, stopped in stage 1.
    assert_eq!(run(&s, &parts, false, Some(70)), 70);
    assert_eq!(run(&s, &parts, true, Some(200)), 200);
    // Killed at step 215: the checkpoint on disk is the one of step 200,
    // thermo.csv already has the row of step 210, and the row of step 220
    // is half written.
    let saved = fs::read(parts.join("checkpoint.bin")).unwrap();
    assert_eq!(run(&s, &parts, true, Some(215)), 215);
    fs::write(parts.join("checkpoint.bin"), saved).unwrap();
    let mut thermo = fs::read_to_string(parts.join("thermo.csv")).unwrap();
    assert!(thermo.contains("\n210,"), "the row of step 210 is there");
    thermo.push_str("220,1.1,0.78");
    fs::write(parts.join("thermo.csv"), thermo).unwrap();
    assert_eq!(run(&s, &parts, true, None), 360);
    for name in ["thermo.csv", "profile.csv", "checkpoint.bin"] {
        let (a, b) = (fs::read(whole.join(name)), fs::read(parts.join(name)));
        assert!(a.unwrap() == b.unwrap(), "{name} differs");
    }
    // The files hold what the scenario asks for: 12 + 24 rows of thermo.csv
    // and 12 of profile.csv, and the barostat's 48 tries in the checkpoint.
    let rows = record::rows(&whole.join("thermo.csv"), 0, 360).unwrap();
    assert_eq!((rows.len(), rows[0][0], rows[35][0]), (36, 10.0, 360.0));
    let profile = record::rows(&whole.join("profile.csv"), 0, 360).unwrap();
    let state = record::load(&whole).unwrap();
    assert_eq!((profile.len(), state.step, state.barostat.0), (12, 360, 48));
    assert!(state.barostat.1 > 10 && state.side != 8.524);
}

#[test]
fn stage_starts_at_its_thermostat_s_temperature() {
    // The atoms start at 1.6 and the stage asks for 0.7867. With a row
    // every step, the first row, after one step, must be near 0.7867:
    // measured 0.7859. Tolerance: 2%. The row must be in the file as soon
    // as the checkpoint is, while the run is still open.
    let s = scenario::parse(&TEXT.replacen("sample 10", "sample 1", 1)).unwrap();
    let folder = folder("drive-start");
    let mut run = Run::new(&s, &folder, 2, false).unwrap();
    run.stage(0, Some(1)).unwrap();
    let first = &record::rows(&folder.join("thermo.csv"), 0, 1).unwrap()[0];
    assert!((first[2] / 0.7867 - 1.0).abs() <= 0.02, "T = {}", first[2]);
    assert_eq!(run.state.step, 1);
}

#[test]
fn averages_come_from_five_blocks_of_rows() {
    // The step column of the rows after step 120: 130, 140, ... 360. Five
    // blocks of 4 rows have means 145, 185, 225, 265 and 305, so the mean
    // is 225 and its error sqrt(4000 / 5) = 28.28. The last 4 rows are
    // left over. And a checkpoint of another run is refused.
    let s = scenario::parse(TEXT).unwrap();
    let folder = folder("drive-rows");
    run(&s, &folder, false, None);
    let rows = record::rows(&folder.join("thermo.csv"), 120, 360).unwrap();
    assert_eq!(rows.len(), 24);
    let (mean, error) = record::average(&rows, 0);
    assert_eq!(mean, 225.0);
    assert!((error - 800f64.sqrt()).abs() <= 1e-12, "error {error}");
    let other = scenario::parse(&TEXT.replace("cells 5", "cells 4")).unwrap();
    assert!(Run::new(&other, &folder, 2, true).is_err());
}

#[test]
fn the_runner_steps_as_the_library_does_by_hand() {
    // 7 free steps, then 20 held, with a row and a barostat move (drawn
    // from seed + 2 in stage 1) every 5 steps of the stage: after steps
    // 12, 17, 22 and 27. The same steps by hand end on the same bits.
    let text = "cells 4\nside 6.8192\ntemperature 1.6\nseed 5\ncutoff 2.25\ndt 0.005\n
stage free\nsteps 7\n\nstage held\nsteps 20\nthermostat 0.7867 0.5\nbarostat 3 0.005 5\nsample 5\n";
    let s = scenario::parse(text).unwrap();
    let folder = folder("drive-hand");
    run(&s, &folder, false, None);
    let rows = record::rows(&folder.join("thermo.csv"), 0, 27).unwrap();
    assert_eq!((rows.len(), rows[0][0], rows[3][0]), (4, 12.0, 27.0));
    let ((r, v), space) = (start::atoms(&s), PeriodicBox { side: s.side });
    let mut system = System::with_forces(r, v, space, s.cutoff, Method::listed(SKIN, 2));
    (0..7).for_each(|_| system.step(s.dt));
    thermo::rescale(&mut system.velocities, 0.7867);
    let mut thermostat = Thermostat::new(0.7867, 0.5);
    let mut barostat = Barostat::new(3.0, 0.7867, 0.005, 5 + 2, 0, 0);
    for k in 1..=20 {
        thermostat.half_step(&mut system.velocities, s.dt);
        system.step(s.dt);
        thermostat.half_step(&mut system.velocities, s.dt);
        if k % 5 == 0 {
            barostat.attempt(&mut system, Method::listed(SKIN, 2));
        }
    }
    let by_hand = Checkpoint::of(&system, 27, Some(&thermostat), Some(&barostat));
    assert!(record::load(&folder).unwrap() == by_hand, "the two differ");
}
