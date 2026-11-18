//! A run that is stopped, saved as bytes and picked up again must end
//! where the same run ends when it never stops: to the last bit.
//!
//! Each run is 400 steps of liquid argon (tests/common/mod.rs) under the
//! thermostat. With the barostat, at a pressure of 1, a move is tried
//! every 5 steps. The stopped run stops after steps 153 and 290, turns
//! into the bytes of a checkpoint each time, and is made again from those
//! bytes alone. The two runs are compared through their last checkpoints,
//! byte for byte: step, box, positions, velocities, thermostat, barostat
//! counts, and the positions the neighbor list was last built at.

mod common;

use common::{T, liquid, step};
use omega::barostat::Barostat;
use omega::checkpoint::Checkpoint;
use omega::forces::Method;
use omega::system::System;
use omega::thermostat::Thermostat;

/// What a test changes: forces from a neighbor list on this many threads
/// (0: all pairs), a barostat or not, and where to stop on the way.
struct Plan {
    threads: usize,
    barostat: bool,
    stops: &'static [u64],
}

impl Plan {
    fn method(&self) -> Method {
        match self.threads {
            0 => Method::AllPairs,
            threads => Method::listed(0.3, threads),
        }
    }

    /// The checkpoint after 400 steps, and the barostat moves kept.
    fn run(&self) -> (Checkpoint, u64) {
        let (liquid, mut thermostat) = liquid(1);
        let (r, v) = (liquid.positions, liquid.velocities);
        let mut system = System::with_forces(r, v, liquid.space, 2.25, self.method());
        let new = |tries, kept| Barostat::new(1.0, T, 0.05, 7, tries, kept);
        let mut barostat = self.barostat.then(|| new(0, 0));
        let mut done = 0;
        for &stop in self.stops.iter().chain(&[400]) {
            while done < stop {
                step(&mut system, &mut thermostat);
                done += 1;
                match &mut barostat {
                    Some(b) if done % 5 == 0 => b.attempt(&mut system, self.method()),
                    _ => false,
                };
            }
            // Stop: only the bytes are kept, and the run is made from them.
            let state = Checkpoint::of(&system, done, Some(&thermostat), barostat.as_ref());
            if self.threads > 0 && stop < 400 {
                let (_, built_at) = state.list.as_ref().expect("a neighbor list");
                assert!(
                    *built_at != state.positions,
                    "atoms moved since the list was built"
                );
            }
            let state = Checkpoint::from_bytes(&state.to_bytes()).unwrap();
            system = state.system(self.threads);
            thermostat = Thermostat::new(T, 0.5);
            (thermostat.friction, thermostat.eta) = state.thermostat;
            barostat = self
                .barostat
                .then(|| new(state.barostat.0, state.barostat.1));
        }
        let kept = barostat.as_ref().map_or(0, |b| b.accepted);
        (
            Checkpoint::of(&system, done, Some(&thermostat), barostat.as_ref()),
            kept,
        )
    }
}

/// Runs `threads` and `barostat` without a stop and with two, and checks
/// that the two runs end the same.
fn stopped_run_matches(threads: usize, barostat: bool) {
    let (whole, kept) = Plan {
        threads,
        barostat,
        stops: &[],
    }
    .run();
    let (parts, _) = Plan {
        threads,
        barostat,
        stops: &[153, 290],
    }
    .run();
    assert_eq!(whole.step, 400);
    assert!(
        whole.to_bytes() == parts.to_bytes(),
        "the stopped run ended elsewhere"
    );
    assert!(
        !barostat || (20..=70).contains(&kept),
        "{kept} of 80 moves kept"
    );
}

#[test]
fn all_pairs_under_the_thermostat() {
    stopped_run_matches(0, false);
}

#[test]
fn neighbor_list_under_the_thermostat() {
    stopped_run_matches(1, false);
}

#[test]
fn neighbor_list_under_thermostat_and_barostat() {
    stopped_run_matches(2, true);
}

#[test]
fn all_pairs_under_thermostat_and_barostat() {
    stopped_run_matches(0, true);
}

#[test]
fn a_list_built_where_the_run_stopped_is_not_enough() {
    // Why the checkpoint keeps the positions of the last build: a list
    // built at the stop holds the same pairs in another order, the forces
    // differ in their last bits, and the runs drift apart.
    let (whole, _) = Plan {
        threads: 1,
        barostat: false,
        stops: &[],
    }
    .run();
    let (liquid, mut thermostat) = liquid(1);
    let (r, v) = (liquid.positions, liquid.velocities);
    let mut system = System::with_forces(r, v, liquid.space, 2.25, Method::listed(0.3, 1));
    for done in 1..=400 {
        step(&mut system, &mut thermostat);
        if done == 153 {
            let mut state = Checkpoint::of(&system, done, None, None);
            state.list = Some((0.3, state.positions.clone()));
            system = state.system(1);
        }
    }
    assert!(
        system.positions != whole.positions,
        "the order did not matter"
    );
}
