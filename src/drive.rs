//! Runs a scenario (scenario.rs): steps its starting atoms (start.rs)
//! through its stages, and keeps its files (record.rs).
//!
//! One step of a stage is half a step of the thermostat, a velocity Verlet
//! step, and the other half step of the thermostat. Then, when they are
//! due, a try of the barostat, a row of thermo.csv, a row of profile.csv
//! and a checkpoint.
//!
//! A stage starts with a new thermostat and barostat, and with the
//! velocities scaled to the thermostat's temperature (thermo::rescale). A
//! Nose-Hoover thermostat that starts far from its temperature swings
//! around it for a long time: started three times too hot, the atoms cool
//! to about a fifth of the target before they come back.

use crate::barostat::Barostat;
use crate::checkpoint::Checkpoint;
use crate::forces::Method;
use crate::periodic::PeriodicBox;
use crate::record::{self, Record};
use crate::scenario::Scenario;
use crate::start;
use crate::system::System;
use crate::thermo;
use crate::thermostat::Thermostat;
use std::io;
use std::path::Path;

/// LAMMPS's default skin in reduced units (neighbors.rs). Our own choice.
pub const SKIN: f64 = 0.3;

/// A scenario in progress.
pub struct Run<'a> {
    scenario: &'a Scenario,
    threads: usize,
    record: Record,
    /// The atoms now.
    pub system: System,
    /// The run at its last checkpoint.
    pub state: Checkpoint,
}

impl<'a> Run<'a> {
    /// Starts `scenario` from its beginning, or with `resume` from the
    /// checkpoint in `folder`, where the files of the run are kept. Forces
    /// use `threads` threads.
    pub fn new(s: &'a Scenario, folder: &Path, threads: usize, resume: bool) -> io::Result<Self> {
        let (r, v) = start::atoms(s);
        let (system, state) = if resume {
            let state = record::load(folder)?;
            if state.positions.len() != r.len() || state.cutoff != s.cutoff {
                return Err(io::Error::other("the checkpoint is not from this scenario"));
            }
            (state.system(threads), state)
        } else {
            let (space, method) = (PeriodicBox { side: s.side }, Method::listed(SKIN, threads));
            let system = System::with_forces(r, v, space, s.cutoff, method);
            let state = Checkpoint::of(&system, 0, None, None);
            (system, state)
        };
        let record = Record::open(folder, state.step)?;
        Ok(Self {
            scenario: s,
            threads,
            record,
            system,
            state,
        })
    }

    /// Runs stage `index` of the scenario, or what is left of it, up to
    /// step `stop` of the whole run if one is given. A checkpoint is saved
    /// every `checkpoint` steps of the run and at the end.
    pub fn stage(&mut self, index: usize, stop: Option<u64>) -> io::Result<()> {
        let (s, system) = (self.scenario, &mut self.system);
        let stage = &s.stages[index];
        let first: u64 = s.stages[..index].iter().map(|st| st.steps as u64).sum();
        let last = first + stage.steps as u64;
        let end = stop.map_or(last, |stop| stop.min(last));
        let mut step = self.state.step;
        if step < first || step >= end {
            return Ok(());
        }
        // A stage picked up part of the way carries on with the saved
        // thermostat and barostat.
        let fresh = step == first;
        let mut thermostat = stage
            .thermostat
            .map(|(t, period)| Thermostat::new(t, period));
        match &mut thermostat {
            Some(t) if fresh => thermo::rescale(&mut system.velocities, t.temperature),
            Some(t) => (t.friction, t.eta) = self.state.thermostat,
            None => {}
        }
        let (tries, kept) = if fresh { (0, 0) } else { self.state.barostat };
        let seed = s.seed.wrapping_add(1 + index as u64);
        let mut barostat = stage.barostat.map(|(p, max_step, every)| {
            let (t, _) = stage.thermostat.expect("a barostat needs a thermostat");
            (Barostat::new(p, t, max_step, seed, tries, kept), every)
        });
        while step < end {
            if let Some(t) = &mut thermostat {
                t.half_step(&mut system.velocities, s.dt);
            }
            system.step(s.dt);
            if let Some(t) = &mut thermostat {
                t.half_step(&mut system.velocities, s.dt);
            }
            step += 1;
            let due = |every: usize| every > 0 && (step - first).is_multiple_of(every as u64);
            match &mut barostat {
                Some((b, every)) if due(*every) => {
                    b.attempt(system, Method::listed(SKIN, self.threads));
                }
                _ => {}
            }
            if due(stage.sample) {
                self.record.sample(step, s.dt, system)?;
            }
            if due(stage.profile) {
                self.record.slices(step, system)?;
            }
            if step == end || (s.checkpoint > 0 && step.is_multiple_of(s.checkpoint as u64)) {
                let b = barostat.as_ref().map(|(b, _)| b);
                self.state = Checkpoint::of(system, step, thermostat.as_ref(), b);
                self.record.save(&self.state)?;
            }
        }
        Ok(())
    }
}
