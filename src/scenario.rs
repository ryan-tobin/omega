//! Scenario files: one text file describes a full run.
//!
//! Each line is a key and its values, separated by spaces. A `#` starts a
//! comment. The lines before the first `stage` line describe the start:
//!
//!   cells 6           fcc cells along each side: 4 * 6^3 lattice sites
//!   side 10.229       box side, in sigma
//!   fill 0.5          optional: atoms only where x is below this fraction
//!                     of the side, which leaves the rest of the box empty
//!   temperature 0.8   the starting velocities are drawn for it, eps / kB
//!   hot 0.5 8 0.9     optional: the atoms at or beyond this fraction of the
//!                     side start at this temperature instead, and only
//!                     this share of them is kept, picked at random
//!   seed 1964         for the starting velocities and the barostat
//!   cutoff 2.5        force cutoff, in sigma
//!   dt 0.005          time step, in tau
//!   checkpoint 5000   optional: save the state every so many steps
//!
//! Each `stage` line starts a stretch of the run, and the lines below it
//! set it up:
//!
//!   stage measure       its name
//!   steps 20000         how many time steps
//!   thermostat 0.8 0.5  optional: temperature and period (thermostat.rs)
//!   barostat 1.0 0.05 25  optional: pressure, largest change of the box
//!                       side, and steps between tries (barostat.rs)
//!   sample 10           optional: a row of thermo.csv every so many steps
//!   profile 100         optional: a row of profile.csv every so many steps
//!
//! A stage without a thermostat moves freely, at constant energy. A
//! barostat needs a thermostat in the same stage.

/// One stretch of a run. Fields as in the file; 0 means never.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Stage {
    pub name: String,
    pub steps: usize,
    /// Temperature in eps / kB and period in tau.
    pub thermostat: Option<(f64, f64)>,
    /// Pressure in eps / sigma^3, largest step in sigma, steps per try.
    pub barostat: Option<(f64, f64, usize)>,
    pub sample: usize,
    pub profile: usize,
}

/// A full run. Fields as in the file, in reduced units.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Scenario {
    pub cells: usize,
    pub side: f64,
    pub fill: f64,
    pub temperature: f64,
    pub hot: Option<(f64, f64, f64)>,
    pub seed: u64,
    pub cutoff: f64,
    pub dt: f64,
    pub checkpoint: usize,
    pub stages: Vec<Stage>,
}

/// Reads the text of a scenario file. An error names the line it is on.
pub fn parse(text: &str) -> Result<Scenario, String> {
    let (mut s, mut seeded) = (Scenario::default(), false);
    s.fill = 1.0;
    for (index, line) in text.lines().enumerate() {
        let line = line.split('#').next().unwrap_or("").trim();
        let mut words = line.split_whitespace();
        let Some(key) = words.next() else { continue };
        let bad = || format!("line {}: cannot use `{line}`", index + 1);
        if key == "stage" {
            let name = words.next().ok_or_else(bad)?.to_string();
            s.stages.push(Stage {
                name,
                ..Stage::default()
            });
            continue;
        }
        let values: Result<Vec<f64>, _> = words.map(str::parse).collect();
        let values = values.map_err(|_| bad())?;
        // A count: a whole number, zero or more.
        let whole = |x: f64| {
            (x >= 0.0 && x.fract() == 0.0)
                .then_some(x as usize)
                .ok_or_else(bad)
        };
        match (key, s.stages.last_mut(), values.as_slice()) {
            ("cells", None, [n]) => s.cells = whole(*n)?,
            ("side", None, [x]) => s.side = *x,
            ("fill", None, [x]) => s.fill = *x,
            ("temperature", None, [t]) => s.temperature = *t,
            ("hot", None, [from, t, keep]) => s.hot = Some((*from, *t, *keep)),
            ("seed", None, [n]) => (s.seed, seeded) = (whole(*n)? as u64, true),
            ("cutoff", None, [x]) => s.cutoff = *x,
            ("dt", None, [x]) => s.dt = *x,
            ("checkpoint", None, [n]) => s.checkpoint = whole(*n)?,
            ("steps", Some(stage), [n]) => stage.steps = whole(*n)?,
            ("thermostat", Some(stage), [t, period]) => stage.thermostat = Some((*t, *period)),
            ("barostat", Some(stage), [p, step, every]) => {
                stage.barostat = Some((*p, *step, whole(*every)?));
            }
            ("sample", Some(stage), [n]) => stage.sample = whole(*n)?,
            ("profile", Some(stage), [n]) => stage.profile = whole(*n)?,
            _ => return Err(bad()),
        }
    }
    check(&s, seeded)?;
    Ok(s)
}

/// Refuses a scenario no run can use.
fn check(s: &Scenario, seeded: bool) -> Result<(), String> {
    let set = s.cells > 0 && s.side > 0.0 && s.temperature > 0.0 && s.cutoff > 0.0 && s.dt > 0.0;
    if !set || !seeded {
        return Err("cells, side, temperature, seed, cutoff and dt must all be set".into());
    }
    let share = |x: f64| x > 0.0 && x <= 1.0;
    if !share(s.fill) || s.hot.is_some_and(|(_, t, keep)| t <= 0.0 || !share(keep)) {
        return Err(
            "fill and the share kept must be above 0 and at most 1, a hot temperature above 0"
                .into(),
        );
    }
    if s.stages.is_empty() {
        return Err("a run needs at least one stage".into());
    }
    for stage in &s.stages {
        let name = &stage.name;
        if stage.steps == 0 {
            return Err(format!("stage {name}: steps must be set"));
        }
        match (stage.thermostat, stage.barostat) {
            (None, Some(_)) => return Err(format!("stage {name}: a barostat needs a thermostat")),
            (Some((t, period)), _) if t <= 0.0 || period <= 0.0 => {
                return Err(format!("stage {name}: thermostat values must be above 0"));
            }
            (_, Some((_, step, every))) if step <= 0.0 || every == 0 => {
                return Err(format!(
                    "stage {name}: barostat step and interval must be above 0"
                ));
            }
            _ => {}
        }
    }
    Ok(())
}
