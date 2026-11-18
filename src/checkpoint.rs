//! Checkpoints: everything a run needs to carry on later as if it had
//! never stopped, as bytes.
//!
//! The file holds every f64 as its 8 bytes, so nothing is rounded. Layout,
//! every number little-endian: the 8 bytes `omega-cp`, the format version
//! as a u32, then nine numbers of 8 bytes each: step, atom count N, box
//! side, cutoff, friction, eta, barostat moves tried and kept, and the
//! skin (-1 without a neighbor list). Then the positions (3N f64), the
//! velocities (3N f64), and with a neighbor list the positions at its last
//! build (3N f64).
//!
//! The forces are not stored. They follow from the positions. With a
//! neighbor list they are summed in the order of the list (forces.rs), and
//! a different order changes the last bits, so the list is built again
//! from the positions it was last built at. The same list gives the same
//! forces to the last bit (tested in tests/restart.rs).

use crate::barostat::Barostat;
use crate::forces::Method;
use crate::neighbors::NeighborList;
use crate::periodic::PeriodicBox;
use crate::system::System;
use crate::thermostat::Thermostat;
use crate::vec3::Vec3;

const MAGIC: &[u8; 8] = b"omega-cp";
const VERSION: u32 = 1;
/// Bytes before the positions: 8 + 4 + 9 * 8.
const HEADER: usize = 84;

/// The state of a run between two steps.
#[derive(Debug, Clone, PartialEq)]
pub struct Checkpoint {
    /// Steps done since the run began. Unit: none.
    pub step: u64,
    /// Box side and force cutoff, in sigma.
    pub side: f64,
    pub cutoff: f64,
    /// The thermostat's friction xi, in 1 / tau, and its eta (thermostat.rs).
    pub thermostat: (f64, f64),
    /// The barostat's moves tried and kept (barostat.rs).
    pub barostat: (u64, u64),
    /// Positions in sigma and velocities in sigma / tau.
    pub positions: Vec<Vec3>,
    pub velocities: Vec<Vec3>,
    /// A neighbor list's skin, in sigma, and the positions at its last build.
    pub list: Option<(f64, Vec<Vec3>)>,
}

impl Checkpoint {
    /// The state of a run after `step` steps: `system` and, if the run has
    /// them, its thermostat and barostat.
    pub fn of(
        system: &System,
        step: u64,
        thermostat: Option<&Thermostat>,
        barostat: Option<&Barostat>,
    ) -> Self {
        let list = match &system.method {
            Method::Listed { list, .. } => Some((list.skin, list.built_at.clone())),
            _ => None,
        };
        Self {
            step,
            side: system.space.side,
            cutoff: system.cutoff,
            thermostat: thermostat.map_or((0.0, 0.0), |t| (t.friction, t.eta)),
            barostat: barostat.map_or((0, 0), |b| (b.tries, b.accepted)),
            positions: system.positions.clone(),
            velocities: system.velocities.clone(),
            list,
        }
    }

    /// The system this checkpoint was taken from, with its forces computed
    /// again. A neighbor list works on `threads` threads.
    pub fn system(&self, threads: usize) -> System {
        let space = PeriodicBox { side: self.side };
        let method = match &self.list {
            None => Method::AllPairs,
            Some((skin, built_at)) => {
                let mut list = NeighborList::new(*skin);
                list.build(built_at, space, self.cutoff, threads);
                Method::Listed { list, threads }
            }
        };
        let (r, v) = (self.positions.clone(), self.velocities.clone());
        System::with_forces(r, v, space, self.cutoff, method)
    }

    /// The checkpoint as bytes, in the layout at the top of this file.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = MAGIC.to_vec();
        out.extend(VERSION.to_le_bytes());
        out.extend(self.step.to_le_bytes());
        out.extend((self.positions.len() as u64).to_le_bytes());
        for x in [self.side, self.cutoff, self.thermostat.0, self.thermostat.1] {
            out.extend(x.to_le_bytes());
        }
        out.extend(self.barostat.0.to_le_bytes());
        out.extend(self.barostat.1.to_le_bytes());
        let skin = self.list.as_ref().map_or(-1.0, |(skin, _)| *skin);
        out.extend(skin.to_le_bytes());
        let mut arrays = vec![&self.positions, &self.velocities];
        arrays.extend(self.list.iter().map(|(_, at)| at));
        for x in arrays.into_iter().flatten().flat_map(|v| [v.x, v.y, v.z]) {
            out.extend(x.to_le_bytes());
        }
        out
    }

    /// Reads a checkpoint back from `bytes`. Refuses bytes that are not a
    /// checkpoint of this version, or are cut short or too long.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, String> {
        let version = VERSION.to_le_bytes();
        if bytes.len() < HEADER || bytes[..8] != MAGIC[..] || bytes[8..12] != version {
            return Err(format!("not an omega checkpoint of version {VERSION}"));
        }
        // The i-th 8-byte number after the version.
        let eight = |i: usize| -> [u8; 8] { bytes[12 + 8 * i..20 + 8 * i].try_into().unwrap() };
        let int = |i: usize| u64::from_le_bytes(eight(i));
        let float = |i: usize| f64::from_le_bytes(eight(i));
        let (n, skin) = (int(1) as usize, float(8));
        let arrays = if skin < 0.0 { 2 } else { 3 };
        let length = n
            .checked_mul(24 * arrays)
            .and_then(|m| m.checked_add(HEADER));
        if length != Some(bytes.len()) {
            return Err("the checkpoint has the wrong length".into());
        }
        // The k-th array of N vectors after the nine numbers.
        let vectors = |k: usize| -> Vec<Vec3> {
            let x = |i: usize| float(9 + 3 * k * n + i);
            let vector = |i: usize| Vec3::new(x(3 * i), x(3 * i + 1), x(3 * i + 2));
            (0..n).map(vector).collect()
        };
        Ok(Self {
            step: int(0),
            side: float(2),
            cutoff: float(3),
            thermostat: (float(4), float(5)),
            barostat: (int(6), int(7)),
            positions: vectors(0),
            velocities: vectors(1),
            list: (skin >= 0.0).then(|| (skin, vectors(2))),
        })
    }
}
