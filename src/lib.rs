//! Omega, layer 1: atoms.
//!
//! What this crate computes today: the size in SI of each reduced unit
//! for argon, three-component vector math, the Lennard-Jones energy and
//! force between two atoms, and the motion of many atoms in a repeating
//! box, stepped forward in time by velocity Verlet.

pub mod lattice;
pub mod lj;
pub mod ln;
pub mod pairs;
pub mod periodic;
pub mod rng;
pub mod system;
pub mod units;
pub mod vec3;
pub mod velocities;
