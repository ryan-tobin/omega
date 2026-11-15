//! Omega, layer 1: atoms.
//!
//! What this crate computes today: the size in SI of each reduced unit
//! for argon, three-component vector math, the Lennard-Jones energy and
//! force between two atoms, and the motion of many atoms in a repeating
//! box, stepped forward in time by velocity Verlet. It measures their
//! temperature, potential energy, pressure, radial distribution function,
//! mean squared displacement and diffusion coefficient, and writes them
//! as CSV and the positions as XYZ.

pub mod csv;
pub mod lattice;
pub mod lj;
pub mod ln;
pub mod msd;
pub mod pairs;
pub mod periodic;
pub mod rdf;
pub mod rng;
pub mod system;
pub mod thermo;
pub mod units;
pub mod vec3;
pub mod velocities;
pub mod xyz;
