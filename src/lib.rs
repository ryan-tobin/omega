//! Omega, layer 1: atoms.
//!
//! What this crate computes today: the size in SI of each reduced unit
//! for argon, three-component vector math, the Lennard-Jones energy and
//! force between two atoms, and the motion of many atoms in a repeating
//! box, stepped forward in time by velocity Verlet. It measures their
//! temperature, potential energy, pressure, radial distribution function,
//! mean squared displacement and diffusion coefficient, and writes them
//! as CSV and the positions as XYZ.
//!
//! It also runs liquid argon as Rahman did in 1964 and compares the result
//! with his paper: the mean squared displacement and velocity
//! autocorrelation over every time origin, error bars from blocks of the
//! run, peaks of curves, the structure factor S(k), and plots drawn as SVG.
//!
//! Since M5 the forces can also come from a neighbor list built from cells,
//! with the work split across threads, so the cost of a step grows in
//! proportion to the number of atoms instead of as its square.

pub mod cells;
pub mod correlate;
pub mod csv;
pub mod forces;
pub mod lattice;
pub mod lj;
pub mod ln;
pub mod msd;
pub mod neighbors;
pub mod pairs;
pub mod parallel;
pub mod peaks;
pub mod periodic;
pub mod plot;
pub mod rahman;
pub mod rdf;
pub mod rng;
pub mod run;
pub mod stats;
pub mod structure;
pub mod svg;
pub mod system;
pub mod thermo;
pub mod units;
pub mod vec3;
pub mod velocities;
pub mod xyz;