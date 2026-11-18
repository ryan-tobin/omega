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
//!
//! Since M6 a run can be saved to a folder and read back exactly, and a
//! saved moment drawn as a picture: the atoms as lit balls colored by
//! speed, with the moment's measurements as text over them. The view
//! program shows those pictures in a window. Only that program uses
//! crates from outside, two of them, and only it is built with them.
//!
//! Since M7, with the `gpu` feature, the forces can also be computed on a
//! GPU in 32-bit floats, through the wgpu crate (gpu.rs)

pub mod barostat;
pub mod cells;
pub mod checkpoint;
pub mod correlate;
pub mod csv;
pub mod font;
pub mod forces;
pub mod hud;
pub mod lattice;
pub mod lj;
pub mod ln;
pub mod msd;
pub mod neighbors;
pub mod pairs;
pub mod parallel;
pub mod peaks;
pub mod periodic;
pub mod picture;
pub mod plot;
pub mod rahman;
pub mod rdf;
pub mod recorder;
pub mod replay;
pub mod rng;
pub mod run;
pub mod scene;
pub mod stats;
pub mod structure;
pub mod svg;
pub mod system;
pub mod thermo;
pub mod thermostat;
pub mod units;
pub mod vec3;
pub mod velocities;
pub mod xyz;

#[cfg(feature = "gpu")]
pub mod gpu;
