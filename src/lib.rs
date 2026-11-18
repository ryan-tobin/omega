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
//!
//! Since M8 a run can be held at a temperature (thermostat.rs) and at a
//! pressure (barostat.rs), saved and picked up again to the last bit
//! (checkpoint.rs), and described from start to end by one text file
//! (scenario.rs), which drive.rs runs.
//!
//! Since M9 it measures how a liquid resists flow and passes heat: the
//! shear stress and heat current at each moment (flux.rs), their time
//! correlations integrated into the viscosity and the thermal
//! conductivity (greenkubo.rs, transport.rs), and the energy and pressure
//! the cutoff leaves out (tail.rs), checked against Viscardy et al.
//! (viscardy.rs) and the NIST benchmark (nist.rs).

pub mod barostat;
pub mod cells;
pub mod checkpoint;
pub mod correlate;
pub mod csv;
pub mod drive;
pub mod flux;
pub mod font;
pub mod forces;
pub mod greenkubo;
pub mod hud;
pub mod lattice;
pub mod lj;
pub mod ln;
pub mod msd;
pub mod neighbors;
pub mod nist;
pub mod pairs;
pub mod parallel;
pub mod peaks;
pub mod periodic;
pub mod picture;
pub mod plot;
pub mod profile;
pub mod rahman;
pub mod rdf;
pub mod record;
pub mod recorder;
pub mod replay;
pub mod rng;
pub mod run;
pub mod scenario;
pub mod scene;
pub mod start;
pub mod stats;
pub mod structure;
pub mod svg;
pub mod system;
pub mod tail;
pub mod thermo;
pub mod thermostat;
pub mod transport;
pub mod units;
pub mod vec3;
pub mod velocities;
pub mod viscardy;
pub mod xyz;

#[cfg(feature = "gpu")]
pub mod gpu;
