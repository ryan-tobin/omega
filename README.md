# omega

Layer 1 of Project Omega: atoms. Simulation code written from the ground up
in Rust, one checked step at a time, with no third-party crates.

## What it computes today

The size in SI of each reduced unit for argon, using sigma = 3.4 angstrom and
eps/kB = 120 K (Rahman, Phys. Rev. 136, A405 (1964)) and a mass of 39.95 u
(CIAAW abridged standard atomic weight). Six tests in src/units.rs check the
results against values worked out by hand from the CODATA 2022 constants.

Three-component vector math: sums, differences, negation, scaling, dot
products and lengths. Eight tests in src/vec3.rs check cases worked out by
hand, using whole numbers and halves so that the answers are exact.

The Lennard-Jones energy and force between two atoms, in reduced units. Seven
tests in src/lj.rs check that the energy is zero at one sigma, that the minimum
sits at 2^(1/6) sigma where the energy is minus eps, and that the force equals
minus the measured slope of the energy to within 1e-6 eps/sigma.

Many atoms in a box that repeats in every direction. The atoms start on a
face-centered cubic lattice, with velocities drawn for a chosen temperature
from our own random number generator. Each atom feels the Lennard-Jones force
of every other atom within the cutoff, and velocity Verlet steps them forward
in time. Tests in tests/conservation.rs check, on 108 atoms, that total
energy holds and total momentum stays zero over 10,000 steps, and that 1000
steps forward and then 1000 back return the atoms to their start.
cargo run --release prints the same quantities for 864 atoms in a box the
size of Rahman's, with starting velocities drawn for his temperature.

Measurements, each checked by tests on a case with a known answer:
temperature from the kinetic energy, pressure from the virial, the radial
distribution function g(r), the mean squared displacement, and the diffusion
coefficient from its slope. Velocity rescaling brings the atoms to a chosen
temperature before anything is measured. Measurements are written as CSV and
positions as extended XYZ, in angstroms, which tools such as ASE can read.
cargo run --release --bin measure melts 864 atoms in Rahman's box, holds
them at his temperature, then measures them for 10,000 free steps and writes
output/thermo.csv, output/msd.csv, output/rdf.csv and output/positions.xyz.
Comparing these with Rahman's published numbers comes next.

## Run it

    cargo run --release
    cargo run --release --bin measure
    cargo test

Use --release: the 864-atom runs are slow without it.

Built and tested with Rust 1.100.0 on macOS (Apple Silicon).

## License

Licensed under either of

 * Apache License, Version 2.0
   ([LICENSE-APACHE](LICENSE-APACHE) or <http://www.apache.org/licenses/LICENSE-2.0>)
 * MIT license
   ([LICENSE-MIT](LICENSE-MIT) or <http://opensource.org/licenses/MIT>)

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.