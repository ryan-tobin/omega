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

Rahman's liquid argon, compared with his paper. cargo run --release --bin
argon melts 864 atoms in his box at his temperature, measures 200 tau of
free motion in 10 blocks, and prints our numbers next to his: the diffusion
coefficient two ways, from the mean squared displacement and from the
velocity autocorrelation function, with error bars from the spread of the
blocks, the peaks of g(r) and of the structure factor S(k), and when the
velocity autocorrelation first turns negative. It writes the curves to
output/argon and four plots to docs. The write-up, with the plots, is
[docs/liquid-argon.md](docs/liquid-argon.md). Five tests in tests/rahman.rs
check the same numbers against his. They take about 90 s with --release,
so a plain cargo test skips them.

Scale. The forces can also come from a neighbor list: each atom's
neighbors within the cutoff plus a skin of 0.3 sigma, found by cutting the
box into cells and searching only the cells around each atom. The list is
built again once some atom has moved more than half the skin, and the
atoms are split across threads. The cost of a step then grows in
proportion to the number of atoms instead of as its square. Tests in
tests/forces.rs and tests/scale.rs check the list against all pairs, at
one moment and along runs, and check that any number of threads gives the
same result to the last bit. tests/listed_liquid.rs runs Rahman's liquid
on the list and checks its temperature and diffusion coefficient against
the all-pairs run; like the Rahman tests it is skipped by a plain cargo
test. cargo run --release --bin scale times a step from 864 to 1,000,188
atoms: all pairs, the list on 1 thread, and the list on every thread the
machine offers. It takes about 4 minutes on 4 threads.

## Run it

    cargo run --release
    cargo run --release --bin measure
    cargo run --release --bin argon
    cargo test
    cargo test --release -- --ignored

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