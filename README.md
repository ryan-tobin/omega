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

That is all so far. Nothing moves yet: there is no time stepping and no
simulation.

## Run it

    cargo run
    cargo test

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