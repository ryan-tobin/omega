# omega

Layer 1 of Project Omega: atoms. Simulation code written from the ground up
in Rust, one checked step at a time, with no third-party crates.

## What it computes today

The size in SI of each reduced unit for argon, using sigma = 3.4 angstrom and
eps/kB = 120 K (Rahman, Phys. Rev. 136, A405 (1964)) and a mass of 39.95 u
(CIAAW abridged standard atomic weight). Six tests in src/units.rs check the
results against values worked out by hand from the CODATA 2022 constants.

That is all so far. There is no simulation yet.

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