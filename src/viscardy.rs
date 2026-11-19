//! Published viscosity and thermal conductivity of the Lennard-Jones
//! liquid at T = 0.722, density 0.8442, cutoff 2.5 sigma, by the number of
//! atoms in the box.
//!
//! S. Viscardy, J. Servantie and P. Gaspard, "Transport and Helfand
//! moments in the Lennard-Jones fluid. I. Shear viscosity" (2007),
//! arXiv:cond-mat/0701253, Table II, and "II. Thermal conductivity",
//! arXiv:cond-mat/0701254, Table I, the rows "This work 2007 HM (MD)".
//! Read in the arXiv texts on 2026-11-18. Their runs: velocity Verlet with
//! dt = 0.003, motion that conserves energy (no thermostat is mentioned),
//! 2e8 to 6e8 steps per box. Their numbers come from Helfand moments, the
//! integrals of our currents, fitted from 5 to 10 tau (viscosity) and 2 to
//! 8 tau (conductivity); they say Green-Kubo agreed (their Fig. 1). They
//! do not say whether their potential was shifted. The forces, and so the
//! motion, are the same either way.
//!
//! Their infinite-box values are fits in 1/N: viscosity 3.291, with
//! +- 0.029 in the table and +- 0.057 in the text, and conductivity 6.990,
//! +- 0.061 in the table and +- 0.030 in the text. Only the boxes are used
//! here.

/// The state: temperature in eps / kB, number density in sigma^-3 and the
/// cutoff in sigma
pub const TEMPERATURE: f64 = 0.722;
pub const DENSITY: f64 = 0.8442;
pub const CUTOFF: f64 = 2.5;

/// One box: the number of atoms, the viscosity and its error in
/// (m eps)^(1/2) / sigma^2, and the conductivity and its error in
/// kB (eps / m)^(1/2) / sigma^2.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Row {
    pub atoms: usize,
    pub viscosity: (f64, f64),
    pub conductivity: (f64, f64),
}

/// Their five boxes: 4 n^3 atoms for n = 3 to 7
pub const ROWS: [Row; 5] = [
    row(108, (3.057, 0.045), (6.663, 0.10)),
    row(256, (3.232, 0.060), (6.867, 0.12)),
    row(500, (3.270, 0.050), (6.932, 0.10)),
    row(864, (3.247, 0.033), (6.948, 0.07)),
    row(1372, (3.268, 0.055), (6.946, 0.12)),
];

const fn row(atoms: usize, viscosity: (f64, f64), conductivity: (f64, f64)) -> Row {
    Row {
        atoms,
        viscosity,
        conductivity,
    }
}

/// THeir row for a box of `atoms` atoms, if they ran one.
pub fn row_for(atoms: usize) -> Option<Row> {
    ROWS.iter().copied().find(|r| r.atoms == atoms)
}
