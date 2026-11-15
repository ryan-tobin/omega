//! The radial distribution function g(r): how many atoms sit a distance r
//! from an atom, compared with an even spread at the same density.
//!
//! Pair distances up to r_max are counted in bins of width dr. Bin k runs
//! from r_k = k dr to r_k + dr and holds the shell of volume
//! V_k = (4 pi / 3) ((r_k + dr)^3 - r_k^3). With H_k unique pairs counted
//! in bin k over S snapshots of N atoms at density rho = N / V,
//! g_k = 2 H_k / (S N rho V_k). Each pair is seen from both of its atoms,
//! hence the 2. An even spread puts rho V_k atoms in the shell, so g = 1
//! means no more and no fewer than that. This is the normalization of the
//! pair_distribution example of Allen and Tildesley, Computer Simulation
//! of Liquids, 2nd ed. (2017). For a thin shell V_k is 4 pi r^2 dr, the
//! form in PHYSICS_REFERENCE.
//!
//! r_max must be at most L/2, so that each pair is counted once, at its
//! nearest copy (periodic.rs).
//!
//! Label: approximated. g is averaged over the width of each bin. Each
//! atom has N - 1 others around it, not N, so atoms placed at random give
//! g = 1 - 1/N.

use crate::periodic::PeriodicBox;
use crate::vec3::Vec3;
use std::f64::consts::PI;

/// A histogram of pair distances, built up one snapshot at a time.
pub struct Rdf {
    space: PeriodicBox,
    dr: f64,
    counts: Vec<u64>,
    snapshots: u64,
    atoms: usize,
}

impl Rdf {
    /// An empty histogram of `bins` bins from 0 to `r_max` , in sigma, for
    /// atoms in `space`. `r_max` must be at most half the box side.
    pub fn new(space: PeriodicBox, r_max: f64, bins: usize) -> Self {
        assert!(r_max <= 0.5 * space.side, "r_max is over half the box");
        Self {
            space,
            dr: r_max / bins as f64,
            counts: vec![0; bins],
            snapshots: 0,
            atoms: 0,
        }
    }

    /// Counts every pair distance in one snapshot of `positions`, in
    /// sigma. Every snapshot must have the same number fo atoms.
    pub fn add(&mut self, positions: &[Vec3]) {
        if self.snapshots == 0 {
            self.atoms = positions.len();
        }
        assert_eq!(positions.len(), self.atoms, "atom count changed");
        let bins = self.counts.len();
        for i in 0..positions.len() {
            for j in (i + 1)..positions.len() {
                let r = self.space.nearest_image(positions[i] - positions[j]);
                let k = (r.norm() / self.dr) as usize;
                if k < bins {
                    self.counts[k] += 1;
                }
            }
        }
        self.snapshots += 1;
    }

    /// H_k: unique pairs counted in each bin, over all snapshots. Unit: none.
    pub fn counts(&self) -> &[u64] {
        &self.counts
    }

    /// One (r, g) for each bin: r at the middle of the bin, in sigma, and
    /// g, a plain number.
    pub fn g(&self) -> Vec<(f64, f64)> {
        let l = self.space.side;
        let n = self.atoms as f64;
        let density = n / (l * l * l);
        let mut out = Vec::with_capacity(self.counts.len());
        for (k, &h) in self.counts.iter().enumerate() {
            let lower = k as f64 * self.dr;
            let upper = lower + self.dr;
            let shell = 4.0 / 3.0 * PI * (upper * upper * upper - lower * lower * lower);
            let g = 2.0 * h as f64 / (self.snapshots as f64 * n * density * shell);
            out.push((lower + 0.5 * self.dr, g));
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lattice;
    use crate::rng::Rng;

    #[test]
    fn perfect_lattice_puts_each_shell_in_one_bin() {
        // 108 atoms on an fcc lattice of cube side a = 1 in a box of side
        // 3, so rho = 4. Up to L/2 = 1.5 each atom has 12, 6, 24 and 12
        // neighbors at sqrt(1/2), 1, sqrt(3/2) and sqrt(2) (lattice.rs).
        // With 40 bins of width 0.0375 these fall in bins 18, 26, 32 and
        // 37, none near a bin edge. Each shell holds 108 * n / 2 unique
        // pairs. Two snapshots double every count but not g.
        let space = PeriodicBox { side: 3.0 };
        let positions = lattice::fcc(3, space.side);
        let mut rdf = Rdf::new(space, 1.5, 40);
        rdf.add(&positions);
        rdf.add(&positions);
        let mut want = [0; 40];
        for (k, n) in [(18, 12), (26, 6), (32, 24), (37, 12)] {
            want[k] = 2 * 108 * n / 2;
        }
        assert_eq!(rdf.counts(), want);
        // First shell: g = 12 / (rho V_18), V_18 from 0.675 to 0.7125.
        // Tolerance 1e-12 relative, for rounding, our own choice.
        let g = rdf.g();
        let shell = 4.0 / 3.0 * PI * (0.7125_f64.powi(3) - 0.675_f64.powi(3));
        let want_g = 12.0 / (4.0 * shell);
        assert!((g[18].1 / want_g - 1.0).abs() <= 1e-12, "g {}", g[18].1);
        assert!((g[18].0 - 0.69375).abs() <= 1e-12, "r {}", g[18].0);
    }

    #[test]
    fn atoms_placed_at_random_give_g_of_1_minus_1_over_n() {
        // 100 atoms placed evenly at random in a box of side 5, 200 times.
        // Each bin from 0.5 to 2.5 expects 3,000 to 60,000 counts, so its
        // relative standard error is 1.8% at most. Tolerance: 5 standard
        // errors in each bin, our own choice.
        let space = PeriodicBox { side: 5.0 };
        let mut rng = Rng::new(21);
        let mut rdf = Rdf::new(space, 2.5, 25);
        for _ in 0..200 {
            let mut positions = Vec::with_capacity(100);
            for _ in 0..100 {
                let x = rng.uniform();
                let y = rng.uniform();
                let z = rng.uniform();
                positions.push(Vec3::new(x, y, z) * space.side);
            }
            rdf.add(&positions);
        }
        let want = 1.0 - 1.0 / 100.0;
        for (k, (r, g)) in rdf.g().into_iter().enumerate().skip(5) {
            let error = 1.0 / (rdf.counts()[k] as f64).sqrt();
            assert!((g / want - 1.0).abs() <= 5.0 * error, "g({r}) = {g}");
        }
    }
}
