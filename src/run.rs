//! A measured free run, cut into blocks: the temperature, g(r), the mean
//! squared displacement (MSD), the velocity autocorrelation function
//! (VACF), and the diffusion coefficient D measured two ways.
//!
//! In each block of BLOCK_STEPS steps the MSD and the VACF are averaged
//! over every time origin (correlate.rs). D comes from each: one sixth of
//! the least-squares slope of the MSD from lag FIT_FROM on (msd.rs), and
//! one third of the integral of the VACF. The spread of the block values
//! gives the error bar of their mean (stats.rs). g(r) spans the whole run.
//!
//! Label: approximated. The blocks are 20 tau long and the VACF dies away
//! within about 2 tau, so the blocks are close to independent.

use crate::correlate::Correlation;
use crate::msd;
use crate::rdf::Rdf;
use crate::stats;
use crate::system::System;
use crate::thermo;

// Our own choices. A block is 4000 steps, sampled every 5. With
// dt = 0.005 tau that is 20 tau, sampled every 0.025 tau. The MSD runs to
// lag 400 (10 tau) and is fitted from lag 80 (2 tau) on, after each atom
// has forgotten its first speed. The VACF runs to lag 200 (5 tau), where
// it has died away. g(r) is gathered every 50 steps out to 5.1 sigma,
// just under half of Rahman's box, in 255 bins of 0.02 sigma.
pub const BLOCK_STEPS: usize = 4000;
pub const SAMPLE_EVERY: usize = 5;
pub const MSD_LAGS: usize = 400;
pub const FIT_FROM: usize = 80;
pub const VACF_LAGS: usize = 200;
pub const RDF_EVERY: usize = 50;
pub const RDF_MAX: f64 = 5.1;
pub const RDF_BINS: usize = 255;

/// What one block measured: the average temperature, in eps / kB, D from
/// the MSD and from the VACF, in sigma^2 / tau, and the MSD and VACF for
/// each lag, in sigma^2 and (sigma / tau)^2.
pub struct Block {
    pub temperature: f64,
    pub d_msd: f64,
    pub d_vacf: f64,
    pub msd: Vec<f64>,
    pub vacf: Vec<f64>,
}

/// The whole run: its blocks, g(r) over all of them, the total energy at
/// the start, in eps, and the largest relative change from it so far.
pub struct Run {
    pub blocks: Vec<Block>,
    pub rdf: Rdf,
    pub start_energy: f64,
    pub energy_change: f64,
}

impl Run {
    /// A run with no blocks yet, starting from `system` as it is now.
    pub fn new(system: &System) -> Self {
        Self {
            blocks: Vec::new(),
            rdf: Rdf::new(system.space, RDF_MAX, RDF_BINS),
            start_energy: system.total_energy(),
            energy_change: 0.0,
        }
    }

    /// Moves `system` freely for one more block, with time step `dt` in
    /// tau, measures it, and returns what the block measured. The time
    /// between samples is SAMPLE_EVERY * dt.
    pub fn add_block(&mut self, system: &mut System, dt: f64) -> &Block {
        let mut msd = Correlation::msd(MSD_LAGS);
        let mut vacf = Correlation::vacf(VACF_LAGS);
        let mut t_sum = 0.0;
        for step in 0..BLOCK_STEPS {
            if step % SAMPLE_EVERY == 0 {
                msd.add(&system.positions);
                vacf.add(&system.velocities);
                t_sum += thermo::temperature(&system.velocities);
            }
            if step % RDF_EVERY == 0 {
                self.rdf.add(&system.positions);
            }
            system.step(dt);
            let change = (system.total_energy() - self.start_energy).abs();
            let change = change / self.start_energy.abs();
            self.energy_change = self.energy_change.max(change);
        }
        let h = SAMPLE_EVERY as f64 * dt;
        let (msd, vacf) = (msd.values(), vacf.values());
        let times: Vec<f64> = (FIT_FROM..=MSD_LAGS).map(|k| k as f64 * h).collect();
        self.blocks.push(Block {
            temperature: t_sum / (BLOCK_STEPS / SAMPLE_EVERY) as f64,
            d_msd: msd::diffusion(&times, &msd[FIT_FROM..]),
            d_vacf: stats::trapezoid(&vacf, h) / 3.0,
            msd,
            vacf,
        });
        &self.blocks[self.blocks.len() - 1]
    }

    /// The mean over blocks of the number `f` picks, and its standard
    /// error, in the unit of that number.
    pub fn mean(&self, f: fn(&Block) -> f64) -> (f64, f64) {
        let values: Vec<f64> = self.blocks.iter().map(f).collect();
        stats::mean_and_error(&values)
    }

    /// The MSD or the VACF, picked by `f`, averaged over blocks lag by lag.
    pub fn curve(&self, f: fn(&Block) -> &Vec<f64>) -> Vec<f64> {
        let n = self.blocks.len() as f64;
        let mean = |k: usize| self.blocks.iter().map(|b| f(b)[k]).sum::<f64>() / n;
        (0..f(&self.blocks[0]).len()).map(mean).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::periodic::PeriodicBox;
    use crate::vec3::Vec3;

    #[test]
    fn atoms_flying_freely_give_known_values_in_every_block() {
        // Two atoms 50 sigma apart in a box of side 100, at the same
        // velocity v, |v|^2 = 0.25, so they never feel each other. T is
        // 2 * 0.25 / (3 * 2 - 3) = 1/6. The VACF is |v|^2 at every lag, so
        // D from the VACF is 0.25 * 5 / 3. The MSD is |v|^2 t^2, and the
        // least-squares slope of t^2 over evenly spaced t from 2 to 10 is
        // 2 + 10 = 12, so D from the MSD is 0.25 * 12 / 6 = 0.5. Measured:
        // within 4.3e-13, from rounding in the 8000 position updates.
        // Tolerance 1e-10 relative, our own choice.
        let v = Vec3::new(0.3, 0.4, 0.0);
        let positions = vec![Vec3::new(1.0, 1.0, 1.0), Vec3::new(51.0, 1.0, 1.0)];
        let space = PeriodicBox { side: 100.0 };
        let mut system = System::new(positions, vec![v, v], space, 2.5);
        let mut run = Run::new(&system);
        run.add_block(&mut system, 0.005);
        run.add_block(&mut system, 0.005);
        let close = |got: f64, want: f64| (got / want - 1.0).abs() <= 1e-10;
        for b in &run.blocks {
            assert!(close(b.temperature, 0.25 / 1.5), "T {}", b.temperature);
            assert!(close(b.d_vacf, 0.25 * 5.0 / 3.0), "D vacf {}", b.d_vacf);
            assert!(close(b.d_msd, 0.5), "D msd {}", b.d_msd);
        }
        assert!(run.curve(|b| &b.vacf).iter().all(|&v| close(v, 0.25)));
        let (d, error) = run.mean(|b| b.d_msd);
        assert!(close(d, 0.5) && error <= 1e-10, "D {d} +- {error}");
    }
}