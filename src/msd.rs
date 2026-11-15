//! How far atoms wander: the mean squared displacement (MSD), and the
//! diffusion coefficient D that comes from it.
//!
//! MSD(t) is the average over atoms of |r_i(t) - r_i(0)|^2. Positions are
//! never wrapped back into the box (system.rs), so an atom that leaves
//! through one face keeps counting the whole distance it went.
//!
//! At first an atom flies freely and the MSD grows as t^2. In a liquid,
//! after a few collisions, each atom wanders like a random walk and the
//! MSD grows in a straight line, MSD = 6 D t: 2 D t for each of the three
//! directions (Einstein, 1905). D is one sixth of the slope.
//!
//! Label: approximated. D is measured from one starting moment over a
//! finite stretch of time, for atoms in a finite repeating box.

use crate::vec3::Vec3;

/// The starting positions r_i(0) that displacements are measured from.
pub struct Msd {
    origin: Vec<Vec3>,
}

impl Msd {
    /// Takes `positions`, in sigma, as the starting positions.
    pub fn new(positions: &[Vec3]) -> Self {
        Self {
            origin: positions.to_vec(),
        }
    }

    /// The mean squared displacement of `positions` from the starting
    /// positions, in sigma^2. Positions are in sigma and must not have
    /// been wrapped back into the box.
    pub fn measure(&self, positions: &[Vec3]) -> f64 {
        assert_eq!(positions.len(), self.origin.len(), "atom count changed");
        let mut sum = 0.0;
        for (r, &r0) in positions.iter().zip(&self.origin) {
            sum += (*r - r0).norm_sq();
        }
        sum / positions.len() as f64
    }
}

/// Diffusion coefficient D, in sigma^2 / tau: one sixth of the slope of
/// the least-squares straight line through the points (t, MSD), with
/// `times` in tau and `msd` in sigma^2.
pub fn diffusion(times: &[f64], msd: &[f64]) -> f64 {
    assert_eq!(times.len(), msd.len());
    assert!(times.len() >= 2, "a slope needs two or more points");
    let n = times.len() as f64;
    let t_mean = times.iter().sum::<f64>() / n;
    let m_mean = msd.iter().sum::<f64>() / n;
    let mut across = 0.0;
    let mut spread = 0.0;
    for (&t, &m) in times.iter().zip(msd) {
        across += (t - t_mean) * (m - m_mean);
        spread += (t - t_mean) * (t - t_mean);
    }
    across / spread / 6.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Rng;

    #[test]
    fn atoms_flying_freely_give_v_squared_t_squared() {
        // Two atoms at constant velocity for t = 2. |v|^2 is
        // 1 + 4 + 0.25 = 5.25 and 0.25 + 0 + 1 = 1.25, averaging 3.25, so
        // MSD = 3.25 * 4 = 13. Every value is exact in an f64.
        let start = [Vec3::new(0.5, 7.5, 3.0), Vec3::new(2.0, 0.5, 7.5)];
        let v = [Vec3::new(1.0, 2.0, -0.5), Vec3::new(-0.5, 0.0, 1.0)];
        let msd = Msd::new(&start);
        let later = [start[0] + v[0] * 2.0, start[1] + v[1] * 2.0];
        assert_eq!(msd.measure(&start), 0.0);
        assert_eq!(msd.measure(&later), 13.0);
    }

    #[test]
    fn straight_line_gives_its_slope_over_6() {
        // MSD = 6 * 0.25 * t + 0.3 at t = 0, 0.5, ..., 10. Tolerance 1e-12,
        // for rounding, our own choice.
        let times: Vec<f64> = (0..=20).map(|i| f64::from(i) * 0.5).collect();
        let msd: Vec<f64> = times.iter().map(|t| 1.5 * t + 0.3).collect();
        let d = diffusion(&times, &msd);
        assert!((d - 0.25).abs() <= 1e-12, "D {d}");
    }

    #[test]
    fn random_walk_gives_its_known_diffusion_coefficient() {
        // 5000 walkers, each taking a step of +H or -H along every axis at
        // every time step DT, chosen at random. One step adds 3 H^2 to the
        // squared displacement on average, so MSD = 3 H^2 t / DT and
        // D = H^2 / (2 DT) = 0.5. Over 40 seeds the measured D had a
        // standard deviation of 1.2% and was never more than 2.7% off.
        // This seed gives 1.0% low. Tolerance: 5%, our own choice.
        const H: f64 = 0.1;
        const DT: f64 = 0.01;
        let mut rng = Rng::new(31);
        let mut positions = vec![Vec3::new(0.0, 0.0, 0.0); 5000];
        let msd = Msd::new(&positions);
        let (mut times, mut values) = (Vec::new(), Vec::new());
        let mut step = || if rng.uniform() < 0.5 { H } else { -H };
        for n in 1..=200 {
            for r in &mut positions {
                *r = *r + Vec3::new(step(), step(), step());
            }
            times.push(f64::from(n) * DT);
            values.push(msd.measure(&positions));
        }
        let d = diffusion(&times, &values);
        assert!((d / 0.5 - 1.0).abs() <= 0.05, "D {d}");
    }
}
