//! Time correlations averaged over many starting moments: the mean
//! squared displacement (MSD) and the velocity autocorrelation function
//! (VACF).
//!
//! M3 measured the MSD from one starting moment (msd.rs). In equilibrium
//! every moment is as good a start as any other, so here every sample is
//! a start, a time origin, and each time lag is averaged over all of them.
//!
//! A correlation keeps the last `lags + 1` frames: the positions or the
//! velocities of every atom, sampled a fixed time h apart. Each new frame
//! a(t) is paired with each kept frame a(t - k h), k = 0 to `lags`, and
//! f(a_i(t), a_i(t - k h)), averaged over the atoms i, is added to lag k.
//! For the MSD f(now, then) = |now - then|^2, and for the VACF
//! f(now, then) = now . then.
//!
//! The VACF says how long an atom keeps moving the way it was moving. Its
//! integral gives the diffusion coefficient a second way, independent of
//! the MSD fit: D = (1/3) * integral from 0 to infinity of VACF(t) dt
//! (Green and Kubo). The MSD is a double integral of the VACF, so at long
//! times its slope over 6 tends to the same D.
//!
//! Label: approximated. Each lag is an average over a finite run, so it
//! carries a statistical error.

use crate::vec3::Vec3;

/// A time correlation, built up one frame at a time.
pub struct Correlation {
    pair: fn(Vec3, Vec3) -> f64,
    frames: Vec<Vec<Vec3>>,
    added: usize,
    sums: Vec<f64>,
    counts: Vec<u64>,
}

impl Correlation {
    /// The MSD for lags 0 to `lags`: lag k is the average of
    /// |r_i(t + k h) - r_i(t)|^2, in sigma^2 for positions in sigma.
    /// Positions must not have been wrapped back into the box.
    pub fn msd(lags: usize) -> Self {
        Self::new(lags, |now, then| (now - then).norm_sq())
    }

    /// The VACF for lags 0 to `lags`: lag k is the average of
    /// v_i(t + k h) . v_i(t), in (sigma / tau)^2 for velocities in
    /// sigma / tau. At lag 0 it is the average v^2, which is 3 T.
    pub fn vacf(lags: usize) -> Self {
        Self::new(lags, |now, then| now.dot(then))
    }

    fn new(lags: usize, pair: fn(Vec3, Vec3) -> f64) -> Self {
        Self {
            pair,
            frames: Vec::with_capacity(lags + 1),
            added: 0,
            sums: vec![0.0; lags + 1],
            counts: vec![0; lags + 1],
        }
    }

    /// Adds the next frame, one vector per atom, taken h after the last.
    /// Every frame must have the same number of atoms.
    pub fn add(&mut self, frame: &[Vec3]) {
        if let Some(first) = self.frames.first() {
            assert_eq!(frame.len(), first.len(), "atom count changed");
        }
        let keep = self.sums.len();
        let slot = self.added % keep;
        if self.frames.len() < keep {
            self.frames.push(frame.to_vec());
        } else {
            self.frames[slot].copy_from_slice(frame);
        }
        let newest = &self.frames[slot];
        for k in 0..keep.min(self.added + 1) {
            let older = &self.frames[(self.added - k) % keep];
            let mut sum = 0.0;
            for (&now, &then) in newest.iter().zip(older) {
                sum += (self.pair)(now, then);
            }
            self.sums[k] += sum / frame.len() as f64;
            self.counts[k] += 1;
        }
        self.added += 1;
    }

    /// The average for each lag k = 0 to `lags`, at time k h. Needs more
    /// than `lags` frames, so that every lag has been seen.
    pub fn values(&self) -> Vec<f64> {
        assert!(self.added > self.sums.len() - 1, "fewer frames than lags");
        let pairs = self.sums.iter().zip(&self.counts);
        pairs.map(|(&s, &c)| s / c as f64).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn atoms_flying_freely() {
        // Two atoms at constant velocity, sampled every h = 0.5. Lag k is
        // k h of flight, so the MSD is the average |v|^2 = 3.25 times
        // (k h)^2, and every v(t + k h) . v(t) is |v|^2, so the VACF is
        // 3.25 at every lag. 10 frames, lags 0 to 4. Every value is exact.
        let v = [Vec3::new(1.0, 2.0, -0.5), Vec3::new(-0.5, 0.0, 1.0)];
        let start = [Vec3::new(0.5, 7.5, 3.0), Vec3::new(2.0, 0.5, 7.5)];
        let mut msd = Correlation::msd(4);
        let mut vacf = Correlation::vacf(4);
        for n in 0..10 {
            let t = f64::from(n) * 0.5;
            msd.add(&[start[0] + v[0] * t, start[1] + v[1] * t]);
            vacf.add(&v);
        }
        let want: Vec<f64> = (0..5)
            .map(|k| 3.25 * (f64::from(k) * 0.5).powi(2))
            .collect();
        assert_eq!(msd.values(), want);
        assert_eq!(vacf.values(), vec![3.25; 5]);
    }

    #[test]
    fn every_origin_counts_once_after_the_oldest_frame_is_dropped() {
        // One atom at x = n^2 for n = 0 to 9, lags 0 to 2, so only 3
        // frames are kept. Lag 1 averages (n^2 - (n-1)^2)^2 = (2n - 1)^2
        // over n = 1 to 9: 969 / 9. Lag 2 averages (4n - 4)^2 over
        // n = 2 to 9: 3264 / 8 = 408. Sums of whole numbers are exact.
        let mut msd = Correlation::msd(2);
        for n in 0..10 {
            msd.add(&[Vec3::new(f64::from(n * n), 0.0, 0.0)]);
        }
        assert_eq!(msd.values(), vec![0.0, 969.0 / 9.0, 408.0]);
        // A velocity that flips every frame, -1, 1, -1, ...: the VACF is 1
        // at even lags and -1 at odd ones.
        let mut vacf = Correlation::vacf(2);
        for n in 0..10 {
            vacf.add(&[Vec3::new(f64::from(n % 2) * 2.0 - 1.0, 0.0, 0.0)]);
        }
        assert_eq!(vacf.values(), vec![1.0, -1.0, 1.0]);
    }

    #[test]
    #[should_panic(expected = "fewer frames than lags")]
    fn too_few_frames_are_refused() {
        let mut vacf = Correlation::vacf(3);
        vacf.add(&[Vec3::new(1.0, 0.0, 0.0)]);
        vacf.values();
    }
}