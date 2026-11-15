//! Our own random number generator: SplitMix64 by Sebastiano Vigna (2015),
//! a fixed-increment version of Java 8's SplittableRandom. Reference code:
//! <https://prng.di.unimi.it/splitmix64.c>, which says it passes BigCrush.
//!
//! The state is one 64-bit integer. Each call adds a fixed odd constant to
//! it and scrambles the result. Only integer arithmetic is used, so a seed
//! gives the same numbers on every machine.

use crate::ln::ln;

/// 2^-53. Multiplying a 53-bit integer by it gives a number in [0, 1).
const TWO_TO_MINUS_53: f64 = 1.0 / (1u64 << 53) as f64;

/// A seeded random number generator.
pub struct Rng {
    state: u64,
    /// The second normal number from the last pair, kept for the next call.
    spare: Option<f64>,
}

impl Rng {
    /// A generator that starts from `seed`, a plain integer. The same seed
    /// always gives the same sequence.
    pub fn new(seed: u64) -> Self {
        Self {
            state: seed,
            spare: None,
        }
    }

    /// The next 64 random bits, as in the reference code. Unit: none.
    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    /// A number drawn evenly from [0, 1): the top 53 bits of `next_u64`
    /// times 2^-53, as the reference site recommends. Unit: none.
    pub fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * TWO_TO_MINUS_53
    }

    /// A number drawn from the standard normal distribution: mean 0,
    /// variance 1. Unit: none.
    ///
    /// Marsaglia's polar method (Marsaglia and Bray, SIAM Review 6, 260
    /// (1964)): pick a point (u, v) evenly in the square from -1 to 1 and
    /// keep it only if s = u^2 + v^2 is inside the unit circle and not 0.
    /// Then u * sqrt(-2 ln(s) / s) and v * sqrt(-2 ln(s) / s) are two
    /// independent standard normal numbers.
    pub fn normal(&mut self) -> f64 {
        if let Some(z) = self.spare.take() {
            return z;
        }
        loop {
            let u = 2.0 * self.uniform() - 1.0;
            let v = 2.0 * self.uniform() - 1.0;
            let s = u * u + v * v;
            if s > 0.0 && s < 1.0 {
                let scale = (-2.0 * ln(s) / s).sqrt();
                self.spare = Some(v * scale);
                return u * scale;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_the_reference_code() {
        // First four outputs of splitmix64.c, compiled and run with the
        // state set to 1234567.
        let mut rng = Rng::new(1234567);
        assert_eq!(rng.next_u64(), 6457827717110365317);
        assert_eq!(rng.next_u64(), 3203168211198807973);
        assert_eq!(rng.next_u64(), 9817491932198370423);
        assert_eq!(rng.next_u64(), 4593380528125082431);
    }

    #[test]
    fn same_seed_gives_same_numbers() {
        let mut a = Rng::new(42);
        let mut b = Rng::new(42);
        for _ in 0..1000 {
            assert_eq!(a.normal(), b.normal());
        }
    }

    // The tests below draw N numbers and compare averages with their exact
    // values. Each tolerance is about 4.5 standard errors, our own choice.

    const N: usize = 200_000;

    #[test]
    fn uniform_stays_in_range_with_mean_one_half() {
        // Standard error of the mean: sqrt(1/12 / N) = 0.00065.
        let mut rng = Rng::new(7);
        let mut sum = 0.0;
        for _ in 0..N {
            let x = rng.uniform();
            assert!((0.0..1.0).contains(&x));
            sum += x;
        }
        let mean = sum / N as f64;
        assert!((mean - 0.5).abs() <= 0.003, "mean {mean}");
    }

    #[test]
    fn normal_has_mean_0_variance_1_and_fourth_moment_3() {
        // For a standard normal, the average of z is 0, of z^2 is 1 and of
        // z^4 is 3. Standard errors: sqrt(1/N) = 0.0022,
        // sqrt(2/N) = 0.0032 and sqrt(96/N) = 0.022.
        let mut rng = Rng::new(11);
        let (mut m1, mut m2, mut m4) = (0.0, 0.0, 0.0);
        for _ in 0..N {
            let z = rng.normal();
            m1 += z;
            m2 += z * z;
            m4 += z * z * z * z;
        }
        let n = N as f64;
        let (m1, m2, m4) = (m1 / n, m2 / n, m4 / n);
        assert!(m1.abs() <= 0.01, "mean {m1}");
        assert!((m2 - 1.0).abs() <= 0.015, "variance {m2}");
        assert!((m4 - 3.0).abs() <= 0.1, "fourth moment {m4}");
    }

    #[test]
    fn successive_normals_are_independent() {
        // The polar method gives numbers in pairs. If the two were related,
        // the average of z_i * z_(i+1) would not be 0. Standard error 0.0022.
        let mut rng = Rng::new(13);
        let mut previous = rng.normal();
        let mut sum = 0.0;
        for _ in 0..N {
            let z = rng.normal();
            sum += previous * z;
            previous = z;
        }
        let mean = sum / N as f64;
        assert!(mean.abs() <= 0.01, "average of z_i * z_(i+1) {mean}");
    }
}
