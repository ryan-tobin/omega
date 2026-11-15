//! The natural logarithm, written with +, -, * and / only.
//!
//! Rust's own `f64::ln` is allowed to differ in the last bit from one
//! platform to another: its documentation says its precision is not
//! specified. That matters here because many atoms moving together are
//! chaotic. A one-bit difference in one starting velocity grows until two
//! runs have nothing in common. IEEE 754 fixes the result of +, -, * and /
//! on every machine, so with this function the same seed gives the same
//! run on a Mac and on Linux.

use std::f64::consts::{LN_2, SQRT_2};

/// How many terms of the series to add. With |t| below 0.172, the first
/// term left out is under 1e-18 of the first term kept, far below what an
/// f64 can resolve.
const TERMS: u32 = 11;

/// Natural logarithm of x, a plain number. x must be positive, finite and
/// at least 2^-1022 (the smallest normal f64). Unit: none.
pub fn ln(x: f64) -> f64 {
    assert!((f64::MIN_POSITIVE..f64::INFINITY).contains(&x), "ln of {x}");
    // An f64 is stored as m * 2^k with m in [1, 2). Read k from the
    // exponent bits, then give m the exponent bits of 1.0.
    let bits = x.to_bits();
    let mut k = ((bits >> 52) & 0x7ff) as i64 - 1023;
    let mut m = f64::from_bits((bits & 0x000f_ffff_ffff_ffff) | 0x3ff0_0000_0000_0000);
    // The series converges fastest for m near 1, so move m into
    // [sqrt(1/2), sqrt(2)). Halving m adds one to k.
    if m >= SQRT_2 {
        m *= 0.5;
        k += 1;
    }
    // ln(x) = ln(m) + k ln(2), and ln(m) = 2 atanh(t) with
    // t = (m - 1) / (m + 1). The series is atanh(t) = t + t^3/3 + t^5/5 + ...
    // It is added in Horner form, smallest term first, to limit rounding.
    let t = (m - 1.0) / (m + 1.0);
    let t2 = t * t;
    let mut sum = 0.0;
    for n in (0..TERMS).rev() {
        sum = sum * t2 + 1.0 / f64::from(2 * n + 1);
    }
    2.0 * t * sum + k as f64 * LN_2
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ln_of_one_is_zero() {
        assert_eq!(ln(1.0), 0.0);
    }

    #[test]
    fn ln_of_powers_of_two_is_k_ln_2() {
        // 2^k has m = 1, so t = 0 and only the k ln(2) part is left. This
        // checks that k is read from the bits correctly.
        for k in -1000..1000 {
            let x = f64::from_bits(((k + 1023) as u64) << 52);
            assert_eq!(ln(x), f64::from(k) * LN_2);
        }
    }

    #[test]
    fn ln_of_e_is_one() {
        // E as an f64 is 1.4e-16 below e, so the exact answer, 1 - 5.3e-17,
        // rounds to 1.0. The tolerance, 2e-16, is one unit in the last place
        // just below 1.0, about 1.1e-16, with a margin. Our own choice.
        let miss = (ln(std::f64::consts::E) - 1.0).abs();
        assert!(miss <= 2e-16, "off by {miss:e}");
    }

    #[test]
    fn ln_agrees_with_the_standard_library() {
        // Rust's ln is within about one unit in the last place of the true
        // value on common platforms, so it serves as a reference here. The
        // sweep covers 1e-300 to 1e300 and, most densely, (0, 1), where the
        // random number generator uses ln. The tolerance is our own choice,
        // 1e-15 relative: 4.5 to 9 units in the last place, depending on the
        // value. Measured worst: 3.9e-16, 2 units in the last place.
        let mut worst: f64 = 0.0;
        for i in 1..20_000 {
            let below_one = f64::from(i) / 20_000.0;
            let wide = 10.0_f64.powf(f64::from(i - 10_000) * 0.03);
            for x in [below_one, wide] {
                let want = x.ln();
                let miss = (ln(x) - want).abs() / want.abs().max(1e-300);
                worst = worst.max(miss);
            }
        }
        assert!(worst <= 1e-15, "worst relative miss {worst:e}");
    }
}
