//! Error bars for averages, and the integral of evenly spaced samples.
//!
//! A long run is cut into blocks, the quantity is measured once in each
//! block, and the spread of the block values gives the error bar of their
//! mean. This works when each block is long compared with the time the
//! motion takes to forget itself, so that the block values are
//! independent (Flyvbjerg and Petersen, J. Chem. Phys. 91, 461 (1989)).

/// The mean of `values` and its standard error s / sqrt(n), where
/// s^2 = sum of (x - mean)^2 / (n - 1). Both are in the unit of `values`.
///
/// Label: approximated. The standard error assumes independent values.
/// With n values it is itself uncertain by about 1 part in sqrt(2 (n - 1)).
pub fn mean_and_error(values: &[f64]) -> (f64, f64) {
    assert!(values.len() >= 2, "an error bar needs two or more values");
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    let squares: f64 = values.iter().map(|x| (x - mean) * (x - mean)).sum();
    (mean, (squares / (n - 1.0) / n).sqrt())
}

/// The trapezoid rule: the integral of samples spaced `h` apart, from the
/// first sample to the last, h (y_0 / 2 + y_1 + ... + y_n-1 + y_n / 2).
/// Unit: the unit of the samples times the unit of `h`.
///
/// Label: approximated. Exact for a straight line. For a smooth curve the
/// error shrinks as h^2.
pub fn trapezoid(values: &[f64], h: f64) -> f64 {
    assert!(values.len() >= 2, "an integral needs two or more samples");
    let inner: f64 = values[1..values.len() - 1].iter().sum();
    h * (0.5 * values[0] + inner + 0.5 * values[values.len() - 1])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mean_and_error_of_two_values() {
        // 1 and 7: mean 4, s^2 = (9 + 9) / 1 = 18, error sqrt(18 / 2) = 3.
        // Every step is exact in an f64.
        assert_eq!(mean_and_error(&[1.0, 7.0]), (4.0, 3.0));
    }

    #[test]
    fn equal_values_have_no_error() {
        assert_eq!(mean_and_error(&[2.5; 6]), (2.5, 0.0));
    }

    #[test]
    fn error_shrinks_as_one_over_root_n() {
        // 1, 7 repeated 50 times: mean 4, each value 3 away, so
        // s^2 = 100 * 9 / 99 and the error is sqrt(9 / 99) = 0.3015.
        // Measured: exact. Tolerance 1e-15, for rounding, our own choice.
        let values: Vec<f64> = (0..100)
            .map(|i| if i % 2 == 0 { 1.0 } else { 7.0 })
            .collect();
        let (mean, error) = mean_and_error(&values);
        assert_eq!(mean, 4.0);
        assert!(
            (error - (9.0_f64 / 99.0).sqrt()).abs() <= 1e-15,
            "error {error}"
        );
    }

    #[test]
    fn trapezoid_is_exact_for_a_straight_line() {
        // y = 1 + 2t at t = 0, 0.5, 1, 1.5: the integral is t + t^2 at 1.5,
        // 3.75. Every value is exact in an f64.
        assert_eq!(trapezoid(&[1.0, 2.0, 3.0, 4.0], 0.5), 3.75);
    }

    #[test]
    fn trapezoid_of_a_falling_exponential() {
        // exp(-t) from 0 to 20 in steps of 0.01. The exact integral is
        // 1 - exp(-20). The trapezoid rule overshoots by h^2 / 12 times the
        // change in slope, 1 - exp(-20), about 8.3e-6, less a next term of
        // about h^4 / 720 = 1.4e-11. Measured: 8.3333e-6, which is
        // 1.39e-11 under the first term. Tolerance: 1e-10 around the first
        // term, our own choice.
        let h = 0.01;
        let values: Vec<f64> = (0..=2000).map(|k| (-(k as f64) * h).exp()).collect();
        let exact = 1.0 - (-20.0_f64).exp();
        let over = trapezoid(&values, h) - exact;
        assert!((over - h * h / 12.0 * exact).abs() <= 1e-10, "over {over}");
    }
}