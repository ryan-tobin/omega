//! Where a sampled curve is highest or lowest, and where it first drops
//! below a level. Used on g(r), for its first peak, first trough and second
//! peak, and on the velocity autocorrelation function, for where it first
//! turns negative and how deep it goes.
//!
//! A curve is a list of points (x, y) with x evenly spaced, as rdf.rs and
//! correlate.rs give them. The highest sample in a window is refined by
//! the parabola through it and its two neighbors, so the answer is not
//! tied to the grid. For samples a, b, c at x - h, x, x + h, the parabola
//! has its vertex at x + s h with s = (a - c) / (2 (a - 2b + c)), where
//! its value is b - s (a - c) / 4.
//!
//! Label: approximated. Exact when the curve is a parabola over three
//! samples. A noisy curve moves the vertex by about the noise divided by
//! the curvature.

/// The top of the curve between x = `from` and x = `to`: the highest
/// sample there, refined by the parabola through its neighbors. Returns
/// (x, y) of the vertex, in the units of the points.
pub fn highest(points: &[(f64, f64)], from: f64, to: f64) -> (f64, f64) {
    vertex(points, from, to, 1.0)
}

/// The bottom of the curve between x = `from` and x = `to`, found the same
/// way as `highest`. Returns (x, y) of the vertex.
pub fn lowest(points: &[(f64, f64)], from: f64, to: f64) -> (f64, f64) {
    let (x, y) = vertex(points, from, to, -1.0);
    (x, -y)
}

/// The highest point of `sign` * y, with sign 1 or -1.
fn vertex(points: &[(f64, f64)], from: f64, to: f64, sign: f64) -> (f64, f64) {
    let mut best: Option<usize> = None;
    for (i, &(x, y)) in points.iter().enumerate() {
        let higher = best.is_none_or(|b| sign * y > sign * points[b].1);
        if x >= from && x <= to && higher {
            best = Some(i);
        }
    }
    let i = best.expect("no point between from and to");
    assert!(
        i > 0 && i + 1 < points.len(),
        "extreme at the end of the curve"
    );
    let (a, b, c) = (
        sign * points[i - 1].1,
        sign * points[i].1,
        sign * points[i + 1].1,
    );
    let h = points[i + 1].0 - points[i].0;
    let s = (a - c) / (2.0 * (a - 2.0 * b + c));
    (points[i].0 + s * h, b - s * (a - c) / 4.0)
}

/// The first x at which the curve drops below `level`, on the straight
/// line between the last sample at or above it and the first below it.
/// None if it never does. Unit: the unit of x.
pub fn first_below(points: &[(f64, f64)], level: f64) -> Option<f64> {
    for pair in points.windows(2) {
        let ((x0, y0), (x1, y1)) = (pair[0], pair[1]);
        if y0 >= level && y1 < level {
            return Some(x0 + (y0 - level) / (y0 - y1) * (x1 - x0));
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Samples of `f` at x = 0, 0.1, ..., 3.
    fn curve(f: impl Fn(f64) -> f64) -> Vec<(f64, f64)> {
        (0..=30)
            .map(|k| f64::from(k) * 0.1)
            .map(|x| (x, f(x)))
            .collect()
    }

    #[test]
    fn parabola_vertex_is_found_off_the_grid() {
        // y = 3 - 2 (x - 1.23)^2 peaks at 1.23 with height 3, between the
        // samples at 1.2 and 1.3. Tolerance 1e-12, for rounding, our own
        // choice.
        let points = curve(|x| 3.0 - 2.0 * (x - 1.23) * (x - 1.23));
        let (x, y) = highest(&points, 0.5, 2.0);
        assert!(
            (x - 1.23).abs() <= 1e-12 && (y - 3.0).abs() <= 1e-12,
            "{x} {y}"
        );
        let points = curve(|x| (x - 0.37) * (x - 0.37) - 1.0);
        let (x, y) = lowest(&points, 0.1, 1.0);
        assert!(
            (x - 0.37).abs() <= 1e-12 && (y + 1.0).abs() <= 1e-12,
            "{x} {y}"
        );
    }

    #[test]
    fn window_picks_the_second_peak() {
        // cos(2 pi x) peaks at x = 0, 1, 2, 3 and is lowest at 0.5, 1.5,
        // 2.5. The window 1.6 to 2.4 holds only the peak at 2. Near a peak
        // cos is not a parabola, but the peak sits on a sample and the
        // neighbors are equal, so s = 0 and the sample itself is returned.
        let points = curve(|x| (2.0 * std::f64::consts::PI * x).cos());
        let (x, y) = highest(&points, 1.6, 2.4);
        assert!(
            (x - 2.0).abs() <= 1e-12 && (y - 1.0).abs() <= 1e-12,
            "{x} {y}"
        );
    }

    #[test]
    fn a_higher_point_outside_the_window_is_ignored() {
        // Bumps of height 1 at x = 1 and 5, and of height 5 at x = 3. Each
        // top has equal neighbors, so s = 0 and the sample itself is
        // returned. Every value is exact.
        let heights = [0.0, 1.0, 0.0, 5.0, 0.0, 1.0, 0.0];
        let points: Vec<(f64, f64)> = heights
            .iter()
            .enumerate()
            .map(|(i, &y)| (i as f64, y))
            .collect();
        assert_eq!(highest(&points, 0.5, 2.5), (1.0, 1.0));
        assert_eq!(highest(&points, 3.5, 5.5), (5.0, 1.0));
    }

    #[test]
    #[should_panic(expected = "extreme at the end of the curve")]
    fn an_extreme_on_the_first_sample_is_refused() {
        // With no sample before it there is no parabola to fit.
        highest(&[(0.0, 5.0), (1.0, 1.0), (2.0, 0.0)], 0.0, 2.0);
    }

    #[test]
    fn first_crossing_below_a_level() {
        // y = 1 - x drops below 0.25 at x = 0.75, between the samples at
        // 0.7 and 0.8, on the straight line between them. Measured: exact.
        // Tolerance 1e-12, for rounding. It never drops below -5 here.
        let points = curve(|x| 1.0 - x);
        let x = first_below(&points, 0.25).unwrap();
        assert!((x - 0.75).abs() <= 1e-12, "{x}");
        assert_eq!(first_below(&points, -5.0), None);
        // A curve that starts below the level counts only once it has
        // risen above it: from 1 down to -1, it crosses 0 at x = 2.5.
        let points = [(0.0, -2.0), (1.0, -1.0), (2.0, 1.0), (3.0, -1.0)];
        assert_eq!(first_below(&points, 0.0), Some(2.5));
    }
}