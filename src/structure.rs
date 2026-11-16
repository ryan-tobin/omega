//! The static structure factor S(k): how strongly the liquid scatters
//! X-rays or neutrons by wavenumber k, which is what diffraction
//! experiments measure. It is the Fourier transform of g(r):
//!
//! S(k) = 1 + 4 pi rho * integral of r^2 (g(r) - 1) sin(k r) / (k r) dr.
//!
//! The peaks of S(k) come from the shells of neighbors in g(r): the first,
//! near k = 2 pi / (nearest-neighbor distance), is the main ring in a
//! diffraction picture. Rahman compared the peaks of his S(k) with X-ray
//! data (rahman.rs).
//!
//! The integral runs over the bins of g(r), each bin's value taken at its
//! middle (the midpoint rule), and stops where g(r) stops, at most L/2.
//! Label: approximated. Stopping the integral at r_max puts ripples of
//! wavelength about 2 pi / r_max on S(k), because g(r) - 1 is not yet zero
//! there.
//!
//! sin is Rust's own `f64::sin`, which may differ in the last bit from one
//! platform to another (ln.rs). That is harmless here: S(k) is worked out
//! after the run and feeds nothing back into the motion, so a difference
//! in the last bit stays in the last bit.

use std::f64::consts::PI;

/// S(k) at wavenumber `k`, in 1 / sigma, from the points (r, g) of a
/// radial distribution function (rdf.rs), r in sigma, at number density
/// `density`, in 1 / sigma^3. A plain number.
pub fn structure_factor(g: &[(f64, f64)], density: f64, k: f64) -> f64 {
    let dr = g[1].0 - g[0].0;
    let mut integral = 0.0;
    for &(r, gr) in g {
        let kr = k * r;
        integral += r * r * (gr - 1.0) * kr.sin() / kr * dr;
    }
    1.0 + 4.0 * PI * density * integral
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gaussian_bump_has_a_gaussian_transform() {
        // g(r) - 1 = A exp(-r^2 / (2 s^2)) has the transform
        // 4 pi * integral of r^2 (g - 1) sin(kr) / (kr) dr
        // = A (2 pi s^2)^(3/2) exp(-k^2 s^2 / 2), the three-dimensional
        // Fourier transform of a Gaussian. With s = 0.5 in bins of 0.01 out
        // to 5 (10 s) the midpoint rule is all but exact for so smooth a
        // curve. Measured: 1.3e-15. Tolerance: 1e-12, our own choice.
        let (a, s, density) = (0.8, 0.5, 0.9);
        let g: Vec<(f64, f64)> = (0..500)
            .map(|i| 0.01 * (f64::from(i) + 0.5))
            .map(|r| (r, 1.0 + a * (-r * r / (2.0 * s * s)).exp()))
            .collect();
        for k in [0.5, 2.0, 5.0, 9.0] {
            let want =
                1.0 + density * a * (2.0 * PI * s * s).powf(1.5) * (-k * k * s * s / 2.0).exp();
            let got = structure_factor(&g, density, k);
            assert!((got - want).abs() <= 1e-12, "S({k}) = {got}, want {want}");
        }
    }

    #[test]
    fn evenly_spread_atoms_scatter_like_no_structure() {
        // g = 1 everywhere: no structure, so S(k) = 1 at every k.
        let g: Vec<(f64, f64)> = (0..100)
            .map(|i| (0.05 * (f64::from(i) + 0.5), 1.0))
            .collect();
        assert_eq!(structure_factor(&g, 0.8, 3.0), 1.0);
    }
}