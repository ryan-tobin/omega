//! The Green-Kubo integrals of greenkubo.rs on currents with known
//! answers.

use omega::flux::Flux;
use omega::greenkubo::GreenKubo;
use omega::vec3::Vec3;

const ZERO: Vec3 = Vec3::new(0.0, 0.0, 0.0);

fn flux(shear: Vec3, normal: Vec3, heat: Vec3) -> Flux {
    Flux {
        shear,
        normal,
        heat,
        pairs: 0,
    }
}

#[test]
fn steady_currents_integrate_to_straight_lines() {
    // Currents that never change correlate fully at every lag. Shear
    // (1, 2, 2): S . S = 9. Diagonal (3, 0, 0), less its third, is
    // (2, -1, -1): D . D = 6. Heat (1, 1, 1): J . J = 3. With V = 2,
    // T = 0.5 and h = 0.1 the running integrals are (2 * 9 + 6) t / 10
    // = 2.4 t for the viscosity, 3 t and 1.5 t for its two parts, and
    // (3 / 3) t / (2 * 0.25) = 2 t for the conductivity. Tolerance
    // 1e-14, for rounding.
    let mut gk = GreenKubo::new(4, 0.1, 2.0);
    let steady = flux(
        Vec3::new(1.0, 2.0, 2.0),
        Vec3::new(3.0, 0.0, 0.0),
        Vec3::new(1.0, 1.0, 1.0),
    );
    for _ in 0..10 {
        gk.add(&steady, 0.5);
    }
    assert_eq!(gk.temperature(), 0.5);
    let (off, normal) = gk.viscosity_parts();
    let lines = [
        (gk.viscosity(), 2.4),
        (off, 3.0),
        (normal, 1.5),
        (gk.conductivity(), 2.0),
    ];
    for (curve, slope) in lines {
        for (k, value) in curve.iter().enumerate() {
            let want = slope * k as f64 * 0.1;
            assert!(
                (value - want).abs() <= 1e-14,
                "lag {k}: {value}, want {want}"
            );
        }
    }
}

#[test]
fn a_current_that_flips_each_sample_integrates_to_nothing() {
    // Shear (1, 0, 0), then (-1, 0, 0), and so on: the correlation is
    // 1, -1, 1, ... The trapezoid sums with h = 1 are 0, then
    // 1/2 - 1/2 = 0, then 1/2 - 1 + 1/2 = 0, and so on. Exact.
    let mut gk = GreenKubo::new(3, 1.0, 1.0);
    for n in 0..8 {
        let s = if n % 2 == 0 { 1.0 } else { -1.0 };
        gk.add(&flux(Vec3::new(s, 0.0, 0.0), ZERO, ZERO), 1.0);
    }
    assert_eq!(gk.viscosity(), vec![0.0; 4]);
}
