//! Viscosity and heat conduction from the fluctuations of an equilibrium
//! run: the Green-Kubo integrals.
//!
//! In a liquid left alone the shear stress and the heat current of
//! flux.rs swing around zero. How long a swing lasts says how slowly the
//! liquid relaxes it, and so how much it resists flow or passes heat:
//!
//!   viscosity    eta    = 1 / (V T)   integral of avg S(0) S(t) dt,
//!   conductivity lambda = 1 / (V T^2) integral of avg J_x(0) J_x(t) dt,
//!
//! with avg the average over the run, S one off-diagonal part of the
//! pressure tensor times V, and J the heat current summed over the box
//! (kB = 1). M. S. Green, J. Chem. Phys. 22, 398 (1954) and R. Kubo, J.
//! Phys. Soc. Jpn. 12, 570 (1957) are credited with them (cited from
//! memory, not opened). Viscardy, Servantie and Gaspard write them in
//! this form (arXiv:cond-mat/0701253, eq. 2, and 0701254, eq. 1).
//! correlate.rs keeps each current as a frame of one vector and averages
//! over every time origin.
//!
//! A box with no favored direction gives more than one sample per moment.
//! The dot product of J with an earlier J is three times one part. For
//! the stress, that of S = (S_xy, S_yz, S_zx) is three times one part,
//! and the diagonal parts less a third of the trace, D_xx = P_xx V -
//! (P_xx + P_yy + P_zz) V / 3 and so on, are shear seen in turned axes:
//! the dot product of D = (D_xx, D_yy, D_zz) is four times one part
//! (derived from the form an isotropic average of four indices must take;
//! tests/transport.rs checks that S and D give the same viscosity). So
//! 2 S . S + D . D is ten times one part. P. J. Daivis and D. J. Evans,
//! J. Chem. Phys. 100, 541 (1994) are credited with using all five
//! (cited from memory, not opened).
//!
//! The integral is stopped at a time t_max, once the correlation has died
//! away. Stopped later, it gathers more noise and no more signal.
//!
//! Label: approximated. A finite run, a finite t_max and a finite box.
//! Units: eta in (m eps)^(1/2) / sigma^2, lambda in kB (eps / m)^(1/2) /
//! sigma^2.

use crate::correlate::Correlation;
use crate::flux::Flux;
use crate::stats;
use crate::vec3::Vec3;

/// The correlations of one stretch of run, and its temperature.
pub struct GreenKubo {
    shear: Correlation,
    normal: Correlation,
    heat: Correlation,
    h: f64,
    volume: f64,
    temperatures: f64,
    samples: usize,
}

impl GreenKubo {
    /// Correlations out to `lags` samples, each `h` tau apart, of a box of
    /// volume `volume`, in sigma^3.
    pub fn new(lags: usize, h: f64, volume: f64) -> Self {
        Self {
            shear: Correlation::vacf(lags),
            normal: Correlation::vacf(lags),
            heat: Correlation::vacf(lags),
            h,
            volume,
            temperatures: 0.0,
            samples: 0,
        }
    }

    /// Adds the currents of the next sample, taken h after the last, and
    /// the temperature at that moment, in eps / kB.
    pub fn add(&mut self, flux: &Flux, temperature: f64) {
        let n = flux.normal;
        let third = (n.x + n.y + n.z) / 3.0;
        self.shear.add(&[flux.shear]);
        self.normal.add(&[n - Vec3::new(third, third, third)]);
        self.heat.add(&[flux.heat]);
        self.temperatures += temperature;
        self.samples += 1;
    }

    /// The average temperature of the samples, in eps / kB.
    pub fn temperature(&self) -> f64 {
        self.temperatures / self.samples as f64
    }

    /// The viscosity with the integral stopped at each lag, 0 to `lags`,
    /// from all five parts of the stress.
    pub fn viscosity(&self) -> Vec<f64> {
        let scale = 1.0 / (10.0 * self.volume * self.temperature());
        let (s, d) = (self.shear.values(), self.normal.values());
        let both: Vec<f64> = s.iter().zip(&d).map(|(s, d)| 2.0 * s + d).collect();
        running(&both, self.h, scale)
    }

    /// The same from the three-off diagonal parts along, and from the
    /// three diagonal parts alone.
    pub fn viscosity_parts(&self) -> (Vec<f64>, Vec<f64>) {
        let scale = 1.0 / (self.volume * self.temperature());
        let off = running(&self.shear.values(), self.h, scale / 3.0);
        (off, running(&self.normal.values(), self.h, scale / 4.0))
    }

    /// The thermal conductivity with the integral stopped at each lag.
    pub fn conductivity(&self) -> Vec<f64> {
        let t = self.temperature();
        let scale = 1.0 / (3.0 * self.volume * t * t);
        running(&self.heat.values(), self.h, scale)
    }
}

/// The trapezoid integral of `values`, spaced `h` apart, from the first
/// sample to each sample in turn, times `scale`.
fn running(values: &[f64], h: f64, scale: f64) -> Vec<f64> {
    let mut sums = vec![0.0];
    for k in 1..values.len() {
        sums.push(scale * stats::trapezoid(&values[..=k], h));
    }
    sums
}
