//! Rahman's liquid argon: A. Rahman, Phys. Rev. 136, A405 (1964). His
//! box, temperature and cutoff in reduced units (PHYSICS_REFERENCE), how
//! we bring the atoms there, and the numbers he published that M4
//! compares with ours.

use crate::lattice;
use crate::periodic::PeriodicBox;
use crate::rng::Rng;
use crate::system::System;
use crate::thermo;
use crate::velocities;

/// His box side, in sigma: 864 atoms at 1.374 g/cm^3 with his atomic
/// mass unit. Label: assumed, following the paper.
pub const SIDE: f64 = 10.229;
/// His temperature, 94.4 K, in eps / kB (120 K).
pub const TEMPERATURE: f64 = 0.7867;
/// His cutoff, in sigma.
pub const CUTOFF: f64 = 2.25;
/// Our time step, in tau. His was 1e-14 s, 0.00465 tau.
pub const DT: f64 = 0.005;
/// Steps of melting at his temperature: 20 tau. Our own choice, as M3.
pub const RESCALE_STEPS: usize = 4000;
/// The velocities are rescaled every this many steps while melting.
pub const RESCALE_EVERY: usize = 10;

// His paper is behind a paywall, so only its abstract was read. The
// numbers below are quoted from it by others, named with each number.
// H. Heidenreich's notes on the paper are at hunterheidenreich.pages.dev.

/// His self-diffusion coefficient at 94.4 K, in cm^2 / s. Quoted by Berne,
/// Boon and Rice, J. Chem. Phys. 45, 1086 (1966), and by Heidenreich.
pub const DIFFUSION_CM2_PER_S: f64 = 2.43e-5;
/// Where his g(r) has its first peak, in angstrom: the distance between
/// nearest neighbors that Berne, Boon and Rice read from his g(r).
/// Heidenreich gives 3.7.
pub const FIRST_PEAK_ANGSTROM: f64 = 3.67;
/// When his velocity autocorrelation function first turns negative, in
/// picoseconds. Quoted by J. C. Modla, PhD thesis, Lehigh (1977), and by
/// Heidenreich.
pub const VACF_ZERO_PS: f64 = 0.33;
/// The first four peaks of his S(k), the Fourier transform of his g(r),
/// as k sigma. Quoted by Heidenreich only.
pub const S_PEAKS: [f64; 4] = [6.8, 12.5, 18.5, 24.8];
/// The same peaks in the X-ray measurement he compared with, Eisenstein
/// and Gingrich, Phys. Rev. 62, 261 (1942), at 91.8 K, as k sigma.
/// Quoted by Heidenreich only.
pub const S_PEAKS_XRAY: [f64; 4] = [6.8, 12.3, 18.4, 24.4];
/// Where we look for each peak of S(k), as (from, to) in k sigma: our
/// own choice, about 1.5 to either side of his peak.
pub const S_PEAK_WINDOWS: [(f64, f64); 4] = [(5.5, 8.0), (11.0, 14.0), (17.0, 20.0), (23.0, 26.5)];

/// 864 atoms on an fcc lattice in Rahman's box, with velocities drawn for
/// his temperature from `seed`, stepped RESCALE_STEPS times with the
/// velocities rescaled to his temperature every RESCALE_EVERY steps. The
/// lattice melts in the first few tau, and the rest of the time holds
/// the liquid at his temperature. The last step is a rescaling, so the
/// temperature on return is his.
pub fn liquid(seed: u64) -> System {
    let space = PeriodicBox { side: SIDE };
    let positions = lattice::fcc(6, SIDE);
    let mut rng = Rng::new(seed);
    let start = velocities::maxwell_boltzmann(positions.len(), TEMPERATURE, &mut rng);
    let mut system = System::new(positions, start, space, CUTOFF);
    for step in 1..=RESCALE_STEPS {
        system.step(DT);
        if step % RESCALE_EVERY == 0 {
            thermo::rescale(&mut system.velocities, TEMPERATURE);
        }
    }
    system
}