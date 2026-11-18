//! The two currents whose fluctuations give viscosity and heat conduction:
//! the stress and the heat current.
//!
//! Stress. The pressure tensor times the volume is
//!
//!   P_ab V = sum over atoms of v_a v_b + sum over pairs of r_a F_b,
//!
//! with r = r_i - r_j and F the force on i from j (each mass is 1). Its
//! trace over 3 is the pressure of thermo.rs times V. The parts off the
//! diagonal, xy, yz and zx, are the shear stress. Each is the rate of
//! change of a Helfand moment, for xy the sum of y_i v_x,i, so a stress
//! that lingers means momentum carried far across the box: viscosity.
//! The diagonal parts less their average are shear too, seen in axes
//! turned by 45 degrees (greenkubo.rs uses all five).
//!
//! Heat current. Each atom carries e_i = v_i^2 / 2 plus half the shifted
//! energy of each pair it is in. The current, summed over the box, is
//!
//!   J = sum over atoms of e_i v_i + 1/2 sum over i, j of r_ij (F_ij . v_i),
//!
//! the rate of change of the energy moment sum r_i e_i, as long as no pair
//! straddles the box edge. The first part is energy that atoms carry as
//! they move, the second energy that forces pass from atom to atom. With
//! total momentum zero no term for the average energy is needed. Viscardy,
//! Servantie and Gaspard write the same currents (arXiv:cond-mat/0701253,
//! eq. 11, and 0701254, eq. A5), with the unshifted energy in e_i, which
//! changes J by U(r_c) / 2 times the sum of n_i v_i for n_i pairs: in two
//! trial blocks of 100 tau at their state that moved the conductivity by
//! 1.4e-4 and 1.0e-3 of itself.
//! tests/flux.rs checks the currents pair by pair, and tests/helfand.rs as
//! the rates of the moments.
//!
//! Label: derived, for the cut and shifted model of pairs.rs.

use crate::lj;
use crate::neighbors::NeighborList;
use crate::parallel;
use crate::periodic::PeriodicBox;
use crate::vec3::Vec3;

const ZERO: Vec3 = Vec3::new(0.0, 0.0, 0.0);

/// The currents at one moment.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Flux {
    /// The xy, yz and zx parts of the pressure tensor times the volume,
    /// in eps.
    pub shear: Vec3,
    /// The xx, yy and zz parts of the same, in eps.
    pub normal: Vec3,
    /// The heat current J, in eps sigma / tau.
    pub heat: Vec3,
    /// Pairs closer than the cutoff. Unit: none.
    pub pairs: usize,
}

/// One atom's share from its pairs: the off-diagonal and diagonal stress,
/// the potential energy, the energy forces pass on, and the pair count.
type Share = (Vec3, Vec3, f64, Vec3, usize);

/// The currents of atoms at `positions` (sigma) with `velocities`
/// (sigma / tau), from the pairs on `list`, which must hold every pair
/// within `cutoff` (sigma), as it does after a force call. The atoms are
/// split across `threads` threads, and the sums run in atom order, so the
/// result is the same for any number of threads.
pub fn measure(
    positions: &[Vec3],
    velocities: &[Vec3],
    space: PeriodicBox,
    cutoff: f64,
    list: &NeighborList,
    threads: usize,
) -> Flux {
    assert_eq!(positions.len(), velocities.len());
    let cutoff_sq = cutoff * cutoff;
    let shift = lj::energy(Vec3::new(cutoff, 0.0, 0.0));
    let mut each: Vec<Share> = vec![(ZERO, ZERO, 0.0, ZERO, 0); positions.len()];
    parallel::for_each(&mut each, threads, |i, out| {
        let v = velocities[i];
        let (mut off, mut diagonal, mut u, mut passed, mut n) = (ZERO, ZERO, 0.0, ZERO, 0);
        for &j in list.of(i) {
            let r = space.nearest_image(positions[i] - positions[j]);
            if r.norm_sq() < cutoff_sq {
                let f = lj::force(r);
                off = off + Vec3::new(r.x * f.y, r.y * f.z, r.z * f.x);
                diagonal = diagonal + Vec3::new(r.x * f.x, r.y * f.y, r.z * f.z);
                u += lj::energy(r) - shift;
                passed = passed + r * f.dot(v);
                n += 1;
            }
        }
        *out = (off, diagonal, u, passed, n);
    });
    let (mut shear, mut normal, mut heat, mut listed) = (ZERO, ZERO, ZERO, 0);
    for (&(off, diagonal, u, passed, n), &v) in each.iter().zip(velocities) {
        // Each pair is on the list of both its atoms, hence the halves.
        shear = shear + Vec3::new(v.x * v.y, v.y * v.z, v.z * v.x) + off * 0.5;
        normal = normal + Vec3::new(v.x * v.x, v.y * v.y, v.z * v.z) + diagonal * 0.5;
        heat = heat + v * (0.5 * v.norm_sq() + 0.5 * u) + passed * 0.5;
        listed += n;
    }
    Flux {
        shear,
        normal,
        heat,
        pairs: listed / 2,
    }
}
