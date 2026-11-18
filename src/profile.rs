//! The density along x of a slab of liquid with vapor on both sides.
//!
//! The box is cut into slices of equal width along x, and the atoms in
//! each slice are counted. A slab wanders slowly along x, and adding up
//! slices of a wandering slab would smear its faces. So each sample is
//! first shifted to put the middle of the slab at the middle of the box:
//! the half of the box with the fewest atoms is found, the box is cut open
//! at the middle of that half, where the vapor is, and the average x of
//! the atoms measured from the cut gives the middle of the slab.
//!
//! Label: derived. It is a count. Only +, -, * and / are used, so the same
//! atoms give the same profile on every machine.

use crate::periodic::PeriodicBox;
use crate::vec3::Vec3;

/// The number density, in 1 / sigma^3, in each of `slices` slices of equal
/// width along x, with the slab moved to the middle of the box. Slice k
/// runs from k L / slices to (k + 1) L / slices. `positions` in sigma.
pub fn along_x(positions: &[Vec3], space: PeriodicBox, slices: usize) -> Vec<f64> {
    let l = space.side;
    let width = l / slices as f64;
    let slice = |x: f64| ((x.rem_euclid(l) / width) as usize).min(slices - 1);
    let count = |shift: f64| {
        let mut counts = vec![0usize; slices];
        for r in positions {
            counts[slice(r.x - shift)] += 1;
        }
        counts
    };
    // The emptiest half: `half` slices in a row, around the ring.
    let counts = count(0.0);
    let half = slices / 2;
    let in_half = |start: usize| {
        (0..half)
            .map(|k| counts[(start + k) % slices])
            .sum::<usize>()
    };
    let start = (0..slices).min_by_key(|&s| in_half(s)).unwrap();
    let cut = (start as f64 + 0.5 * half as f64) * width;
    // The average x measured from the cut is the middle of the slab.
    let from_cut: f64 = positions.iter().map(|r| (r.x - cut).rem_euclid(l)).sum();
    let middle = cut + from_cut / positions.len() as f64;
    let volume = width * l * l;
    let centered = count(middle - 0.5 * l);
    centered.into_iter().map(|c| c as f64 / volume).collect()
}

/// The mean of a centered profile over its middle tenth, and over the
/// fifth farthest from the middle, a tenth at each end. For a slab of
/// liquid in its vapor these are the densities of the liquid and of the
/// vapor, if both reach that far from the two faces of the slab. The
/// vapor gets the wider part because it holds far fewer atoms to count.
pub fn plateaus(profile: &[f64]) -> (f64, f64) {
    let (n, k) = (profile.len(), profile.len() / 20);
    let mean = |slices: &[f64]| slices.iter().sum::<f64>() / slices.len() as f64;
    let ends = [&profile[..2 * k], &profile[n - 2 * k..]].concat();
    (mean(&profile[n / 2 - k..n / 2 + k]), mean(&ends))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Rng;

    #[test]
    fn slab_across_the_box_edge_lands_in_the_middle() {
        // Box of side 10, 10 slices 1 wide. Four planes of 9 atoms at
        // x = 8.5, 9.5, 0.5 and 1.5: a slab around the edge of the box.
        // Slices 2 to 6 are the emptiest half, the cut is at 4.5, the
        // atoms are 4, 5, 6 and 7 from it, so the slab's middle is at 10
        // and it moves by 5: to slices 3, 4, 5 and 6. Each holds 9 atoms
        // in a volume of 1 * 10 * 10. Exact.
        let mut positions = Vec::new();
        for x in [8.5, 9.5, 0.5, 1.5] {
            for k in 0..9 {
                positions.push(Vec3::new(x, f64::from(k), f64::from(k % 3)));
            }
        }
        let profile = along_x(&positions, PeriodicBox { side: 10.0 }, 10);
        let want = [0.0, 0.0, 0.0, 0.09, 0.09, 0.09, 0.09, 0.0, 0.0, 0.0];
        assert_eq!(profile, want);
    }

    #[test]
    fn profile_does_not_depend_on_where_the_slab_sits() {
        // 400 atoms at random in a slab 4 wide, in a box of side 12 with
        // 24 slices. Moving every atom by the same amount along x, even
        // out of the box, must give the same profile: the same counts,
        // unless an atom sits within rounding of a slice edge (none does
        // for this seed). The profile adds up to 400 atoms.
        let space = PeriodicBox { side: 12.0 };
        let mut rng = Rng::new(3);
        let mut random = |scale: f64| scale * rng.uniform();
        let slab: Vec<Vec3> = (0..400)
            .map(|_| Vec3::new(random(4.0), random(12.0), random(12.0)))
            .collect();
        let first = along_x(&slab, space, 24);
        let atoms: f64 = first.iter().map(|rho| rho * 0.5 * 144.0).sum();
        assert!((atoms - 400.0).abs() <= 1e-9, "{atoms} atoms");
        for shift in [3.7, -8.2, 29.31] {
            let moved: Vec<Vec3> = slab
                .iter()
                .map(|&r| r + Vec3::new(shift, 0.0, 0.0))
                .collect();
            assert_eq!(along_x(&moved, space, 24), first, "moved by {shift}");
        }
    }

    #[test]
    fn plateaus_are_the_middle_and_the_two_ends() {
        // 40 slices, slice k holding k^2. The middle tenth is slices 18
        // to 21: (324 + 361 + 400 + 441) / 4 = 381.5. The ends are slices
        // 0 to 3 and 36 to 39: (0 + 1 + 4 + 9 + 1296 + 1369 + 1444 + 1521)
        // / 8 = 705.5. Exact.
        let profile: Vec<f64> = (0..40).map(|k| f64::from(k * k)).collect();
        assert_eq!(plateaus(&profile), (381.5, 705.5));
    }
}
