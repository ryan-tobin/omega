//! A measurement made on moving atoms: the mean squared displacement of
//! atoms that cross the faces of the box. The tolerance is our own choice
//! and says why.

use omega::msd::Msd;
use omega::periodic::PeriodicBox;
use omega::system::System;
use omega::vec3::Vec3;

#[test]
fn msd_counts_the_whole_path_across_the_box() {
    // Two atoms far apart in a box of side 8, both moving at the same
    // velocity, so they never come within the cutoff and nothing pushes
    // them. After 1000 steps of 0.01 tau, t = 10 and each has gone
    // (10, 5, 2.5) sigma, out through a face. Positions are never wrapped,
    // so MSD = |v|^2 t^2 = 1.3125 * 100 = 131.25. Rounding in the 1000
    // position updates is measured at 2.5e-12. Tolerance: 1e-9.
    let v = Vec3::new(1.0, 0.5, 0.25);
    let positions = vec![Vec3::new(1.0, 1.0, 1.0), Vec3::new(5.0, 5.0, 5.0)];
    let mut system = System::new(positions, vec![v, v], PeriodicBox { side: 8.0 }, 2.5);
    let msd = Msd::new(&system.positions);
    for _ in 0..1000 {
        system.step(0.01);
    }
    let m = msd.measure(&system.positions);
    assert!((m - 131.25).abs() <= 1e-9, "MSD {m}");
}
