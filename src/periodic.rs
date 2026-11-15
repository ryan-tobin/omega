//! A cubic box of side L that repeats in every direction.
//!
//! An atom that leaves through one face comes back through the opposite
//! face, so there are no walls. Each atom feels the nearest copy of every
//! other atom: for each component d of a separation, the nearest copy is at
//! d - L * round(d / L), which lies between -L/2 and L/2. In a cube this
//! always finds the nearest copy. But if the force cutoff were over L/2, a
//! second copy of the same atom could also be inside it and would be
//! missed. So the cutoff must be at most L/2.
//!
//! Label: approximated. A few hundred atoms in a repeating box stand in for
//! a large amount of liquid. Nothing longer than the box can exist in it,
//! such as a sound wave longer than L.

use crate::vec3::Vec3;

/// The repeating box. `side` is L, in sigma.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PeriodicBox {
    pub side: f64,
}

impl PeriodicBox {
    /// The separation `d` moved to its nearest copy, so that each component
    /// lies in [-L/2, L/2]. Unit: sigma, the same as `d`.
    pub fn nearest_image(self, d: Vec3) -> Vec3 {
        let l = self.side;
        Vec3::new(
            d.x - l * (d.x / l).round(),
            d.y - l * (d.y / l).round(),
            d.z - l * (d.z / l).round(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // A box of side 8. Every number below is a whole number of eighths,
    // which an f64 holds exactly, so the tests compare with ==.
    const BOX: PeriodicBox = PeriodicBox { side: 8.0 };

    #[test]
    fn short_separation_is_unchanged() {
        let d = Vec3::new(1.5, -2.0, 3.5);
        assert_eq!(BOX.nearest_image(d), d);
    }

    #[test]
    fn atoms_near_opposite_faces_are_close() {
        // Atoms at x = 0.5 and x = 7.5. Inside the box they are 7 apart,
        // but the copy of the second one at 7.5 - 8 = -0.5 is only 1 away.
        let d = Vec3::new(0.5, 0.0, 0.0) - Vec3::new(7.5, 0.0, 0.0);
        assert_eq!(BOX.nearest_image(d), Vec3::new(1.0, 0.0, 0.0));
    }

    #[test]
    fn separations_several_boxes_long_come_back() {
        // 17 = 2 * 8 + 1 and -12.5 = -2 * 8 + 3.5.
        let d = Vec3::new(17.0, -12.5, 0.0);
        assert_eq!(BOX.nearest_image(d), Vec3::new(1.0, 3.5, 0.0));
    }

    #[test]
    fn every_component_ends_within_half_a_box() {
        for i in -100..=100 {
            let d = f64::from(i) * 0.25;
            let e = BOX.nearest_image(Vec3::new(d, -d, 0.5 * d));
            for c in [e.x, e.y, e.z] {
                assert!(c.abs() <= 4.0, "{d} gave {c}");
            }
        }
    }
}
