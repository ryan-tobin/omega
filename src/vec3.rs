//! Three-component vectors: positions, separations, velocities and forces.
//!
//! A vector carries whatever unit its components are given in. Inside the
//! engine that is a reduced unit (see units.rs): sigma for a position,
//! eps / sigma for a force.
//!
//! Label: assumed. Space is flat and three-dimensional, so the length of a
//! vector follows Pythagoras: |a|^2 = x^2 + y^2 + z^2.

use std::ops::{Add, Mul, Neg, Sub};

/// A vector with components x, y and z, all in the same unit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    /// Builds a vector from its three components.
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Dot product a . b = ax*bx + ay*by + az*bz.
    /// Unit: the unit of `self` times the unit of `other`.
    pub fn dot(self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Squared length |a|^2 = a . a. Unit: the unit of `self`, squared.
    pub fn norm_sq(self) -> f64 {
        self.dot(self)
    }

    /// Length |a| = sqrt(a . a). Unit: the unit of `self`
    pub fn norm(self) -> f64 {
        self.norm_sq().sqrt()
    }
}

/// a + b, component by component. Both must be in the same unit
impl Add for Vec3 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }
}

/// a - b, component by component. Both must be in the same unit.
/// For two positions, r_i - r_j is the separation pointing from j to i.
impl Sub for Vec3 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

/// -a: the same length, pointing the opposite way.
impl Neg for Vec3 {
    type Output = Self;

    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }
}

/// a * s: every component multiplied by the plain number s
impl Mul<f64> for Vec3 {
    type Output = Self;

    fn mul(self, s: f64) -> Self {
        Self::new(self.x * s, self.y * s, self.z * s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Every input and every expected value here is a whole number or a
    // half, which an f64 holds exactly, and Rust's +, -, * and sqrt return
    // the correctly rounded result. So the tests compare with == and the
    // tolerance is zero.
    const A: Vec3 = Vec3::new(1.0, 2.0, 3.0);
    const B: Vec3 = Vec3::new(4.0, -5.0, 7.0);
    const C: Vec3 = Vec3::new(3.0, 4.0, 12.0);

    #[test]
    fn add_sums_each_component() {
        assert_eq!(A + B, Vec3::new(5.0, -3.0, 10.0));
    }

    #[test]
    fn sub_takes_each_component_away() {
        assert_eq!(A - B, Vec3::new(-3.0, 7.0, -4.0));
    }

    #[test]
    fn neg_flips_each_component() {
        assert_eq!(-B, Vec3::new(-4.0, 5.0, -7.0));
    }

    #[test]
    fn mul_scales_each_component() {
        assert_eq!(B * 0.5, Vec3::new(2.0, -2.5, 3.5));
    }

    #[test]
    fn dot_of_a_and_b_is_15() {
        // 1*4 + 2*(-5) + 3*7 = 4 - 10 + 21 = 15, in either order.
        assert_eq!(A.dot(B), 15.0);
        assert_eq!(B.dot(A), 15.0);
    }

    #[test]
    fn dot_of_perpendicular_vectors_is_zero() {
        // (1, 2, 3) . (3, 0, -1) = 3 + 0 - 3 = 0.
        assert_eq!(A.dot(Vec3::new(3.0, 0.0, -1.0)), 0.0);
    }

    #[test]
    fn norm_sq_of_3_4_12_is_169() {
        // 9 + 16 + 144 = 169.
        assert_eq!(C.norm_sq(), 169.0);
    }

    #[test]
    fn norm_of_3_4_12_is_13() {
        // sqrt(169) = 13.
        assert_eq!(C.norm(), 13.0);
    }
}
