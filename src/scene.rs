//! The atoms of one frame, drawn into a picture.
//!
//! Each atom is wrapped into the box, the box is turned in front of the
//! eye by the camera, and the atom is drawn as a ball (picture.rs). The
//! picture is scaled so that the box fits whichever way it is turned:
//! its corners are sqrt(3) L / 2 from its middle.
//!
//! Color shows speed. An atom at rest is blue, one at the root-mean-square
//! speed is white, and one at twice that speed or more is red, with even
//! steps in between. In equilibrium the speeds follow the Maxwell-Boltzmann
//! distribution, whose integral above the root-mean-square speed is
//! erfc(sqrt(3/2)) + sqrt(6/pi) exp(-3/2) = 0.3916 (worked out by hand and
//! checked by adding up the distribution). So 39% of the atoms should be
//! redder than white, whatever the temperature (tested).
//!
//! A ball also gets darker with distance, down to half its color at the
//! farthest corner of the box, so that near and far can be told apart.
//!
//! Label: assumed. The radius, the colors and the darkening are drawing
//! choices, our own. The camera's sines and cosines are Rust's, which may
//! differ in the last bit between machines. Nothing here feeds back into
//! the motion.

use crate::periodic::PeriodicBox;
use crate::picture::Picture;
use crate::replay::Snapshot;
use crate::vec3::Vec3;

/// Radius of a drawn atom, in sigma. Two balls touch when their atoms
/// are one sigma apart, where the pair energy crosses zero (lj.rs).
pub const RADIUS: f64 = 0.5;
/// The share of the picture's shorter side that the box may fill.
const FILL: f64 = 0.9;
/// The share of its color a ball keeps at the farthest corner of the box.
const FAR: f64 = 0.5;
/// The colors of an atom at rest, at the root-mean-square speed, and at
/// twice that speed or more, as red, green and blue.
const SLOW: [u8; 3] = [32, 96, 240];
const MIDDLE: [u8; 3] = [240, 240, 240];
const FAST: [u8; 3] = [240, 64, 32];

/// From where the box is seen: turned by `yaw` about the upright axis,
/// then tipped by `pitch` about the level one, both in radians.
pub struct Camera {
    pub yaw: f64,
    pub pitch: f64,
}

impl Camera {
    /// `r` as the eye sees it: x to the right, y up, z toward the eye,
    /// in the unit of `r`. A turn keeps every length.
    pub fn turn(&self, r: Vec3) -> Vec3 {
        let (sin_yaw, cos_yaw) = self.yaw.sin_cos();
        let (sin_pitch, cos_pitch) = self.pitch.sin_cos();
        let x = cos_yaw * r.x + sin_yaw * r.z;
        let z = cos_yaw * r.z - sin_yaw * r.x;
        let y = cos_pitch * r.y - sin_pitch * z;
        Vec3::new(x, y, sin_pitch * r.y + cos_pitch * z)
    }
}

/// The color `t` of the way from `from` to `to`, t from 0 to 1.
fn blend(from: [u8; 3], to: [u8; 3], t: f64) -> [u8; 3] {
    let mix = |k: usize| from[k] as f64 + t * (to[k] as f64 - from[k] as f64);
    [0, 1, 2].map(|k| mix(k).round() as u8)
}

/// The color of an atom moving at `speed`, where `full` is the speed
/// that is drawn red, both in sigma / tau. Half of `full` is white, and
/// an atom at rest is blue whatever `full` is.
pub fn color_of_speed(speed: f64, full: f64) -> [u8; 3] {
    if speed == 0.0 {
        return SLOW;
    }
    let s = (2.0 * speed / full).min(2.0);
    if s < 1.0 {
        blend(SLOW, MIDDLE, s)
    } else {
        blend(MIDDLE, FAST, s - 1.0)
    }
}

/// Twice the root-mean-square speed of `velocities`, in sigma / tau: the
/// speed that `color_of_speed` draws red.
pub fn full_speed(velocities: &[Vec3]) -> f64 {
    let sum: f64 = velocities.iter().map(|v| v.norm_sq()).sum();
    2.0 * (sum / velocities.len() as f64).sqrt()
}

/// Draws the atoms of `frame`, in the box `space`, as `camera` sees them,
/// with the middle of the box at the middle of `picture`. Each atom is
/// colored by its speed, with `full`, in sigma / tau, drawn red.
pub fn draw(
    picture: &mut Picture,
    frame: &Snapshot,
    space: PeriodicBox,
    camera: &Camera,
    full: f64,
) {
    let l = space.side;
    let reach = 0.75f64.sqrt() * l;
    let shorter = picture.width.min(picture.height) as f64;
    let scale = FILL * 0.5 * shorter / (reach + RADIUS);
    let (cx, cy) = (0.5 * picture.width as f64, 0.5 * picture.height as f64);
    let middle = Vec3::new(0.5 * l, 0.5 * l, 0.5 * l);
    for (r, v) in frame.positions.iter().zip(&frame.velocities) {
        let wrapped = Vec3::new(r.x.rem_euclid(l), r.y.rem_euclid(l), r.z.rem_euclid(l));
        let p = camera.turn(wrapped - middle);
        let near = FAR + (1.0 - FAR) * 0.5 * (1.0 + p.z / reach);
        let color = color_of_speed(v.norm(), full).map(|c| (c as f64 * near).round() as u8);
        let (x, y) = (cx + scale * p.x, cy - scale * p.y);
        picture.ball(x, y, scale * p.z, scale * RADIUS, color);
    }
}