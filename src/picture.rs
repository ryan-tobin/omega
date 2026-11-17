//! A picture: a grid of pixels that drawing fills in and a window shows.
//!
//! Each pixel is one u32 holding red, green and blue, 8 bits each, as
//! 0x00RRGGBB, row by row from the top left. That is the layout the
//! softbuffer crate asks for (softbuffer 0.4.8, src/lib.rs, "Data
//! representation"), so the viewer copies a picture into its window as
//! it is.
//!
//! An atom is drawn as a lit ball seen straight on from far away, so its
//! size on screen does not change with its distance. A ball of radius R
//! covers the pixels less than R from its center. At a pixel a distance d
//! from the center, its surface is sqrt(R^2 - d^2) nearer to the eye than
//! its center is. The picture keeps, for each pixel, how near the nearest
//! surface drawn there so far is. A ball is drawn only where it is nearer
//! than that, so a ball behind another stays hidden whichever of the two
//! is drawn first.
//!
//! The light comes from the upper left, in front: from the direction
//! (-1, -2, 2) / 3 in (right, down, toward the eye), a vector of length 1.
//! A point of the surface that faces the light head on gets the full
//! color, and one that faces away from it a quarter of it:
//! light = 1/4 + 3/4 max(0, n . l), with n the outward direction of the
//! surface there.
//!
//! Label: assumed. An atom has no surface and no color. The ball and its
//! shading are chosen to make positions readable. Nothing here feeds
//! back into the motion.

/// The share of its color a surface keeps where no light falls on it.
const AMBIENT: f64 = 0.25;
/// The direction toward the light, of length 1, as its parts to the
/// right, down and toward the eye: it stands to the upper left, in
/// front of the picture.
const LIGHT: [f64; 3] = [-1.0 / 3.0, -2.0 / 3.0, 2.0 / 3.0];

/// Packs red, green and blue, each 0 to 255, into one pixel.
pub fn rgb(color: [u8; 3]) -> u32 {
    (color[0] as u32) << 16 | (color[1] as u32) << 8 | color[2] as u32
}

/// A grid of `width` by `height` pixels.
pub struct Picture {
    pub width: usize,
    pub height: usize,
    /// The pixels, row by row from the top left.
    pub pixels: Vec<u32>,
    /// How near the nearest surface drawn at each pixel is, in pixels.
    depth: Vec<f64>,
}

impl Picture {
    /// A picture with every pixel `color` and nothing drawn in it.
    pub fn new(width: usize, height: usize, color: u32) -> Self {
        Self {
            width,
            height,
            pixels: vec![color; width * height],
            depth: vec![f64::NEG_INFINITY; width * height],
        }
    }

    /// Draws a ball of `radius` pixels with its center at (x, y), counted
    /// in pixels from the top left pixel, x to the right and y down, and
    /// at depth `z`, also in pixels, larger being nearer to the eye.
    /// `color` is its red, green and blue where it faces the light.
    /// Parts outside the picture are left out.
    pub fn ball(&mut self, x: f64, y: f64, z: f64, radius: f64, color: [u8; 3]) {
        let first = |c: f64| (c - radius).ceil().max(0.0) as usize;
        let end = |c: f64, size: usize| ((c + radius).floor() + 1.0).clamp(0.0, size as f64);
        for j in first(y)..end(y, self.height) as usize {
            for i in first(x)..end(x, self.width) as usize {
                let (dx, dy) = (i as f64 - x, j as f64 - y);
                let under = radius * radius - dx * dx - dy * dy;
                let k = j * self.width + i;
                if under <= 0.0 || z + radius <= self.depth[k] {
                    continue;
                }
                let dz = under.sqrt();
                if z + dz <= self.depth[k] {
                    continue;
                }
                self.depth[k] = z + dz;
                let facing = (dx * LIGHT[0] + dy * LIGHT[1] + dz * LIGHT[2]) / radius;
                let light = AMBIENT + (1.0 - AMBIENT) * facing.max(0.0);
                // Adding a half and dropping the fraction rounds to nearest.
                self.pixels[k] = rgb(color.map(|c| (c as f64 * light + 0.5) as u8));
            }
        }
    }

    /// Fills with `color` the rectangle `width` pixels wide and `height`
    /// high whose top left pixel is (x, y), in front of everything drawn
    /// before or after. Parts outside the picture are left out.
    pub fn block(&mut self, x: usize, y: usize, width: usize, height: usize, color: u32) {
        for j in y..(y + height).min(self.height) {
            for i in x..(x + width).min(self.width) {
                self.pixels[j * self.width + i] = color;
                self.depth[j * self.width + i] = f64::INFINITY;
            }
        }
    }
}