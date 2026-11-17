//! What the viewer holds and how it draws a frame into its window.
//!
//! The window code follows softbuffer's winit example (softbuffer 0.4.8,
//! examples/winit.rs): set the buffer to the window's size in pixels,
//! fill it, present it.

use omega::hud;
use omega::picture::Picture;
use omega::replay::Replay;
use omega::scene::{self, Camera};
use softbuffer::{Context, Surface};
use std::num::NonZeroU32;
use std::rc::Rc;
use std::time::Instant;
use winit::dpi::PhysicalPosition;
use winit::event_loop::OwnedDisplayHandle;
use winit::window::Window;

// Our own choices: a dark background, text cells 2 points wide (4 pixels
// on a screen with 2 pixels to the point), and the box first seen from
// the left and from above, turned by 0.6 and tipped by 0.4 radians.
const BACKGROUND: u32 = 0x0010_1418;
const CELL_POINTS: f64 = 2.0;
const START: Camera = Camera {
    yaw: 0.6,
    pitch: 0.4,
};

/// What the viewer keeps: the run, what is on screen, and the window.
pub struct Viewer {
    pub replay: Replay,
    /// Each frame measured again from its atoms: the numbers shown.
    pub shown: Vec<[f64; 5]>,
    /// The speed drawn red, in sigma / tau: twice the root-mean-square
    /// speed of the first frame.
    pub full: f64,
    /// The frame on screen.
    pub frame: usize,
    pub playing: bool,
    /// When the next frame is due while playing.
    pub next_tick: Instant,
    pub camera: Camera,
    /// True while the left mouse button is held down.
    pub dragging: bool,
    /// Where the mouse was last seen, in pixels.
    pub cursor: PhysicalPosition<f64>,
    pub context: Context<OwnedDisplayHandle>,
    /// The window's pixels, once the window exists.
    pub surface: Option<Surface<OwnedDisplayHandle, Rc<Window>>>,
}

impl Viewer {
    /// A viewer of `replay` that starts playing at its first frame.
    pub fn new(replay: Replay, shown: Vec<[f64; 5]>, context: Context<OwnedDisplayHandle>) -> Self {
        Self {
            full: scene::full_speed(&replay.frames[0].velocities),
            replay,
            shown,
            frame: 0,
            playing: true,
            next_tick: Instant::now(),
            camera: START,
            dragging: false,
            cursor: PhysicalPosition::new(0.0, 0.0),
            context,
            surface: None,
        }
    }

    /// Moves `by` frames on, around the end of the run, and asks the
    /// window to draw again.
    pub fn go(&mut self, by: usize) {
        self.frame = (self.frame + by) % self.replay.frames.len();
        if let Some(surface) = &self.surface {
            surface.window().request_redraw();
        }
    }

    /// Draws the current frame and its numbers, and shows them.
    pub fn draw(&mut self) {
        let Some(surface) = &mut self.surface else {
            return;
        };
        let size = surface.window().inner_size();
        let (Some(w), Some(h)) = (NonZeroU32::new(size.width), NonZeroU32::new(size.height)) else {
            return;
        };
        // The width of a text cell in whole pixels.
        let points = surface.window().scale_factor();
        let cell = (CELL_POINTS * points).round().max(1.0) as usize;
        let mut picture = Picture::new(size.width as usize, size.height as usize, BACKGROUND);
        let (frame, space) = (&self.replay.frames[self.frame], self.replay.space);
        scene::draw(&mut picture, frame, space, &self.camera, self.full);
        let lines = hud::lines(self.shown[self.frame], !self.playing);
        hud::draw(&mut picture, &lines, cell);
        surface.resize(w, h).expect("cannot size the pixels");
        let mut buffer = surface.buffer_mut().expect("cannot get the pixels");
        buffer.copy_from_slice(&picture.pixels);
        buffer.present().expect("cannot show the picture");
    }
}