//! What the viewer does when the window, the keys, the mouse or the
//! clock ask for something. winit calls these functions (winit 0.30.13,
//! the ApplicationHandler trait in src/application.rs).

use crate::viewer::Viewer;
use softbuffer::Surface;
use std::f64::consts::FRAC_PI_2;
use std::rc::Rc;
use std::time::{Duration, Instant};
use winit::application::ApplicationHandler;
use winit::dpi::LogicalSize;
use winit::event::{ElementState, KeyEvent, MouseButton, StartCause, WindowEvent};
use winit::event_loop::{ActiveEventLoop, ControlFlow};
use winit::keyboard::{Key, NamedKey};
use winit::window::{Window, WindowId};

// Our own choices. A new frame every 33 ms: frames saved 0.05 tau apart
// then play 1.5 tau each second, which is 3.3 ps of argon (units.rs).
// Dragging 300 pixels, which is 150 points on a screen with 2 pixels to
// the point, turns the box a quarter turn. The window starts 1100 by
// 700 points.
const TICK: Duration = Duration::from_millis(33);
const TURN_PER_PIXEL: f64 = FRAC_PI_2 / 300.0;
const SIZE: LogicalSize<f64> = LogicalSize::new(1100.0, 700.0);

impl ApplicationHandler for Viewer {
    /// Called once the program may open its window.
    fn resumed(&mut self, event_loop: &ActiveEventLoop) {
        let attributes = Window::default_attributes()
            .with_title("omega: replay of a saved run")
            .with_inner_size(SIZE);
        let window = event_loop.create_window(attributes);
        let window = Rc::new(window.expect("cannot open a window"));
        let surface = Surface::new(&self.context, window);
        self.surface = Some(surface.expect("cannot draw into the window"));
        self.next_tick = Instant::now() + TICK;
        self.go(0);
    }

    /// Called each time the program wakes up. While playing, and once
    /// the window is open, moves one frame on when its time has come.
    fn new_events(&mut self, _: &ActiveEventLoop, _: StartCause) {
        let now = Instant::now();
        if self.playing && self.surface.is_some() && now >= self.next_tick {
            self.next_tick = now + TICK;
            self.go(1);
        }
    }

    fn window_event(&mut self, event_loop: &ActiveEventLoop, _: WindowId, event: WindowEvent) {
        match event {
            WindowEvent::CloseRequested => event_loop.exit(),
            WindowEvent::RedrawRequested => self.draw(),
            WindowEvent::MouseInput {
                state,
                button: MouseButton::Left,
                ..
            } => self.dragging = state == ElementState::Pressed,
            WindowEvent::CursorMoved { position, .. } => {
                if self.dragging {
                    self.camera.yaw += (position.x - self.cursor.x) * TURN_PER_PIXEL;
                    self.camera.pitch += (position.y - self.cursor.y) * TURN_PER_PIXEL;
                    self.go(0);
                }
                self.cursor = position;
            }
            WindowEvent::KeyboardInput { event, .. } => self.key(event_loop, event),
            _ => {}
        }
    }

    /// Called before the program sleeps: until the next frame is due
    /// while playing, or until something happens while paused.
    fn about_to_wait(&mut self, event_loop: &ActiveEventLoop) {
        event_loop.set_control_flow(match self.playing {
            true => ControlFlow::WaitUntil(self.next_tick),
            false => ControlFlow::Wait,
        });
    }
}

impl Viewer {
    /// Acts on a key going down: space, the arrows or escape. A held
    /// arrow keeps stepping. A held space bar counts once.
    fn key(&mut self, event_loop: &ActiveEventLoop, event: KeyEvent) {
        if event.state != ElementState::Pressed {
            return;
        }
        match event.logical_key {
            Key::Named(NamedKey::Escape) => event_loop.exit(),
            Key::Named(NamedKey::Space) if !event.repeat => {
                self.playing = !self.playing;
                self.go(0);
            }
            Key::Named(NamedKey::ArrowRight) => self.go(1),
            Key::Named(NamedKey::ArrowLeft) => self.go(self.replay.frames.len() - 1),
            _ => {}
        }
    }
}