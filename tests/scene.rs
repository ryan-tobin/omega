//! The camera, the colors and the drawing of atoms (scene.rs), checked
//! against numbers known in advance. Exact unless a tolerance is given,
//! and each tolerance is our own choice and says why.

use omega::periodic::PeriodicBox;
use omega::picture::{Picture, rgb};
use omega::replay::Snapshot;
use omega::rng::Rng;
use omega::scene::{self, Camera};
use omega::vec3::Vec3;
use omega::velocities;
use std::f64::consts::FRAC_PI_2;

const ZERO: Vec3 = Vec3::new(0.0, 0.0, 0.0);

/// One atom at rest at `position` in a box of side 10, drawn into a
/// picture of 400 by 300 as seen with the given `yaw` and `pitch`.
fn one_atom(position: Vec3, yaw: f64, pitch: f64) -> Picture {
    let frame = Snapshot {
        positions: vec![position],
        velocities: vec![ZERO],
    };
    let space = PeriodicBox { side: 10.0 };
    let mut picture = Picture::new(400, 300, 0);
    scene::draw(&mut picture, &frame, space, &Camera { yaw, pitch }, 1.0);
    picture
}

#[test]
fn a_turn_sends_the_axes_where_it_should_and_keeps_lengths() {
    // A quarter turn about the upright axis points x away from the eye
    // and z to the right. A quarter tip points y at the eye and leaves x
    // alone. 1e-15 allows for the rounding of pi / 2, our own choice. A
    // vector of length 13 keeps its length to 1e-14, 5.6 units in the last
    // place of 13 (measured: 1), our own choice.
    let close = |a: Vec3, b: Vec3| (a - b).norm() <= 1e-15;
    let x = Vec3::new(1.0, 0.0, 0.0);
    let y = Vec3::new(0.0, 1.0, 0.0);
    let z = Vec3::new(0.0, 0.0, 1.0);
    let (yaw, pitch) = (FRAC_PI_2, 0.0);
    let turned = Camera { yaw, pitch };
    assert!(close(turned.turn(x), -z) && close(turned.turn(z), x));
    let (yaw, pitch) = (0.0, FRAC_PI_2);
    let tipped = Camera { yaw, pitch };
    assert!(close(tipped.turn(y), z) && close(tipped.turn(x), x));
    let (yaw, pitch) = (0.6, 0.4);
    let length = Camera { yaw, pitch }.turn(Vec3::new(3.0, 4.0, 12.0)).norm();
    assert!((length - 13.0).abs() <= 1e-14, "{length}");
}

#[test]
fn color_runs_from_blue_through_white_to_red() {
    // With red at speed 4: blue at rest, white at 2, halfway between
    // them at 1 and at 3, and red from 4 on.
    let color = |speed: f64| scene::color_of_speed(speed, 4.0);
    assert_eq!(color(0.0), [32, 96, 240]);
    assert_eq!(color(1.0), [136, 168, 240]);
    assert_eq!(color(2.0), [240, 240, 240]);
    assert_eq!(color(3.0), [240, 152, 136]);
    assert_eq!(color(4.0), [240, 64, 32]);
    assert_eq!(color(9.0), [240, 64, 32]);
    // Atoms all at rest have no fastest speed. They are still blue.
    assert_eq!(scene::color_of_speed(0.0, 0.0), [32, 96, 240]);
    // Speeds 5 and 5: root mean square 5, so red at 10.
    let v = [Vec3::new(3.0, 4.0, 0.0), Vec3::new(0.0, 0.0, 5.0)];
    assert_eq!(scene::full_speed(&v), 10.0);
}

#[test]
fn thirty_nine_percent_of_atoms_are_redder_than_white() {
    // 100,000 velocities drawn for T = 1. An atom is redder than white,
    // with less blue than 240, if it is faster than the root-mean-square
    // speed: 39.16% of the Maxwell-Boltzmann distribution (scene.rs).
    // Standard error sqrt(0.39 * 0.61 / 100,000) = 0.0015. Tolerance
    // 0.008, about 5 standard errors, our own choice.
    let v = velocities::maxwell_boltzmann(100_000, 1.0, &mut Rng::new(11));
    let full = scene::full_speed(&v);
    let warm = |v: &&Vec3| scene::color_of_speed(v.norm(), full)[2] < 240;
    let share = v.iter().filter(warm).count() as f64 / v.len() as f64;
    assert!((share - 0.3916).abs() <= 0.008, "share {share}");
}

#[test]
fn an_atom_in_the_middle_of_the_box_is_drawn_in_the_middle_of_the_picture() {
    // The middle of the box, given 3 boxes to the right and 2 down, seen
    // from an odd angle. The scale is 0.9 * 150 / (sqrt(75) + 0.5) =
    // 14.74 pixels per sigma, so the ball's radius is 7.37 and it covers
    // the 177 pixels with dx^2 + dy^2 <= 53 (OEIS A000328 gives 149
    // within 7; 50, 52 and 53 add 12, 8 and 8). Its middle is blue,
    // darkened to 3/4 for being halfway back, and lit at 3/4.
    let picture = one_atom(Vec3::new(35.0, -15.0, 5.0), 0.6, 0.4);
    assert_eq!(picture.pixels.iter().filter(|&&p| p != 0).count(), 177);
    assert_eq!(picture.pixels[150 * 400 + 200], rgb([18, 54, 135]));
}

#[test]
fn a_near_atom_is_drawn_brighter_than_a_far_one() {
    // An atom 3 sigma in front of the middle, seen straight on, and one
    // 3 behind. The near one is given a whole box nearer still, which
    // wrapping takes back. The corners are sqrt(75) = 8.66 from the
    // middle, so the near one keeps 1/2 + 1/4 (1 + 3 / 8.66) = 0.8366 of
    // its color and the far one 0.6634. Blue at rest is 240: 201 and 159,
    // and at 3/4 light 151 and 119.
    let blue = |z: f64| one_atom(Vec3::new(5.0, 5.0, z), 0.0, 0.0).pixels[150 * 400 + 200] & 0xff;
    assert_eq!((blue(18.0), blue(2.0)), (151, 119));
}

#[test]
fn atoms_keep_right_up_and_near_and_show_their_own_speed() {
    // Seen straight on, at 14.74 pixels per sigma. An atom 2 sigma right
    // of the middle is drawn 29 pixels right of it, and one 2 sigma up
    // is drawn 29 pixels up, with nothing opposite either. A red atom at
    // full speed sits half a sigma right of and in front of a blue one
    // at rest in the middle. 3 pixels right of the middle the red
    // surface is 13.3 pixels toward the eye and the blue one 6.7, so red
    // shows. 3 pixels left there is only blue.
    let offsets = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0; 3], [0.5, 0.0, 0.5]];
    let at = |o: [f64; 3]| Vec3::new(5.0 + o[0], 5.0 + o[1], 5.0 + o[2]);
    let frame = Snapshot {
        positions: offsets.map(at).to_vec(),
        velocities: vec![ZERO, ZERO, ZERO, Vec3::new(0.0, 3.0, 0.0)],
    };
    let (space, yaw, pitch) = (PeriodicBox { side: 10.0 }, 0.0, 0.0);
    let mut picture = Picture::new(400, 300, 0);
    scene::draw(&mut picture, &frame, space, &Camera { yaw, pitch }, 3.0);
    let pixel = |x: usize, y: usize| picture.pixels[y * 400 + x];
    assert!(pixel(229, 150) != 0 && pixel(171, 150) == 0);
    assert!(pixel(200, 121) != 0 && pixel(200, 179) == 0);
    assert!(pixel(203, 150) >> 16 > pixel(203, 150) & 0xff);
    assert!(pixel(197, 150) >> 16 < pixel(197, 150) & 0xff);
}