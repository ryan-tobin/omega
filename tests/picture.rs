//! Balls and blocks drawn into a picture (picture.rs), checked pixel by
//! pixel against counts and colors worked out by hand. All exact.

use omega::picture::{Picture, rgb};

const COLOR: [u8; 3] = [200, 100, 40];

#[test]
fn ball_covers_the_pixels_inside_its_circle() {
    // Radius 10 around pixel (20, 20). 317 whole-number points lie within
    // 10 of a point (OEIS A000328), 12 of them exactly at 10: (10, 0),
    // (6, 8) and (8, 6) with every sign and order. That leaves 305.
    let mut picture = Picture::new(41, 41, 0);
    picture.ball(20.0, 20.0, 0.0, 10.0, COLOR);
    assert_eq!(picture.pixels.iter().filter(|&&p| p != 0).count(), 305);
}

#[test]
fn ball_is_lit_from_the_upper_left() {
    // Radius 9. The middle faces the eye: n . l = 2/3, so the light there
    // is 1/4 + 3/4 * 2/3 = 3/4 and the color (150, 75, 30). The point 3
    // left and 6 up is 6 nearer than the center (9 + 36 + 36 = 81), so it
    // faces the light head on and has the full color. The point 3 right
    // and 6 down faces away, n . l = -1/9, and keeps a quarter.
    let mut picture = Picture::new(41, 41, 0);
    picture.ball(20.0, 20.0, 0.0, 9.0, COLOR);
    assert_eq!(picture.pixels[20 * 41 + 20], rgb([150, 75, 30]));
    assert_eq!(picture.pixels[14 * 41 + 17], rgb(COLOR));
    assert_eq!(picture.pixels[26 * 41 + 23], rgb([50, 25, 10]));
}

#[test]
fn nearer_ball_hides_the_one_behind_whichever_is_drawn_first() {
    // Two balls on one line of sight, the red one 30 pixels nearer. The
    // middle pixel is red at 3/4 light, in either order.
    for order in [[0, 1], [1, 0]] {
        let balls = [(30.0, [200, 0, 0]), (0.0, [0, 0, 200])];
        let mut picture = Picture::new(41, 41, 0);
        for k in order {
            picture.ball(20.0, 20.0, balls[k].0, 9.0, balls[k].1);
        }
        assert_eq!(picture.pixels[20 * 41 + 20], rgb([150, 0, 0]));
    }
}

#[test]
fn balls_that_run_into_each_other_are_cut_where_their_surfaces_cross() {
    // A red ball, and a blue one 6 pixels to its right and 3 nearer,
    // both of radius 10. At the red ball's middle the blue surface is at
    // 3 + 8 = 11, in front of red's 10. Three pixels to the left it is
    // at 3 + sqrt(19) = 7.4, behind red's sqrt(91) = 9.5. So the first
    // pixel is blue and the second red, in either order.
    for order in [[0, 1], [1, 0]] {
        let balls = [(20.0, 0.0, [200, 0, 0]), (26.0, 3.0, [0, 0, 200])];
        let mut picture = Picture::new(41, 41, 0);
        for k in order {
            picture.ball(balls[k].0, 20.0, balls[k].1, 10.0, balls[k].2);
        }
        let (middle, left) = (picture.pixels[20 * 41 + 20], picture.pixels[20 * 41 + 17]);
        assert!(middle >> 16 == 0 && middle & 0xff > 0, "{middle:x}");
        assert!(left >> 16 > 0 && left & 0xff == 0, "{left:x}");
    }
}

#[test]
fn drawing_past_the_edges_is_left_out() {
    // A ball around the top left pixel, one far outside, and a block that
    // runs off the bottom right: 2 x 2 of its pixels are inside, and a
    // ball drawn over them afterwards stays behind them.
    let mut picture = Picture::new(8, 8, 0);
    picture.ball(0.0, 0.0, 0.0, 5.0, COLOR);
    picture.ball(-50.0, 90.0, 0.0, 5.0, COLOR);
    assert_eq!(picture.pixels[0], rgb([150, 75, 30]));
    picture.block(6, 6, 5, 5, 7);
    picture.ball(7.0, 7.0, 0.0, 3.0, COLOR);
    assert_eq!(picture.pixels.iter().filter(|&&p| p == 7).count(), 4);
}