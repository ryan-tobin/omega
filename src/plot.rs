//! Line plots of measured curves, as SVG documents (svg.rs).
//!
//! A plot is 640 by 400 pixels: a title, a frame with grid lines at the
//! tick marks, axis labels, and series of points drawn as solid lines,
//! dashed lines or dots, named in a legend at the top right of the frame.
//! The frame's top left corner is at (LEFT, TOP), and a data point (x, y)
//! lands (x - x_min) / (x_max - x_min) of the frame's width to the right of
//! it and (y_max - y) / (y_max - y_min) of its height below it: SVG counts
//! y downward. Anything outside the frame is cut off.

use crate::svg::{self, Point};

mod axis;
pub use axis::Axis;

/// How a series is drawn.
pub enum Mark {
    Line,
    Dashed,
    Dots,
}

/// Points (x, y), in the units of the axes, drawn in one color (an SVG
/// color name such as "black") and named in the legend.
pub struct Series<'a> {
    pub name: &'a str,
    pub color: &'a str,
    pub mark: Mark,
    pub points: &'a [Point],
}

impl<'a> Series<'a> {
    /// The series with these settings.
    pub fn new(name: &'a str, color: &'a str, mark: Mark, points: &'a [Point]) -> Self {
        Self {
            name,
            color,
            mark,
            points,
        }
    }
}

// The whole plot, and the corner and size of the frame, in pixels.
const SIZE: Point = (640.0, 400.0);
const LEFT: f64 = 70.0;
const TOP: f64 = 40.0;
const FRAME: Point = (550.0, 310.0);

/// Where the data point `p` lands, in pixels.
fn pixel(x: &Axis, y: &Axis, p: Point) -> Point {
    let across = x.fraction(p.0) * FRAME.0;
    let down = (1.0 - y.fraction(p.1)) * FRAME.1;
    (LEFT + across, TOP + down)
}

/// The pixel points `at` drawn the way `s` is drawn. Lines are 1.5 pixels
/// wide, dashes 6 pixels on and 4 off, and dots 3 pixels in radius.
fn draw(s: &Series, at: &[Point]) -> String {
    let stroke = format!("stroke=\"{}\" stroke-width=\"1.5\"", s.color);
    match s.mark {
        Mark::Line => svg::polyline(at, &stroke),
        Mark::Dashed => svg::polyline(at, &format!("{stroke} stroke-dasharray=\"6 4\"")),
        Mark::Dots => at.iter().map(|&p| svg::dot(p, 3.0, s.color)).collect(),
    }
}

/// A whole SVG document: `series` plotted against the axes `x` and `y`,
/// under `title`.
pub fn svg(title: &str, x: &Axis, y: &Axis, series: &[Series]) -> String {
    let (right, bottom) = (LEFT + FRAME.0, TOP + FRAME.1);
    let grid = "stroke=\"#dddddd\"";
    let mut parts = vec![svg::begin(SIZE)];
    parts.push(svg::text((LEFT, 24.0), "start", 0.0, title));
    for value in x.ticks() {
        let (px, _) = pixel(x, y, (value, y.min));
        parts.push(svg::line((px, TOP), (px, bottom), grid));
        let label = x.tick_label(value);
        parts.push(svg::text((px, bottom + 18.0), "middle", 0.0, &label));
    }
    for value in y.ticks() {
        let (_, py) = pixel(x, y, (x.min, value));
        parts.push(svg::line((LEFT, py), (right, py), grid));
        let label = y.tick_label(value);
        parts.push(svg::text((LEFT - 6.0, py + 4.0), "end", 0.0, &label));
    }
    let edge = "fill=\"none\" stroke=\"black\"";
    parts.push(svg::rect((LEFT, TOP), FRAME, edge));
    let middle = (LEFT + FRAME.0 / 2.0, TOP + FRAME.1 / 2.0);
    parts.push(svg::text((middle.0, SIZE.1 - 10.0), "middle", 0.0, x.label));
    parts.push(svg::text((20.0, middle.1), "middle", -90.0, y.label));
    parts.push(svg::clipped_group((LEFT, TOP), FRAME));
    for s in series {
        let at: Vec<Point> = s.points.iter().map(|&p| pixel(x, y, p)).collect();
        parts.push(draw(s, &at));
    }
    parts.push("</g>".to_string());
    for (i, s) in series.iter().enumerate() {
        // The legend: one row per series, a short sample of the series,
        // a line or one dot, with its name to the left.
        let row = TOP + 20.0 + 20.0 * i as f64;
        parts.push(svg::text((right - 46.0, row + 4.0), "end", 0.0, s.name));
        let sample = match s.mark {
            Mark::Dots => vec![(right - 26.0, row)],
            _ => vec![(right - 40.0, row), (right - 12.0, row)],
        };
        parts.push(draw(s, &sample));
    }
    parts.push(svg::end());
    parts.join("\n")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn axes() -> (Axis<'static>, Axis<'static>) {
        (
            Axis::new("r", (0.0, 5.0), 1.0, 0),
            Axis::new("g", (-1.0, 3.0), 0.5, 1),
        )
    }

    #[test]
    fn corners_of_the_data_land_on_the_corners_of_the_frame() {
        // The frame runs from (70, 40) to (620, 350). Every value is exact.
        let (x, y) = axes();
        assert_eq!(pixel(&x, &y, (0.0, -1.0)), (70.0, 350.0));
        assert_eq!(pixel(&x, &y, (5.0, 3.0)), (620.0, 40.0));
        assert_eq!(pixel(&x, &y, (2.5, 1.0)), (345.0, 195.0));
    }

    #[test]
    fn a_series_is_drawn_at_its_pixels_and_named() {
        let (x, y) = axes();
        let points = [(0.0, -1.0), (5.0, 3.0)];
        let line = Series::new("ours", "black", Mark::Line, &points);
        let dashes = Series::new("his", "red", Mark::Dashed, &points);
        let text = svg("test", &x, &y, &[line, dashes]);
        let want = "<polyline points=\"70.0,350.0 620.0,40.0\" fill=\"none\" \
            stroke=\"black\" stroke-width=\"1.5\"/>";
        assert!(text.contains(want), "{text}");
        let dash = "stroke=\"red\" stroke-width=\"1.5\" stroke-dasharray=\"6 4\"/>";
        assert!(text.contains(dash), "{text}");
        assert!(text.contains(">ours</text>"), "{text}");
        assert!(text.starts_with(&svg::begin(SIZE)) && text.ends_with(&svg::end()));
    }
}
