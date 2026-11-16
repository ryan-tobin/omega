//! SVG, the W3C text format for drawings, which web browsers and GitHub
//! show as images (Scalable Vector Graphics 1.1, Second Edition, 2011).
//!
//! A document is one svg element holding shapes and text. Each function
//! here returns one element as text. Places are in pixels from the top
//! left corner, x to the right and y downward, written with one decimal.
//! The characters & and < cannot appear in text as themselves, so they
//! are written as &amp; and &lt;.

/// A place or a size in pixels, (x, y) or (width, height).
pub type Point = (f64, f64);

/// The name of the outermost element, which holds the whole drawing.
pub const ROOT: &str = "svg";

/// The opening of a document `size` pixels big, with a white background
/// so it reads the same on a dark page. Text is sans-serif, 13 pixels
/// high. The document ends with `end()`.
pub fn begin(size: Point) -> String {
    let (w, h) = size;
    let open = format!(r#"<{ROOT} xmlns="http://www.w3.org/2000/svg""#);
    let fonts = r#"font-family="sans-serif" font-size="13""#;
    let view = format!(r#"width="{w:.1}" height="{h:.1}" viewBox="0 0 {w:.1} {h:.1}""#);
    let white = rect((0.0, 0.0), size, r#"fill="white""#);
    format!("{open} {view} {fonts}>\n{white}")
}

/// The end of a document, closing what `begin` opened.
pub fn end() -> String {
    format!("</{ROOT}>\n")
}

/// A straight line from `a` to `b`. `style` holds presentation attributes,
/// such as `stroke="black"`, here and below.
pub fn line(a: Point, b: Point, style: &str) -> String {
    let ((x1, y1), (x2, y2)) = (a, b);
    format!(r#"<line x1="{x1:.1}" y1="{y1:.1}" x2="{x2:.1}" y2="{y2:.1}" {style}/>"#)
}

/// Straight lines joining `points` in order.
pub fn polyline(points: &[Point], style: &str) -> String {
    let pairs: Vec<String> = points
        .iter()
        .map(|(x, y)| format!("{x:.1},{y:.1}"))
        .collect();
    format!(
        r#"<polyline points="{}" fill="none" {style}/>"#,
        pairs.join(" ")
    )
}

/// A dot of radius `r` at `center`, filled with `color`.
pub fn dot(center: Point, r: f64, color: &str) -> String {
    let (x, y) = center;
    format!(r#"<circle cx="{x:.1}" cy="{y:.1}" r="{r:.1}" fill="{color}"/>"#)
}

/// A rectangle with its top left corner at `corner`.
pub fn rect(corner: Point, size: Point, style: &str) -> String {
    let ((x, y), (w, h)) = (corner, size);
    format!(r#"<rect x="{x:.1}" y="{y:.1}" width="{w:.1}" height="{h:.1}" {style}/>"#)
}

/// Text whose baseline passes through `at`. `anchor` is start, middle or
/// end: the part of the text that sits at `at`. `degrees` turns it about
/// `at`, clockwise, so -90 reads upward.
pub fn text(at: Point, anchor: &str, degrees: f64, words: &str) -> String {
    let (x, y) = at;
    let words = words.replace('&', "&amp;").replace('<', "&lt;");
    let turn = format!(r#"transform="rotate({degrees:.1} {x:.1} {y:.1})""#);
    let place = format!(r#"x="{x:.1}" y="{y:.1}" text-anchor="{anchor}""#);
    if degrees == 0.0 {
        format!("<text {place}>{words}</text>")
    } else {
        format!("<text {place} {turn}>{words}</text>")
    }
}

/// The opening of a group whose contents are cut off outside the
/// rectangle at `corner` of `size`. The group ends with `</g>`.
pub fn clipped_group(corner: Point, size: Point) -> String {
    let area = rect(corner, size, "");
    format!("<clipPath id=\"area\">{area}</clipPath>\n<g clip-path=\"url(#area)\">")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn document_opens_with_the_svg_namespace_on_white() {
        let want = format!(
            "<{ROOT} xmlns=\"http://www.w3.org/2000/svg\" width=\"640.0\" \
            height=\"400.0\" viewBox=\"0 0 640.0 400.0\" font-family=\"sans-serif\" \
            font-size=\"13\">\n<rect x=\"0.0\" y=\"0.0\" width=\"640.0\" \
            height=\"400.0\" fill=\"white\"/>"
        );
        assert_eq!(begin((640.0, 400.0)), want);
        assert_eq!(end(), "</svg>\n");
    }

    #[test]
    fn shapes_round_to_a_tenth_of_a_pixel() {
        let want = "<polyline points=\"1.0,2.2 3.1,4.0\" fill=\"none\" stroke=\"red\"/>";
        assert_eq!(
            polyline(&[(1.0, 2.24), (3.06, 4.0)], "stroke=\"red\""),
            want
        );
        let want = "<circle cx=\"5.0\" cy=\"6.0\" r=\"3.0\" fill=\"blue\"/>";
        assert_eq!(dot((5.0, 6.0), 3.0, "blue"), want);
    }

    #[test]
    fn text_escapes_ampersand_and_less_than() {
        let want = "<text x=\"10.0\" y=\"20.0\" text-anchor=\"end\" \
            transform=\"rotate(-90.0 10.0 20.0)\">r &lt; L/2 &amp; g</text>";
        assert_eq!(text((10.0, 20.0), "end", -90.0, "r < L/2 & g"), want);
        let want = "<text x=\"1.0\" y=\"2.0\" text-anchor=\"start\">g(r)</text>";
        assert_eq!(text((1.0, 2.0), "start", 0.0, "g(r)"), want);
    }
}