//! The numbers the viewer shows over the atoms, as lines of text.
//!
//! The lines are drawn from the top left of the picture with the glyphs
//! of font.rs, one line under the other. `read` gets a line back from the
//! pixels, which is how tests/replay.rs checks that the numbers on screen
//! are the numbers of the run.

use crate::font::{self, ADVANCE, HEIGHT, LINE, WIDTH};
use crate::picture::Picture;

/// The color of the text: white.
pub const INK: u32 = 0x00ff_ffff;
/// The empty edge left of the text and above it, in cells.
const MARGIN: usize = 3;

/// The lines of text for one frame. `numbers` are a row of thermo.csv
/// (recorder.rs): the time in tau, the temperature in eps / kB, the
/// potential energy per atom in eps, the pressure in eps / sigma^3 and
/// the total energy in eps. The time is shown to 2 decimals, and the
/// rest to 4 as the measure program prints them. A last line says PAUSED
/// while the replay stands still.
pub fn lines(numbers: [f64; 5], paused: bool) -> Vec<String> {
    let [time, t, u, p, e] = numbers;
    let mut lines = vec![
        format!("TIME {time:10.2}"),
        format!("T    {t:10.4}"),
        format!("U/N  {u:10.4}"),
        format!("P    {p:10.4}"),
        format!("E    {e:10.4}"),
    ];
    if paused {
        lines.push(String::from("PAUSED"));
    }
    lines
}

/// Draws `lines` from the top left of `picture`, one under the other,
/// each cell of a glyph a square of `size` pixels.
pub fn draw(picture: &mut Picture, lines: &[String], size: usize) {
    for (k, line) in lines.iter().enumerate() {
        let y = (MARGIN + k * LINE) * size;
        font::draw(picture, line, MARGIN * size, y, size, INK);
    }
}

/// Reads back the first `count` characters of line `line`, counted from
/// 0, from a picture that `draw` drew into with the same `size`. A cell
/// counts as inked if the pixel at its middle has the color of the text.
pub fn read(picture: &Picture, line: usize, count: usize, size: usize) -> String {
    let inked = |column: usize, row: usize| {
        let (x, y) = (MARGIN + column, MARGIN + line * LINE + row);
        let pixel = (y * size + size / 2) * picture.width + x * size + size / 2;
        picture.pixels[pixel] == INK
    };
    let glyph = |k: usize| {
        let mut cells = [[false; WIDTH]; HEIGHT];
        for (row, cells) in cells.iter_mut().enumerate() {
            for (column, cell) in cells.iter_mut().enumerate() {
                *cell = inked(k * ADVANCE + column, row);
            }
        }
        font::read(cells)
    };
    (0..count).map(glyph).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    // The measure program's row at 5 tau (STATUS, M3), with more digits.
    const ROW: [f64; 5] = [5.0, 0.7918501, -4.6003391, 0.9249162, -2949.6429469];

    #[test]
    fn lines_show_the_numbers_as_the_measure_program_prints_them() {
        let want = [
            "TIME       5.00",
            "T        0.7919",
            "U/N     -4.6003",
            "P        0.9249",
            "E    -2949.6429",
        ];
        assert_eq!(lines(ROW, false), want);
        assert_eq!(lines(ROW, true)[..5], want);
        assert_eq!(lines(ROW, true)[5], "PAUSED");
    }

    #[test]
    fn text_starts_3_cells_in_and_lines_are_9_cells_apart() {
        // At size 1 the T of TIME has its top left cell at pixel (3, 3),
        // and the T of the line below is 9 pixels under it. Exact.
        let mut picture = Picture::new(110, 60, 0);
        draw(&mut picture, &lines(ROW, false), 1);
        let inked = |x: usize, y: usize| picture.pixels[y * 110 + x] == INK;
        assert!(inked(3, 3) && !inked(2, 3) && !inked(3, 2));
        assert!(inked(3, 12) && !inked(3, 11) && !inked(3, 10));
    }

    #[test]
    fn text_drawn_is_read_back_from_the_pixels() {
        // Every character of every line, at three sizes. Exact. A
        // character without a glyph would come back as a space.
        let lines = lines(ROW, true);
        for size in [1, 2, 5] {
            let mut picture = Picture::new(110 * size, 60 * size, 0x0010_1418);
            draw(&mut picture, &lines, size);
            for (k, line) in lines.iter().enumerate() {
                assert_eq!(&read(&picture, k, line.len(), size), line, "size {size}");
            }
        }
    }
}