//! Letters and digits for the numbers the viewer shows, drawn by us.
//!
//! A glyph is a grid of 5 by 7 cells, each inked or empty. Text is drawn
//! by filling a square of pixels for every inked cell. Each character
//! takes 6 cells across, 5 for its glyph and 1 empty, and a line of text
//! takes 9 cells down, 7 and 2 empty.
//!
//! Only the characters the viewer shows have a glyph: the digits, the
//! point, the minus sign, the slash and ten capital letters. Any other
//! character, the space among them, is drawn as an empty cell grid. No
//! two glyphs are alike, so text can be read back from a picture cell by
//! cell (tests/replay.rs does).

use crate::picture::Picture;

/// Cells across and down one glyph.
pub const WIDTH: usize = 5;
pub const HEIGHT: usize = 7;
/// Cells from one character to the next, and from one line to the next.
pub const ADVANCE: usize = 6;
pub const LINE: usize = 9;

/// Each glyph as its 7 rows, top to bottom, with # for an inked cell.
const GLYPHS: [(char, &str); 23] = [
    ('0', ".###. #...# #..## #.#.# ##..# #...# .###."),
    ('1', "..#.. .##.. ..#.. ..#.. ..#.. ..#.. .###."),
    ('2', ".###. #...# ....# ...#. ..#.. .#... #####"),
    ('3', ".###. #...# ....# ..##. ....# #...# .###."),
    ('4', "...#. ..##. .#.#. #..#. ##### ...#. ...#."),
    ('5', "##### #.... ####. ....# ....# #...# .###."),
    ('6', "..##. .#... #.... ####. #...# #...# .###."),
    ('7', "##### ....# ...#. ..#.. .#... .#... .#..."),
    ('8', ".###. #...# #...# .###. #...# #...# .###."),
    ('9', ".###. #...# #...# .#### ....# ...#. .##.."),
    ('.', "..... ..... ..... ..... ..... .##.. .##.."),
    ('-', "..... ..... ..... ##### ..... ..... ....."),
    ('/', "....# ....# ...#. ..#.. .#... #.... #...."),
    ('A', ".###. #...# #...# ##### #...# #...# #...#"),
    ('D', "####. #...# #...# #...# #...# #...# ####."),
    ('E', "##### #.... #.... ####. #.... #.... #####"),
    ('I', ".###. ..#.. ..#.. ..#.. ..#.. ..#.. .###."),
    ('M', "#...# ##.## #.#.# #.#.# #...# #...# #...#"),
    ('N', "#...# ##..# #.#.# #..## #...# #...# #...#"),
    ('P', "####. #...# #...# ####. #.... #.... #...."),
    ('S', ".#### #.... #.... .###. ....# ....# ####."),
    ('T', "##### ..#.. ..#.. ..#.. ..#.. ..#.. ..#.."),
    ('U', "#...# #...# #...# #...# #...# #...# .###."),
];

/// True if the glyph of `c` has ink in the cell `column` cells from its
/// left and `row` cells from its top. A character without a glyph has
/// ink nowhere.
pub fn ink(c: char, column: usize, row: usize) -> bool {
    let glyph = GLYPHS.iter().find(|(name, _)| *name == c);
    glyph.is_some_and(|(_, rows)| rows.as_bytes()[row * (WIDTH + 1) + column] == b'#')
}

/// The character whose glyph has exactly the inked cells of `cells`,
/// given row by row from the top, or a space if there is none.
pub fn read(cells: [[bool; WIDTH]; HEIGHT]) -> char {
    let same = |c: char| (0..HEIGHT).all(|row| (0..WIDTH).all(|x| ink(c, x, row) == cells[row][x]));
    let found = GLYPHS.iter().map(|(name, _)| *name).find(|&c| same(c));
    found.unwrap_or(' ')
}

/// Draws `text` in `color`, each cell a square of `size` pixels, with the
/// top left cell of its first character at pixel (x, y).
pub fn draw(picture: &mut Picture, text: &str, x: usize, y: usize, size: usize, color: u32) {
    for (k, c) in text.chars().enumerate() {
        for row in 0..HEIGHT {
            for column in (0..WIDTH).filter(|&column| ink(c, column, row)) {
                let left = x + (k * ADVANCE + column) * size;
                picture.block(left, y + row * size, size, size, color);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_glyph_is_7_rows_of_5_cells_and_no_two_are_alike() {
        for (k, (c, rows)) in GLYPHS.iter().enumerate() {
            let cells: Vec<&str> = rows.split(' ').collect();
            let cell = |b: u8| b == b'#' || b == b'.';
            let good = |r: &&str| r.len() == WIDTH && r.bytes().all(cell);
            assert!(cells.len() == HEIGHT && cells.iter().all(good), "{c}");
            let earlier = &GLYPHS[..k];
            assert!(earlier.iter().all(|(_, other)| other != rows), "{c}");
        }
    }

    #[test]
    fn a_one_is_drawn_cell_by_cell() {
        // "1" at size 2 with its corner at (3, 1). Its glyph has 10 inked
        // cells, so 40 pixels. The top row has ink in its third cell
        // only: pixels 7 and 8 of rows 1 and 2. Exact.
        let mut picture = Picture::new(20, 20, 0);
        draw(&mut picture, "1", 3, 1, 2, 9);
        assert_eq!(picture.pixels.iter().filter(|&&p| p == 9).count(), 40);
        let row: Vec<usize> = (0..20).filter(|i| picture.pixels[20 + i] == 9).collect();
        assert_eq!(row, [7, 8]);
        assert_eq!(picture.pixels[..20], [0; 20]);
        // A second character starts 6 cells after the first.
        let mut picture = Picture::new(20, 20, 0);
        draw(&mut picture, "11", 0, 0, 1, 9);
        let row: Vec<usize> = (0..20).filter(|&i| picture.pixels[i] == 9).collect();
        assert_eq!(row, [2, 8]);
    }

    #[test]
    fn glyphs_read_back_as_their_characters() {
        // Every glyph's cells give back its character, and cells that
        // match no glyph give a space.
        for (c, _) in GLYPHS {
            let mut cells = [[false; WIDTH]; HEIGHT];
            for (row, line) in cells.iter_mut().enumerate() {
                for (column, cell) in line.iter_mut().enumerate() {
                    *cell = ink(c, column, row);
                }
            }
            assert_eq!(read(cells), c);
        }
        assert_eq!(read([[true; WIDTH]; HEIGHT]), ' ');
        assert_eq!(read([[false; WIDTH]; HEIGHT]), ' ');
    }
}