//! One axis of a plot: a range of data values laid along one side of the
//! frame, with tick marks at even steps.

/// One axis: its label, the range of data values it spans, the spacing of
/// its tick marks, in the unit of the data, and how many decimals the tick
/// labels show.
pub struct Axis<'a> {
    pub label: &'a str,
    pub min: f64,
    pub max: f64,
    pub step: f64,
    pub decimals: usize,
}

impl<'a> Axis<'a> {
    /// The axis from `range.0` to `range.1`, ticked every `step`.
    pub fn new(label: &'a str, range: (f64, f64), step: f64, decimals: usize) -> Self {
        let (min, max) = range;
        Self {
            label,
            min,
            max,
            step,
            decimals,
        }
    }

    /// How far along the axis `value` lies: 0 at min, 1 at max, and
    /// outside 0 to 1 beyond them. A plain number.
    pub fn fraction(&self, value: f64) -> f64 {
        (value - self.min) / (self.max - self.min)
    }

    /// The tick values min, min + step, and so on up to max.
    pub fn ticks(&self) -> Vec<f64> {
        let count = ((self.max - self.min) / self.step).round() as usize;
        (0..=count)
            .map(|i| self.min + i as f64 * self.step)
            .collect()
    }

    /// The label of the tick at `value`, with the axis's decimals.
    pub fn tick_label(&self, value: f64) -> String {
        format!("{value:.*}", self.decimals)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn values_map_to_fractions_of_the_axis() {
        // Every value is exact in an f64.
        let axis = Axis::new("g", (-1.0, 3.0), 0.5, 1);
        assert_eq!(axis.fraction(-1.0), 0.0);
        assert_eq!(axis.fraction(3.0), 1.0);
        assert_eq!(axis.fraction(0.0), 0.25);
        assert_eq!(axis.fraction(5.0), 1.5);
    }

    #[test]
    fn ticks_run_from_min_to_max_with_their_labels() {
        // 0.1 has no exact f64, so 3 * 0.1 is 0.30000000000000004, and the
        // label rounds it to the axis's one decimal.
        let axis = Axis::new("t", (0.0, 0.5), 0.1, 1);
        let ticks = axis.ticks();
        assert_eq!(ticks.len(), 6);
        assert_eq!((ticks[0], ticks[5]), (0.0, 0.5));
        assert_eq!(ticks[3], 0.30000000000000004);
        assert_eq!(axis.tick_label(ticks[3]), "0.3");
    }
}