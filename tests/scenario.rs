//! The scenario reader (scenario.rs): every key is read, optional keys
//! have their defaults, and a mistake is refused with its line.

use omega::scenario::{Scenario, Stage, parse};

const TEXT: &str = "# A test run.
cells 3
side 5.5   # sigma
fill 0.5
temperature 0.8
hot 0.25 3 0.9
seed 7
cutoff 2.5
dt 0.005
checkpoint 50

stage settle
steps 100

stage measure
steps 200
thermostat 0.8 0.5
barostat 1 0.05 25
sample 10
profile 20
";

#[test]
fn every_key_is_read() {
    let settle = Stage {
        name: "settle".into(),
        steps: 100,
        thermostat: None,
        barostat: None,
        sample: 0,
        profile: 0,
    };
    let measure = Stage {
        name: "measure".into(),
        steps: 200,
        thermostat: Some((0.8, 0.5)),
        barostat: Some((1.0, 0.05, 25)),
        sample: 10,
        profile: 20,
    };
    let want = Scenario {
        cells: 3,
        side: 5.5,
        fill: 0.5,
        temperature: 0.8,
        hot: Some((0.25, 3.0, 0.9)),
        seed: 7,
        cutoff: 2.5,
        dt: 0.005,
        checkpoint: 50,
        stages: vec![settle, measure],
    };
    assert_eq!(parse(TEXT), Ok(want));
}

#[test]
fn optional_keys_have_defaults() {
    let text = "cells 3\nside 5.5\ntemperature 1\nseed 0\ncutoff 2\ndt 0.01\nstage a\nsteps 5";
    let s = parse(text).unwrap();
    assert_eq!((s.fill, s.hot, s.checkpoint), (1.0, None, 0));
}

#[test]
fn mistakes_are_refused_with_their_line() {
    let swap = |from: &str, to: &str| parse(&TEXT.replace(from, to));
    let line = |n: usize, text: &str| Err(format!("line {n}: cannot use `{text}`"));
    assert_eq!(swap("fill 0.5", "full 0.5"), line(4, "full 0.5"));
    assert_eq!(swap("side 5.5 ", "side 5.5 6"), line(3, "side 5.5 6"));
    assert_eq!(swap("cells 3", "cells 3.5"), line(2, "cells 3.5"));
    assert_eq!(swap("seed 7", "seed x"), line(7, "seed x"));
    assert_eq!(swap("sample 10", "cutoff 3"), line(19, "cutoff 3"));
    assert_eq!(swap("checkpoint 50", "steps 50"), line(10, "steps 50"));
    assert!(swap("seed 7\n", "").unwrap_err().contains("seed"));
    assert!(swap("steps 100", "").unwrap_err().contains("stage settle"));
    assert!(
        swap("thermostat 0.8 0.5\n", "")
            .unwrap_err()
            .contains("barostat")
    );
    assert!(swap("fill 0.5", "fill 1.5").unwrap_err().contains("fill"));
    assert!(
        swap("hot 0.25 3 0.9", "hot 0.25 3 0")
            .unwrap_err()
            .contains("share")
    );
    assert_eq!(swap("hot 0.25 3 0.9", "hot 0.25 3"), line(6, "hot 0.25 3"));
    assert!(swap("barostat 1 0.05 25", "barostat 1 0.05 0").is_err());
    assert!(parse("cells 3\nside 5.5\ntemperature 1\nseed 0\ncutoff 2\ndt 0.01").is_err());
}
