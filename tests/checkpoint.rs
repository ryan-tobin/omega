//! The checkpoint format (checkpoint.rs): bytes read back to the same
//! checkpoint, the layout is the one the file describes, and damaged bytes
//! are refused.

use omega::checkpoint::Checkpoint;
use omega::vec3::Vec3;

fn small(list: Option<(f64, Vec<Vec3>)>) -> Checkpoint {
    Checkpoint {
        step: 7,
        side: 8.0,
        cutoff: 2.5,
        thermostat: (-0.125, 0.1),
        barostat: (3, 2),
        positions: vec![Vec3::new(1.0, 2.0, 3.0), Vec3::new(0.1, -0.0, 1e-300)],
        velocities: vec![
            Vec3::new(-1.5, 0.5, 0.25),
            Vec3::new(f64::MAX, 0.3, 1.0 / 3.0),
        ],
        list,
    }
}

#[test]
fn bytes_read_back_to_the_same_checkpoint() {
    // Every f64 must come back as the same bits, 0.1, 1/3, -0.0 and
    // the largest f64 included. With and without a neighbor list.
    let built_at = vec![Vec3::new(0.9, 2.0, 3.0), Vec3::new(0.2, 0.0, 0.0)];
    for list in [None, Some((0.3, built_at))] {
        let checkpoint = small(list);
        let back = Checkpoint::from_bytes(&checkpoint.to_bytes()).unwrap();
        assert_eq!(back, checkpoint);
        assert!(back.positions[1].y.is_sign_negative(), "-0.0 kept its sign");
    }
}

#[test]
fn layout_matches_the_one_written_at_the_top() {
    // 2 atoms without a list: 12 + 9 * 8 + 2 * 2 * 24 = 180 bytes.
    // With one: 228. The step is at byte 12, N at 20, the side at 28,
    // the skin at 76, and the first position at 84.
    let bytes = small(None).to_bytes();
    assert_eq!(bytes.len(), 180);
    assert_eq!(&bytes[..12], b"omega-cp\x01\0\0\0");
    assert_eq!(bytes[12..20], 7u64.to_le_bytes());
    assert_eq!(bytes[20..28], 2u64.to_le_bytes());
    assert_eq!(bytes[28..36], 8.0f64.to_le_bytes());
    assert_eq!(bytes[76..84], (-1.0f64).to_le_bytes());
    assert_eq!(bytes[84..92], 1.0f64.to_le_bytes());
    let with_list = small(Some((0.3, vec![Vec3::new(0.0, 0.0, 0.0); 2])));
    assert_eq!(with_list.to_bytes().len(), 228);
}

#[test]
fn damaged_bytes_are_refused() {
    let good = small(None).to_bytes();
    assert!(Checkpoint::from_bytes(&good).is_ok());
    assert!(Checkpoint::from_bytes(&good[..179]).is_err(), "cut short");
    assert!(Checkpoint::from_bytes(&good[..40]).is_err(), "header only");
    let mut long = good.clone();
    long.push(0);
    assert!(Checkpoint::from_bytes(&long).is_err(), "one byte too many");
    let mut magic = good.clone();
    magic[0] = b'O';
    assert!(Checkpoint::from_bytes(&magic).is_err(), "wrong first bytes");
    let mut version = good.clone();
    version[8] = 2;
    assert!(Checkpoint::from_bytes(&version).is_err(), "version 2");
}
