//! Brings 864 atoms in Rahman's box to his temperature by velocity
//! rescaling, then lets them move freely and measures them: temperature,
//! potential energy, pressure, g(r), the mean squared displacement and the
//! diffusion coefficient. Writes the measurements as CSV files and the
//! positions as an XYZ file, in the folder `output`.

use omega::csv;
use omega::lattice;
use omega::msd::{self, Msd};
use omega::periodic::PeriodicBox;
use omega::rdf::Rdf;
use omega::rng::Rng;
use omega::system::System;
use omega::thermo;
use omega::units::ARGON;
use omega::velocities;
use omega::xyz;
use std::fs::{self, File};
use std::io::{self, BufWriter};

// Rahman's box, temperature and cutoff in reduced units (PHYSICS_REFERENCE).
// The time step and the seed are ours, as in M2.
const SIDE: f64 = 10.229;
const TEMPERATURE: f64 = 0.7867;
const CUTOFF: f64 = 2.25;
const DT: f64 = 0.005;
const SEED: u64 = 1964;

// Our own choices. 4000 steps (20 tau), rescaled every 10 steps, melt the
// lattice and hold the liquid at the target temperature. Then 10,000 free
// steps (50 tau) are measured every 10 steps, and the positions written
// every 100. g(r) runs out to half the box in 102 bins of 0.05 sigma. D
// is fitted from 5 tau on, after each atom has forgotten its first speed.
const RESCALE_STEPS: usize = 4000;
const RESCALE_EVERY: usize = 10;
const MEASURE_STEPS: usize = 10_000;
const SAMPLE_EVERY: usize = 10;
const XYZ_EVERY: usize = 100;
const RDF_BINS: usize = 102;
const FIT_FROM: f64 = 5.0;

/// The lattice in Rahman's box, melted and brought to his temperature.
fn rahman_liquid() -> System {
    let space = PeriodicBox { side: SIDE };
    let positions = lattice::fcc(6, SIDE);
    let mut rng = Rng::new(SEED);
    let start = velocities::maxwell_boltzmann(positions.len(), TEMPERATURE, &mut rng);
    let mut system = System::new(positions, start, space, CUTOFF);
    for step in 1..=RESCALE_STEPS {
        system.step(DT);
        if step % RESCALE_EVERY == 0 {
            thermo::rescale(&mut system.velocities, TEMPERATURE);
        }
    }
    system
}

fn main() -> io::Result<()> {
    let mut system = rahman_liquid();
    let space = system.space;
    let n = system.positions.len() as f64;
    fs::create_dir_all("output")?;
    let mut thermo_csv = BufWriter::new(File::create("output/thermo.csv")?);
    let mut msd_csv = BufWriter::new(File::create("output/msd.csv")?);
    let mut xyz_file = BufWriter::new(File::create("output/positions.xyz")?);
    // Columns, in reduced units: time (tau), temperature (eps / kB),
    // potential energy per atom (eps), pressure (eps / sigma^3), total
    // energy (eps); MSD in sigma^2; g(r) at r in sigma.
    let names = ["time", "temperature", "potential", "pressure", "energy"];
    csv::write_header(&mut thermo_csv, &names)?;
    csv::write_header(&mut msd_csv, &["time", "msd"])?;

    println!("864 atoms, box side {SIDE}, cutoff {CUTOFF}, dt {DT}, seed {SEED}");
    println!("rescaled to T = {TEMPERATURE} for {RESCALE_STEPS} steps, then free");
    println!("  time       T     U/N        P          E");
    let msd = Msd::new(&system.positions);
    let mut rdf = Rdf::new(space, 0.5 * SIDE, RDF_BINS);
    let (mut fit_t, mut fit_msd) = (Vec::new(), Vec::new());
    let mut sums = [0.0; 3];
    let e0 = system.total_energy();
    let mut worst: f64 = 0.0;
    for step in 0..=MEASURE_STEPS {
        if step > 0 {
            system.step(DT);
        }
        let e = system.total_energy();
        worst = worst.max((e - e0).abs() / e0.abs());
        if step % SAMPLE_EVERY != 0 {
            continue;
        }
        let time = step as f64 * DT;
        let t = thermo::temperature(&system.velocities);
        let u = system.potential / n;
        let p = thermo::pressure(&system.velocities, system.virial, space);
        let m = msd.measure(&system.positions);
        csv::write_row(&mut thermo_csv, &[time, t, u, p, e])?;
        csv::write_row(&mut msd_csv, &[time, m])?;
        rdf.add(&system.positions);
        sums = [sums[0] + t, sums[1] + u, sums[2] + p];
        if time >= FIT_FROM {
            fit_t.push(time);
            fit_msd.push(m);
        }
        if step % XYZ_EVERY == 0 {
            xyz::write_frame(&mut xyz_file, "Ar", &system.positions, space, time, ARGON)?;
        }
        if step % 1000 == 0 {
            println!("{time:6.1} {t:7.4} {u:7.4} {p:8.4} {e:10.4}");
        }
    }
    let mut rdf_csv = BufWriter::new(File::create("output/rdf.csv")?);
    csv::write_header(&mut rdf_csv, &["r", "g"])?;
    for (r, g) in rdf.g() {
        csv::write_row(&mut rdf_csv, &[r, g])?;
    }
    let samples = (MEASURE_STEPS / SAMPLE_EVERY + 1) as f64;
    let [t, u, p] = sums.map(|s| s / samples);
    let kelvin = t * ARGON.temperature_k();
    let mpa = p * ARGON.pressure_pa() / 1e6;
    println!();
    println!("averages over the free run, reduced units and SI");
    println!("temperature            {t:.4} = {kelvin:.2} K");
    println!("potential energy/atom  {u:.4}");
    println!("pressure               {p:.4} = {mpa:.2} MPa");
    println!("largest relative change in total energy  {worst:.3e}");
    report(&rdf, msd::diffusion(&fit_t, &fit_msd));
    Ok(())
}

/// Prints the highest point of g(r) and the diffusion coefficient `d`,
/// in sigma^2 / tau, in reduced units and in SI.
fn report(rdf: &Rdf, d: f64) {
    let mut peak = (0.0, 0.0);
    for (r, g) in rdf.g() {
        if g > peak.1 {
            peak = (r, g);
        }
    }
    let angstrom = ARGON.sigma_m * 1e10;
    let cm = ARGON.sigma_m * 100.0;
    let d_si = d * cm * cm / ARGON.time_s();
    let (r, g) = peak;
    let r_a = r * angstrom;
    println!("highest point of g(r)  {g:.4} at r = {r:.3} = {r_a:.3} angstrom");
    println!("diffusion coefficient  {d:.4} = {d_si:.3e} cm^2/s");
}
