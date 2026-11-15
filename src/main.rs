//! Prints the size in SI of each reduced unit for argon, the Lennard-Jones
//! energy and force between two atoms, and the energy minimum in SI.

use omega::lj;
use omega::units::ARGON;
use omega::vec3::Vec3;

fn main() {
    println!("Argon reduced units in SI");
    println!("length       {:.4e} m", ARGON.sigma_m);
    println!("energy       {:.4e} J", ARGON.eps_j);
    println!("mass         {:.4e} kg", ARGON.mass_kg);
    println!("time         {:.4e} s", ARGON.time_s());
    println!("temperature  {:.4e} K", ARGON.temperature_k());
    println!("velocity     {:.4e} m/s", ARGON.velocity_m_per_s());
    println!("pressure     {:.4e} Pa", ARGON.pressure_pa());

    // Two atoms on the x axis, a distance r apart, in reduced units: r in
    // sigma, U in eps, F in eps / sigma. A positive force pushes the atoms
    // apart. A negative force pulls them together.
    println!();
    println!("Two atoms, Lennard-Jones, reduced units");
    println!("     r         U         F");
    for r in [0.95, 1.0, 1.1, 1.2, 1.5, 2.0, 2.5] {
        let r_ij = Vec3::new(r, 0.0, 0.0);
        let u = lj::energy(r_ij);
        let f = lj::force(r_ij).x;
        println!("{r:6.2} {u:9.4} {f:9.4}");
    }

    // The bottom of the energy well at r = 2^(1/6), converted to SI for
    // argon by multiplying by the length unit and the energy unit.
    let r_min = 2.0_f64.powf(1.0 / 6.0);
    let u_min = lj::energy(Vec3::new(r_min, 0.0, 0.0));
    println!();
    println!("Lennard-Jones energy minimum, in SI for argon");
    println!("separation   {:.4e} m", r_min * ARGON.sigma_m);
    println!("energy       {:.4e} J", u_min * ARGON.eps_j);
}
