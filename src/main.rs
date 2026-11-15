//! Prints the size in SI of each reduced unit for argon.

use omega::units::ARGON;

fn main() {
    println!("Argon reduced units in SI");
    println!("length       {:.4e} m", ARGON.sigma_m);
    println!("energy       {:.4e} J", ARGON.eps_j);
    println!("mass         {:.4e} kg", ARGON.mass_kg);
    println!("time         {:.4e} s", ARGON.time_s());
    println!("temperature  {:.4e} K", ARGON.temperature_k());
    println!("velocity     {:.4e} m/s", ARGON.velocity_m_per_s());
    println!("pressure     {:.4e} Pa", ARGON.pressure_pa());
}
