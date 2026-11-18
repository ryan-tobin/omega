//! Measures forces from the GPU (gpu/forces.rs) against the CPU: how far
//! apart the two are, and how long a step takes with each.
//!
//! The sizes and settings are those of the scale program (scale.rs): fcc
//! lattices at Rahman's density with velocities for his temperature, his
//! cutoff and our time step. Each size is first stepped STEPS times on the
//! CPU, with the neighbor list on every thread. At the positions that run
//! ends on, the GPU's forces, energy and virial are compared with the
//! CPU's. Then the same start is stepped STEPS times on the GPU's forces.
//!
//! A pair within rounding of the cutoff can count on the GPU and not on the
//! CPU (gpu/forces.rs), and among a million atoms there are a few. So each
//! atom's force is compared with the CPU's at a cutoff EDGE shorter and at
//! a cutoff EDGE longer, and the nearer of the two counts.
//!
//! A time is the wall-clock time of the steps divided by STEPS, after the
//! system is built and the shader compiled. The GPU's time is also split
//! into its four parts (gpu/forces.rs): sorting the atoms by cell on the
//! CPU, sending them, waiting for the GPU, and reading the answer back.
//! What is left of a step is the velocity Verlet update on the CPU.

use omega::forces::Method;
use omega::gpu::Gpu;
use omega::gpu::forces::GpuForces;
use omega::lattice;
use omega::periodic::PeriodicBox;
use omega::rahman::{CUTOFF, DT, SIDE, TEMPERATURE};
use omega::rng::Rng;
use omega::system::System;
use omega::vec3::Vec3;
use omega::velocities;
use std::thread;
use std::time::Instant;

// Our own choices, the same as in scale.rs.
const SEED: u64 = 1964;
const SKIN: f64 = 0.3;
const STEPS: usize = 100;
const CELLS: [usize; 6] = [6, 10, 16, 25, 40, 63];
const EDGE: f64 = 1e-6;
const ZERO: Vec3 = Vec3::new(0.0, 0.0, 0.0);

/// 4 * cells^3 atoms on an fcc lattice at Rahman's density, with
/// velocities for his temperature, and forces computed by `method`.
fn argon(cells: usize, method: Method) -> System {
    let space = PeriodicBox {
        side: SIDE / 6.0 * cells as f64,
    };
    let positions = lattice::fcc(cells, space.side);
    let mut rng = Rng::new(SEED);
    let velocities = velocities::maxwell_boltzmann(positions.len(), TEMPERATURE, &mut rng);
    System::with_forces(positions, velocities, space, CUTOFF, method)
}

/// Steps `system` STEPS times. Returns the wall-clock time per step, in
/// milliseconds.
fn time(system: &mut System) -> f64 {
    let start = Instant::now();
    for _ in 0..STEPS {
        system.step(DT);
    }
    start.elapsed().as_secs_f64() * 1e3 / STEPS as f64
}

/// The CPU's forces at the positions of `system`, in eps / sigma, with the
/// cutoff `by` sigma longer.
fn reference(system: &System, by: f64, threads: usize) -> Vec<Vec3> {
    let mut forces = vec![ZERO; system.positions.len()];
    let mut method = Method::listed(SKIN, threads);
    method.compute(&system.positions, system.space, CUTOFF + by, &mut forces);
    forces
}

fn main() {
    let Some(gpu) = Gpu::new() else {
        println!("no GPU that wgpu can use on this machine");
        return;
    };
    let threads = thread::available_parallelism().map_or(1, |n| n.get());
    println!("GPU: {}", gpu.name);
    println!("fcc lattice at Rahman's density, box side {SIDE} per 864 atoms");
    println!("cutoff {CUTOFF}, dt {DT}, seed {SEED}, {STEPS} steps");
    println!("CPU: neighbor list, skin {SKIN}, {threads} threads");
    println!();
    println!("GPU against CPU after {STEPS} steps on the CPU, eps / sigma and relative");
    println!("    atoms  largest force  largest miss  miss / largest    energy    virial");
    let mut times = Vec::new();
    for cells in CELLS {
        let mut cpu = argon(cells, Method::listed(SKIN, threads));
        let n = cpu.positions.len();
        let t_cpu = time(&mut cpu);
        let mut on_gpu = GpuForces::new(&gpu);
        let mut forces = vec![ZERO; n];
        let (u, w) = on_gpu.compute(&cpu.positions, cpu.space, CUTOFF, &mut forces);
        let inner = reference(&cpu, -EDGE, threads);
        let outer = reference(&cpu, EDGE, threads);
        let miss = |i: usize| {
            (forces[i] - inner[i])
                .norm()
                .min((forces[i] - outer[i]).norm())
        };
        let miss = (0..n).map(miss).fold(0.0, f64::max);
        let most = cpu.forces.iter().map(|f| f.norm()).fold(0.0, f64::max);
        let du = (u / cpu.potential - 1.0).abs();
        let dw = (w / cpu.virial - 1.0).abs();
        let ratio = miss / most;
        println!("{n:9} {most:14.3} {miss:13.1e} {ratio:15.1e} {du:9.1e} {dw:9.1e}");
        let mut system = argon(cells, Method::Gpu(Box::new(on_gpu)));
        if let Method::Gpu(on_gpu) = &mut system.method {
            on_gpu.seconds = [0.0; 4];
        }
        let t_gpu = time(&mut system);
        if let Method::Gpu(on_gpu) = &system.method {
            let [sort, send, wait, read] = on_gpu.seconds.map(|s| s * 1e3 / STEPS as f64);
            times.push(format!(
                "{n:9} {t_cpu:10.3} {t_gpu:10.3} {sort:8.3} {send:8.3} {wait:8.3} {read:8.3}"
            ));
        }
    }
    println!();
    println!("milliseconds per step");
    println!("    atoms   CPU list GPU forces     sort     send     wait     read");
    for line in times {
        println!("{line}");
    }
}
