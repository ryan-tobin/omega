//! Lennard-Jones forces computed on the GPU.
//!
//! Each call sorts the atoms by cell on the CPU (sorted.rs), sends them to
//! the GPU, runs forces.wgsl once for every atom, and reads back each
//! atom's force, energy and virial. The time steps stay on the CPU in
//! 64-bit floats (system.rs). Only the forces are computed in 32-bit.
//!
//! The atoms are sorted again at every call, with no skin (neighbors.rs).
//! The forces then depend on the positions alone, not on when the atoms
//! were last sorted, so a run can be retraced and restarted
//! (tests/gpu_runs.rs).
//!
//! Label: approximated. Cut and shifted at r_c as in pairs.rs, and worked
//! out in 32-bit floats: each separation is rounded to about 7 digits, and
//! each atom's sums are rounded as they grow. Measured against all pairs on
//! a software GPU (tests/gpu.rs): each force within 1.6e-6 of the largest
//! force in the system, the energy within 1.5e-7 and the virial within
//! 6.2e-7. A pair within about 1e-7 sigma of r_c can count on the GPU and
//! not on the CPU, or the other way round. That moves two forces by 0.081
//! eps / sigma, the force at r_c, and leaves the energy alone, which is
//! shifted to zero there. The two atoms of a pair get exactly opposite
//! forces (forces.wgsl), but each atom rounds its own sum, so the forces on
//! all atoms no longer add up to exactly zero.

use crate::gpu::Gpu;
use crate::gpu::buffers::{Buffers, PARAM_BYTES};
use crate::gpu::sorted::{RESULT_BYTES, Sorted, TICKS};
use crate::lj;
use crate::periodic::PeriodicBox;
use crate::vec3::Vec3;
use std::time::Instant;

/// Atoms the GPU takes in one group: the workgroup size in forces.wgsl.
const GROUP: usize = 64;

/// Forces from the GPU, for `forces::Method::Gpu`.
pub struct GpuForces {
    device: wgpu::Device,
    queue: wgpu::Queue,
    pipeline: wgpu::ComputePipeline,
    buffers: Option<Buffers>,
    /// Wall-clock seconds spent so far in each part of `compute`: sorting
    /// on the CPU, sending to the GPU, waiting for the GPU, and reading its
    /// answer back into atom order.
    pub seconds: [f64; 4],
}

impl GpuForces {
    /// Compiles forces.wgsl for `gpu`.
    pub fn new(gpu: &Gpu) -> Self {
        let source = wgpu::ShaderSource::Wgsl(include_str!("forces.wgsl").into());
        let module = gpu
            .device
            .create_shader_module(wgpu::ShaderModuleDescriptor {
                label: None,
                source,
            });
        // No layout given: wgpu works it out from the shader's bindings.
        let pipeline = gpu
            .device
            .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: None,
                layout: None,
                module: &module,
                entry_point: Some("forces"),
                compilation_options: Default::default(),
                cache: None,
            });
        Self {
            device: gpu.device.clone(),
            queue: gpu.queue.clone(),
            pipeline,
            buffers: None,
            seconds: [0.0; 4],
        }
    }

    /// Fills `forces` with the total force on each atom, in eps / sigma,
    /// and returns the total potential energy and the virial W, both in
    /// eps, as pairs::compute does. `positions` are in sigma, and `cutoff`
    /// is r_c in sigma, at most half the box side.
    pub fn compute(
        &mut self,
        positions: &[Vec3],
        space: PeriodicBox,
        cutoff: f64,
        forces: &mut [Vec3],
    ) -> (f64, f64) {
        assert_eq!(forces.len(), positions.len());
        assert!(cutoff <= 0.5 * space.side, "cutoff is over half the box");
        let mut times = vec![Instant::now()];
        let sorted = Sorted::new(positions, space, cutoff);
        let (n, m) = (positions.len(), sorted.per_side);
        times.push(Instant::now());
        let buffers = match self.buffers.take() {
            Some(buffers) if buffers.sizes == (n, m) => buffers,
            _ => Buffers::new(&self.device, &self.pipeline, n, m),
        };
        let shift = lj::energy(Vec3::new(cutoff, 0.0, 0.0));
        let mut params = Vec::with_capacity(PARAM_BYTES);
        for x in [space.side / TICKS, cutoff * cutoff, shift] {
            params.extend((x as f32).to_le_bytes());
        }
        for count in [m, n] {
            params.extend((count as u32).to_le_bytes());
        }
        self.queue.write_buffer(&buffers.params, 0, &params);
        self.queue.write_buffer(&buffers.atoms, 0, &sorted.atoms);
        self.queue.write_buffer(&buffers.starts, 0, &sorted.starts);
        let mut encoder = self.device.create_command_encoder(&Default::default());
        {
            let mut pass = encoder.begin_compute_pass(&Default::default());
            pass.set_pipeline(&self.pipeline);
            pass.set_bind_group(0, &buffers.bind_group, &[]);
            pass.dispatch_workgroups(n.div_ceil(GROUP) as u32, 1, 1);
        }
        let size = (n * RESULT_BYTES) as u64;
        encoder.copy_buffer_to_buffer(&buffers.results, 0, &buffers.download, 0, size);
        self.queue.submit([encoder.finish()]);
        times.push(Instant::now());
        // Ask to read the copy, then wait until the GPU has finished.
        let answer = buffers.download.slice(..);
        answer.map_async(wgpu::MapMode::Read, |mapped| {
            mapped.expect("the GPU's answer could not be read")
        });
        let waited = self.device.poll(wgpu::PollType::wait_indefinitely());
        waited.expect("the GPU did not finish");
        times.push(Instant::now());
        let bytes = answer.get_mapped_range().expect("the answer is not ready");
        let sums = sorted.unpack(&bytes, forces);
        drop(bytes);
        buffers.download.unmap();
        self.buffers = Some(buffers);
        times.push(Instant::now());
        for (total, pair) in self.seconds.iter_mut().zip(times.windows(2)) {
            *total += (pair[1] - pair[0]).as_secs_f64();
        }
        sums
    }
}
