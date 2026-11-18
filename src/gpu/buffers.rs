//! The memory on the GPU that forces.wgsl reads and writes.
//!
//! A buffer is a block of GPU memory. The shader reads three: its
//! parameters, the atoms, and where each cell starts. It writes one, the
//! results. The CPU may not read a buffer that a shader writes, so the
//! results are copied on the GPU into a fifth buffer, `download`, that the
//! CPU maps and reads. The bind group ties the first four to the binding
//! numbers 0 to 3 in the shader.
//!
//! The byte layouts are those of sorted.rs and of the Params of forces.wgsl.

use crate::gpu::sorted::{ATOM_BYTES, RESULT_BYTES};
use wgpu::BufferUsages as Usage;

/// Bytes of the Params of forces.wgsl: three f32 and two u32.
pub const PARAM_BYTES: usize = 20;

/// The buffers for one number of atoms and of cells.
pub struct Buffers {
    pub params: wgpu::Buffer,
    pub atoms: wgpu::Buffer,
    pub starts: wgpu::Buffer,
    pub results: wgpu::Buffer,
    /// The copy of `results` that the CPU reads.
    pub download: wgpu::Buffer,
    pub bind_group: wgpu::BindGroup,
    /// The number of atoms, and of cells along a side, they were made for.
    pub sizes: (usize, usize),
}

impl Buffers {
    /// Buffers for `n` atoms in `m` cells a side, bound to the shader of
    /// `pipeline`.
    pub fn new(
        device: &wgpu::Device,
        pipeline: &wgpu::ComputePipeline,
        n: usize,
        m: usize,
    ) -> Self {
        let make = |size: usize, usage| {
            device.create_buffer(&wgpu::BufferDescriptor {
                label: None,
                size: size as u64,
                usage,
                mapped_at_creation: false,
            })
        };
        let params = make(PARAM_BYTES, Usage::UNIFORM | Usage::COPY_DST);
        let atoms = make(n * ATOM_BYTES, Usage::STORAGE | Usage::COPY_DST);
        let starts = make((m * m * m + 1) * 4, Usage::STORAGE | Usage::COPY_DST);
        let results = make(n * RESULT_BYTES, Usage::STORAGE | Usage::COPY_SRC);
        let download = make(n * RESULT_BYTES, Usage::COPY_DST | Usage::MAP_READ);
        let bound = [&params, &atoms, &starts, &results];
        let entries: Vec<_> = bound
            .iter()
            .enumerate()
            .map(|(binding, buffer)| wgpu::BindGroupEntry {
                binding: binding as u32,
                resource: buffer.as_entire_binding(),
            })
            .collect();
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &pipeline.get_bind_group_layout(0),
            entries: &entries,
        });
        Self {
            params,
            atoms,
            starts,
            results,
            download,
            bind_group,
            sizes: (n, m),
        }
    }
}
