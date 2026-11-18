// Lennard-Jones forces on the GPU, in reduced units (sigma, eps), written in
// WGSL. The GPU runs forces() once for each atom, many atoms at a time.
//
// The atoms arrive sorted by cell (sorted.rs): the box is cut into m x m x m
// cells, each at least r_c wide, so every atom within r_c of an atom sits in
// the atom's own cell or in one of the 26 around it. Cells across a face of
// the box are reached by counting cells modulo m.
//
// Each position is three whole numbers of ticks, a tick being 2^-32 of the
// box side. Whole numbers wrap around at 2^32 as the box does at its side,
// so the difference of two positions, read as signed, is the separation to
// the nearest copy. It is exact, and swapping the two atoms gives exactly
// the opposite, so the two atoms of a pair get exactly opposite forces, on
// any GPU that rounds a number and its negative alike. Nothing after that
// difference is exact: a GPU has no 64-bit floats here, so the rest is in
// 32-bit floats, about 7 digits.
//
// Each pair is worked out twice, once for each of its atoms, as on the CPU
// with a neighbor list (forces.rs).

struct Params {
    tick: f32,      // one tick, in sigma: the box side / 2^32
    cutoff_sq: f32, // r_c squared, in sigma^2
    shift: f32,     // U(r_c), in eps: the energy is cut and shifted (pairs.rs)
    per_side: u32,  // cells along each side of the box, m
    atoms: u32,     // number of atoms
}

struct Atom {
    r: vec3<u32>, // position, in ticks
    cell: u32,    // its cell, (x * m + y) * m + z
}

@group(0) @binding(0) var<uniform> params: Params;
@group(0) @binding(1) var<storage, read> atoms: array<Atom>;
// The atoms of cell c are atoms[starts[c]] up to, not including,
// atoms[starts[c + 1]].
@group(0) @binding(2) var<storage, read> starts: array<u32>;
// Five numbers for each atom: the force on it (x, y, z, in eps / sigma), and
// its sums of pair energies (eps) and of r . F (eps) over its neighbors.
@group(0) @binding(3) var<storage, read_write> results: array<f32>;

@compute @workgroup_size(64)
fn forces(@builtin(global_invocation_id) id: vec3<u32>) {
    let i = id.x;
    // Atoms are handed out in groups of 64, so the last group can run past
    // the end.
    if (i >= params.atoms) {
        return;
    }
    let m = params.per_side;
    let me = atoms[i];
    let c = vec3<u32>(me.cell / (m * m), me.cell / m % m, me.cell % m);
    // With 2 cells a side, the cell one step down is also the cell one step
    // up. It is met once: 2 cells along each axis, not 3.
    let reach = min(m, 3u);
    var f = vec3<f32>(0.0, 0.0, 0.0);
    var u = 0.0;
    var w = 0.0;
    for (var dx = 0u; dx < reach; dx++) {
        for (var dy = 0u; dy < reach; dy++) {
            for (var dz = 0u; dz < reach; dz++) {
                // The neighbor cell is c - 1, c or c + 1 along each axis,
                // modulo m, kept unsigned.
                let n = (c + m + vec3<u32>(dx, dy, dz) - 1u) % m;
                let cell = (n.x * m + n.y) * m + n.z;
                for (var j = starts[cell]; j < starts[cell + 1u]; j++) {
                    // r_i - r_j to the nearest copy: in ticks, then sigma.
                    let ticks = vec3<i32>(me.r - atoms[j].r);
                    let r = vec3<f32>(ticks) * params.tick;
                    let r2 = dot(r, r);
                    if (r2 < params.cutoff_sq && j != i) {
                        // lj.rs: F = 24 (2 / r^12 - 1 / r^6) r / r^2 and
                        // U = 4 (1 / r^12 - 1 / r^6).
                        let inv2 = 1.0 / r2;
                        let inv6 = inv2 * inv2 * inv2;
                        let pull = 24.0 * (2.0 * inv6 * inv6 - inv6) * inv2;
                        f += r * pull;
                        u += 4.0 * (inv6 * inv6 - inv6) - params.shift;
                        w += r2 * pull;
                    }
                }
            }
        }
    }
    let at = 5u * i;
    results[at] = f.x;
    results[at + 1u] = f.y;
    results[at + 2u] = f.z;
    results[at + 3u] = u;
    results[at + 4u] = w;
}