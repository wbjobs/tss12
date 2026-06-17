struct Params {
    n: u32,
    k: u32,
    iter: u32,
    padding: u32,
}

@group(0) @binding(0)
var<uniform> params: Params;

@group(0) @binding(1)
var<storage, read> matrix_a: array<f32>;

@group(0) @binding(2)
var<storage, read> vector_x: array<f32>;

@group(0) @binding(3)
var<storage, read_write> vector_y: array<f32>;

@compute @workgroup_size(256, 1, 1)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i: u32 = gid.x;
    let n: u32 = params.n;

    if (i >= n) {
        return;
    }

    var sum: f32 = 0.0;
    let row_base: u32 = i * n;

    var j: u32 = 0u;
    loop {
        if (j >= n) {
            break;
        }
        let remaining: u32 = n - j;
        let block: u32 = select(32u, remaining, remaining >= 32u);

        var block_sum: f32 = 0.0;
        var jj: u32 = 0u;
        loop {
            if (jj >= block) {
                break;
            }
            let idx: u32 = j + jj;
            block_sum = block_sum + matrix_a[row_base + idx] * vector_x[idx];
            jj = jj + 1u;
        }
        sum = sum + block_sum;

        j = j + block;
    }

    vector_y[i] = sum;
}
