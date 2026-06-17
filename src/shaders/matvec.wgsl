struct Params {
    n: u32,
    matrix_len: u32,
    x_len: u32,
    y_len: u32,
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
    let mlen: u32 = params.matrix_len;
    let xlen: u32 = params.x_len;
    let ylen: u32 = params.y_len;

    if (i >= n || i >= ylen) {
        return;
    }

    var sum: f32 = 0.0;
    var row_base: u32 = i * n;

    if (row_base >= mlen && n > 0u) {
        vector_y[i] = 0.0;
        return;
    }

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
            let mat_idx: u32 = row_base + idx;

            if (mat_idx < mlen && idx < xlen) {
                let a: f32 = matrix_a[mat_idx];
                let x: f32 = vector_x[idx];
                if (isFinite(a) && isFinite(x)) {
                    block_sum = block_sum + a * x;
                }
            }

            jj = jj + 1u;
        }
        if (isFinite(block_sum)) {
            sum = sum + block_sum;
        }

        j = j + block;
    }

    if (!isFinite(sum)) {
        sum = 0.0;
    }

    if (i < ylen) {
        vector_y[i] = sum;
    }
}
