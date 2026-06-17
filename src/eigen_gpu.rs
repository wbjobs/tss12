use crate::types::{EigenResult, Algorithm};
use crate::eigen_cpu::{symmetrize, sanitize_matrix, regularize_positive_definite, solve_generalized_eigen, solve_eigen_jacobi};
use wgpu::util::DeviceExt;
use bytemuck::{Pod, Zeroable};
use std::time::Instant;

#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
struct GPUParams {
    n: u32,
    matrix_len: u32,
    x_len: u32,
    y_len: u32,
}

pub struct GPUEigenSolver {
    device: wgpu::Device,
    queue: wgpu::Queue,
}

impl GPUEigenSolver {
    pub async fn new() -> Result<Self, String> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: wgpu::Backends::all(),
            dx12_shader_compiler: Default::default(),
            gles_minor_version: Default::default(),
        });

        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                compatible_surface: None,
                force_fallback_adapter: false,
            })
            .await
            .map_err(|e| format!("No GPU adapter found: {}", e))?;

        let (device, queue) = adapter
            .request_device(
                &wgpu::DeviceDescriptor {
                    label: Some("Huckel GPU Device"),
                    features: wgpu::Features::empty(),
                    limits: wgpu::Limits::default(),
                },
                None,
            )
            .await
            .map_err(|e| format!("Failed to create GPU device: {}", e))?;

        Ok(Self { device, queue })
    }

    fn create_compute_pipeline(&self, shader_src: &str, entry_point: &str) -> wgpu::ComputePipeline {
        let shader = self
            .device
            .create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("Eigen Shader"),
                source: wgpu::ShaderSource::Wgsl(shader_src.into()),
            });

        self.device
            .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some("Eigen Pipeline"),
                layout: None,
                module: &shader,
                entry_point,
                compilation_options: Default::default(),
                cache: None,
            })
    }

    pub fn solve_eigen(
        &self,
        matrix: &Vec<Vec<f64>>,
        overlap: Option<&Vec<Vec<f64>>>,
        num_eigen: Option<usize>,
    ) -> EigenResult {
        let n = matrix.len();

        if let Some(s) = overlap {
            if n <= 300 {
                return solve_generalized_eigen(matrix, s, 2000, 1e-10);
            }
            println!("[GPU] Generalized eigenvalue problem (>300) falling back to CPU for stability");
            return solve_generalized_eigen(matrix, s, 2000, 1e-10);
        }

        let k = num_eigen.unwrap_or(n.min(120)).min(n);
        if n <= 250 {
            return solve_eigen_jacobi(matrix, 1000, 1e-10);
        }

        pollster::block_on(self.solve_eigen_gpu(matrix, n, k))
            .unwrap_or_else(|e| {
                eprintln!("[GPU] Failed ({}), falling back to CPU", e);
                solve_eigen_jacobi(matrix, 1000, 1e-10)
            })
    }

    async fn solve_eigen_gpu(
        &self,
        matrix: &Vec<Vec<f64>>,
        n: usize,
        k: usize,
    ) -> Result<EigenResult, String> {
        println!("[GPU] Lanczos eigensolver: {}x{} matrix, requesting {} eigenpairs", n, n, k);
        let start = Instant::now();

        let mut m_clone = matrix.clone();
        sanitize_matrix(&mut m_clone);
        symmetrize(&mut m_clone);

        let mut a_flat: Vec<f32> = Vec::with_capacity(n * n);
        let mut max_val = 0.0f32;
        for row in &m_clone {
            for &val in row {
                let v = val as f32;
                if v.abs() > max_val { max_val = v.abs(); }
                a_flat.push(v);
            }
        }
        if max_val < 1e-20 { max_val = 1.0; }
        let inv_max = 1.0 / max_val;
        for v in &mut a_flat { *v *= inv_max; }
        let trace: f32 = (0..n).map(|i| a_flat[i * n + i]).sum();
        let avg = trace / n as f32;
        for i in 0..n { a_flat[i * n + i] -= avg; }

        let actual_matrix_len = a_flat.len() as u32;
        let (eigenvalues, eigenvectors) = self.lanczos_iteration(&a_flat, n, k, actual_matrix_len).await?;

        let mut ev_scaled: Vec<f64> = eigenvalues.iter().map(|&e| e as f64 * max_val as f64 + avg as f64).collect();

        println!("[GPU] Eigensolver done in {:.2?}", start.elapsed());

        let mut pairs: Vec<(f64, Vec<f64>)> = ev_scaled
            .drain(..)
            .zip(eigenvectors.into_iter())
            .collect();
        pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        Ok(EigenResult {
            eigenvalues: pairs.iter().map(|x| x.0).collect(),
            eigenvectors: pairs.into_iter().map(|x| x.1).collect(),
            algorithm: Algorithm::SimpleHuckel,
        })
    }

    async fn lanczos_iteration(
        &self,
        a_flat: &Vec<f32>,
        n: usize,
        k: usize,
        actual_matrix_len: u32,
    ) -> Result<(Vec<f64>, Vec<Vec<f64>>), String> {
        let m = (k * 3).min(n).max(k + 20);
        println!("[GPU] Lanczos: {} iterations for {} eigenvalues", m, k);

        let mut v: Vec<Vec<f32>> = vec![vec![0.0f32; n]; m + 1];
        let mut alpha = vec![0.0f32; m];
        let mut beta = vec![0.0f32; m];

        for j in 0..n {
            v[0][j] = (1.0 / (n as f32)).sqrt();
        }

        for j in 1..=m {
            let w = self.gpu_mat_vec_safe(a_flat, &v[j - 1], n, actual_matrix_len).await
                .unwrap_or_else(|_| self.cpu_mat_vec(a_flat, &v[j - 1], n));

            let mut w_vec = w;
            if w_vec.len() != n { w_vec.resize(n, 0.0); }
            if j > 1 {
                for i in 0..n {
                    let b = if j - 2 < beta.len() { beta[j - 2] } else { 0.0 };
                    w_vec[i] -= b * v[j - 2][i];
                }
            }

            let mut aj = 0.0f32;
            for i in 0..n {
                aj += v[j - 1][i] * w_vec[i];
            }
            if j - 1 < alpha.len() { alpha[j - 1] = aj; }
            for i in 0..n { w_vec[i] -= aj * v[j - 1][i]; }

            if j < m {
                let bj2: f32 = w_vec.iter().map(|x| x * x).sum();
                let bj = bj2.sqrt();
                beta[j - 1] = bj;
                if bj < 1e-12 {
                    return self.extract_ritz(&v, &alpha, &beta, j.min(m), k, n);
                }
                let inv = 1.0 / bj;
                for i in 0..n { v[j][i] = w_vec[i] * inv; }
            }
        }
        self.extract_ritz(&v, &alpha, &beta, m, k, n)
    }

    fn cpu_mat_vec(&self, a: &Vec<f32>, x: &Vec<f32>, n: usize) -> Vec<f32> {
        let mut y = vec![0.0f32; n];
        for i in 0..n {
            let mut sum = 0.0f32;
            let rb = i * n;
            for j in 0..n {
                if rb + j < a.len() && j < x.len() {
                    sum += a[rb + j] * x[j];
                }
            }
            y[i] = if sum.is_finite() { sum } else { 0.0 };
        }
        y
    }

    async fn gpu_mat_vec_safe(
        &self,
        a: &Vec<f32>,
        x: &Vec<f32>,
        n: usize,
        actual_matrix_len: u32,
    ) -> Result<Vec<f32>, String> {
        if n < 256 || actual_matrix_len != (n * n) as u32 {
            return Ok(self.cpu_mat_vec(a, x, n));
        }
        self.gpu_mat_vec(a, x, n, actual_matrix_len).await
            .map_err(|e| format!("GPU matvec: {}", e))
    }

    async fn gpu_mat_vec(
        &self,
        a: &Vec<f32>,
        x: &Vec<f32>,
        n: usize,
        actual_matrix_len: u32,
    ) -> Result<Vec<f32>, String> {
        let shader_src = include_str!("shaders/matvec.wgsl");
        let pipeline = self.create_compute_pipeline(shader_src, "main");

        let a_len = a.len() as u32;
        let x_len = x.len() as u32;
        let y_len = n as u32;

        let a_buf = self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Matrix A"),
            contents: bytemuck::cast_slice(a),
            usage: wgpu::BufferUsages::STORAGE,
        });

        let x_buf = self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Vector X"),
            contents: bytemuck::cast_slice(x),
            usage: wgpu::BufferUsages::STORAGE,
        });

        let y_bytes = (n * std::mem::size_of::<f32>()) as u64;
        let y_buf = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Vector Y"),
            size: y_bytes,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });

        let read_buf = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Read Buffer"),
            size: y_bytes,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });

        let params = GPUParams {
            n: n as u32,
            matrix_len: a_len.min(actual_matrix_len),
            x_len,
            y_len,
        };
        let params_buf = self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("Params"),
            contents: bytemuck::bytes_of(&params),
            usage: wgpu::BufferUsages::UNIFORM,
        });

        let bind_group = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("MatVec Bind Group"),
            layout: &pipeline.get_bind_group_layout(0),
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: params_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: a_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: x_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: y_buf.as_entire_binding() },
            ],
        });

        let mut encoder = self.device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("MatVec Encoder"),
        });
        {
            let mut cpass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("MatVec Pass"),
                timestamp_writes: None,
            });
            cpass.set_pipeline(&pipeline);
            cpass.set_bind_group(0, &bind_group, &[]);
            let workgroups = ((n as u32) + 255) / 256;
            cpass.dispatch_workgroups(workgroups.max(1), 1, 1);
        }
        encoder.copy_buffer_to_buffer(&y_buf, 0, &read_buf, 0, y_bytes);
        self.queue.submit(Some(encoder.finish()));

        let slice = read_buf.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |v| tx.send(v).unwrap());
        self.device.poll(wgpu::Maintain::wait()).panic_on_timeout();
        rx.recv().unwrap().map_err(|e| format!("Map error: {}", e))?;

        let data = slice.get_mapped_range();
        let result: Vec<f32> = bytemuck::cast_slice(&data).to_vec();
        drop(data);
        read_buf.unmap();

        Ok(result)
    }

    fn extract_ritz(
        &self,
        v: &Vec<Vec<f32>>,
        alpha: &Vec<f32>,
        beta: &Vec<f32>,
        m: usize,
        k: usize,
        n: usize,
    ) -> Result<(Vec<f64>, Vec<Vec<f64>>), String> {
        let m_actual = m.min(alpha.len()).max(1);
        let mut t = vec![vec![0.0f64; m_actual]; m_actual];
        for i in 0..m_actual {
            t[i][i] = alpha[i] as f64;
            if i + 1 < m_actual && i + 1 < beta.len() {
                t[i][i + 1] = beta[i] as f64;
                t[i + 1][i] = beta[i] as f64;
            }
        }

        let small = solve_eigen_jacobi(&t, 1500, 1e-12);

        let num_extract = k.min(small.eigenvalues.len()).max(1);
        let mut ritz = vec![vec![0.0f64; n]; num_extract];
        for eig in 0..num_extract {
            for atom_idx in 0..n {
                let mut sum = 0.0f64;
                for li in 0..m_actual {
                    if li < v.len() && atom_idx < v[li].len() && eig < small.eigenvectors.len() && li < small.eigenvectors[eig].len() {
                        sum += v[li][atom_idx] as f64 * small.eigenvectors[eig][li];
                    }
                }
                ritz[eig][atom_idx] = sum;
            }
            let norm_sq: f64 = ritz[eig].iter().map(|x| x * x).sum();
            if norm_sq > 1e-20 {
                let inv = 1.0 / norm_sq.sqrt();
                for val in ritz[eig].iter_mut() { *val *= inv; }
            }
        }

        let mut ev: Vec<f64> = small.eigenvalues.iter().take(num_extract).copied().collect();
        while ev.len() < k { ev.push(0.0); }
        while ritz.len() < k { ritz.push(vec![0.0f64; n]); }

        Ok((ev, ritz))
    }

    pub fn batch_solve_eigen(
        &self,
        matrices: &Vec<Vec<Vec<f64>>>,
        overlaps: &Vec<Option<Vec<Vec<f64>>>>,
        use_gpu: bool,
    ) -> Vec<EigenResult> {
        let n_items = matrices.len();
        println!("[Batch] Computing {} eigenvalue problems", n_items);

        let mut results = Vec::with_capacity(n_items);
        for (idx, matrix) in matrices.iter().enumerate() {
            let n = matrix.len();
            let result = if let Some(ref s) = overlaps.get(idx).and_then(|o| o.as_ref()) {
                solve_generalized_eigen(matrix, s, 2000, 1e-10)
            } else if use_gpu && n > 250 {
                match pollster::block_on(self.solve_eigen_gpu(matrix, n, n.min(80))) {
                    Ok(r) => r,
                    Err(_) => solve_eigen_jacobi(matrix, 1000, 1e-10),
                }
            } else {
                solve_eigen_jacobi(matrix, 1000, 1e-10)
            };
            results.push(result);
            if (idx + 1) % 5 == 0 || idx + 1 == n_items {
                println!("  [Batch] {}/{} steps completed", idx + 1, n_items);
            }
        }
        results
    }
}
