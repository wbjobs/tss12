use crate::types::EigenResult;
use wgpu::util::DeviceExt;
use bytemuck::{Pod, Zeroable};
use std::time::Instant;

#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
struct GPUParams {
    n: u32,
    k: u32,
    iter: u32,
    padding: u32,
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

        let pipeline = self
            .device
            .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some("Eigen Pipeline"),
                layout: None,
                module: &shader,
                entry_point,
                compilation_options: Default::default(),
                cache: None,
            });
        pipeline
    }

    pub fn solve_eigen(
        &self,
        matrix: &Vec<Vec<f64>>,
        num_eigen: Option<usize>,
    ) -> EigenResult {
        let n = matrix.len();
        let k = num_eigen.unwrap_or(n.min(100)).min(n);

        if n <= 200 {
            return crate::eigen_cpu::solve_eigen_jacobi(matrix, 500, 1e-10);
        }

        pollster::block_on(self.solve_eigen_gpu(matrix, n, k))
            .unwrap_or_else(|_| crate::eigen_cpu::solve_eigen_jacobi(matrix, 500, 1e-10))
    }

    async fn solve_eigen_gpu(
        &self,
        matrix: &Vec<Vec<f64>>,
        n: usize,
        k: usize,
    ) -> Result<EigenResult, String> {
        println!("[GPU] Running Lanczos eigensolver for {}x{} matrix, requesting {} eigenpairs", n, n, k);
        let start = Instant::now();

        let mut a_flat: Vec<f32> = matrix.iter().flatten().map(|&x| x as f32).collect();
        self.lanczos_normalize(&mut a_flat, n)?;

        let (eigenvalues, eigenvectors) = self.lanczos_iteration(&a_flat, n, k).await?;

        println!("[GPU] Eigensolver completed in {:.2?}", start.elapsed());

        let mut sorted_pairs: Vec<(f64, Vec<f64>)> = eigenvalues
            .into_iter()
            .zip(eigenvectors.into_iter())
            .collect();
        sorted_pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        Ok(EigenResult {
            eigenvalues: sorted_pairs.iter().map(|x| x.0).collect(),
            eigenvectors: sorted_pairs.into_iter().map(|x| x.1).collect(),
        })
    }

    fn lanczos_normalize(&self, a: &mut Vec<f32>, n: usize) -> Result<(), String> {
        let mut max_val = 0.0f32;
        for &val in a.iter() {
            max_val = max_val.max(val.abs());
        }
        if max_val > 0.0 {
            for val in a.iter_mut() {
                *val /= max_val;
            }
        }
        let trace: f32 = (0..n).map(|i| a[i * n + i]).sum();
        let avg = trace / n as f32;
        for i in 0..n {
            a[i * n + i] -= avg;
        }
        Ok(())
    }

    async fn lanczos_iteration(
        &self,
        a_flat: &Vec<f32>,
        n: usize,
        k: usize,
    ) -> Result<(Vec<f64>, Vec<Vec<f64>>), String> {
        let m = (k * 3).min(n).max(k + 10);
        println!("[GPU] Lanczos: using {} iterations for {} eigenvalues", m, k);

        let mut v = vec![vec![0.0f32; n]; m + 1];
        let mut alpha = vec![0.0f32; m];
        let mut beta = vec![0.0f32; m];

        for j in 0..n {
            v[0][j] = (1.0 / (n as f32)).sqrt();
        }

        for j in 1..=m {
            let w = self.gpu_mat_vec(a_flat, &v[j - 1], n).await?;

            let mut w_vec = w;
            if j > 1 {
                for i in 0..n {
                    w_vec[i] -= beta[j - 2] * v[j - 2][i];
                }
            }

            let mut aj = 0.0f32;
            for i in 0..n {
                aj += v[j - 1][i] * w_vec[i];
            }
            alpha[j - 1] = aj;

            for i in 0..n {
                w_vec[i] -= aj * v[j - 1][i];
            }

            if j < m {
                let bj2: f32 = w_vec.iter().map(|x| x * x).sum();
                let bj = bj2.sqrt();
                beta[j - 1] = bj;
                if bj < 1e-12 {
                    return self.reorthogonalize_and_extract(&v, &alpha, &beta, j, k, n);
                }
                let inv_bj = 1.0 / bj;
                for i in 0..n {
                    v[j][i] = w_vec[i] * inv_bj;
                }
            }
        }

        self.reorthogonalize_and_extract(&v, &alpha, &beta, m, k, n)
    }

    async fn gpu_mat_vec(
        &self,
        a: &Vec<f32>,
        x: &Vec<f32>,
        n: usize,
    ) -> Result<Vec<f32>, String> {
        if n < 128 {
            let mut y = vec![0.0f32; n];
            for i in 0..n {
                let mut sum = 0.0f32;
                for j in 0..n {
                    sum += a[i * n + j] * x[j];
                }
                y[i] = sum;
            }
            return Ok(y);
        }

        let shader_src = include_str!("shaders/matvec.wgsl");
        let pipeline = self.create_compute_pipeline(shader_src, "main");

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

        let y_buf = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Vector Y"),
            size: (n * std::mem::size_of::<f32>()) as u64,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });

        let read_buf = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Read Buffer"),
            size: (n * std::mem::size_of::<f32>()) as u64,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });

        let params = GPUParams {
            n: n as u32,
            k: 0,
            iter: 0,
            padding: 0,
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
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: params_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: a_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: x_buf.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: y_buf.as_entire_binding(),
                },
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
            cpass.dispatch_workgroups(workgroups, 1, 1);
        }

        encoder.copy_buffer_to_buffer(&y_buf, 0, &read_buf, 0, (n * std::mem::size_of::<f32>()) as u64);
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

    fn reorthogonalize_and_extract(
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
            if i + 1 < m_actual {
                t[i][i + 1] = beta[i] as f64;
                t[i + 1][i] = beta[i] as f64;
            }
        }

        let small_result = crate::eigen_cpu::solve_eigen_jacobi(&t, 1000, 1e-12);

        let num_extract = k.min(small_result.eigenvalues.len());
        let mut ritz_vectors = vec![vec![0.0f64; n]; num_extract];

        for eig_idx in 0..num_extract {
            for atom_idx in 0..n {
                let mut sum = 0.0f64;
                for lanczos_idx in 0..m_actual {
                    if lanczos_idx < v.len() && atom_idx < v[lanczos_idx].len() {
                        sum += v[lanczos_idx][atom_idx] as f64 * small_result.eigenvectors[eig_idx][lanczos_idx];
                    }
                }
                ritz_vectors[eig_idx][atom_idx] = sum;
            }

            let norm: f64 = ritz_vectors[eig_idx].iter().map(|x| x * x).sum::<f64>().sqrt();
            if norm > 1e-15 {
                let inv_norm = 1.0 / norm;
                for val in ritz_vectors[eig_idx].iter_mut() {
                    *val *= inv_norm;
                }
            }
        }

        let mut eigenvalues: Vec<f64> = small_result.eigenvalues
            .iter()
            .take(num_extract)
            .copied()
            .collect();

        while eigenvalues.len() < k {
            eigenvalues.push(0.0);
        }
        while ritz_vectors.len() < k {
            ritz_vectors.push(vec![0.0f64; n]);
        }

        Ok((eigenvalues, ritz_vectors))
    }
}
