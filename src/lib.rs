pub mod types;
pub mod huckel;
pub mod eigen_cpu;
pub mod eigen_gpu;
pub mod svg_gen;
pub mod output;
pub mod html_template;

use pyo3::prelude::*;
use pyo3::exceptions::PyRuntimeError;
use std::fs;
use std::path::Path;
use std::time::Instant;

use types::{Molecule, CalculationOutput, EigenResult};
use huckel::{build_hamiltonian, build_overlap_matrix, compute_electron_density};
use eigen_cpu::solve_eigen_jacobi;
use eigen_gpu::GPUEigenSolver;
use svg_gen::generate_energy_level_svg;
use output::{write_density_binary, density_to_base64, molecule_to_json};
use html_template::get_html_template;

#[pyclass]
pub struct HuckelEngine {
    use_gpu: bool,
    gpu_solver: Option<GPUEigenSolver>,
}

#[pymethods]
impl HuckelEngine {
    #[new]
    #[pyo3(signature = (use_gpu=true))]
    fn new(use_gpu: bool) -> Self {
        Self {
            use_gpu,
            gpu_solver: None,
        }
    }

    fn initialize_gpu<'py>(&mut self, py: Python<'py>) -> PyResult<()> {
        if !self.use_gpu {
            return Ok(());
        }
        if self.gpu_solver.is_some() {
            return Ok(());
        }
        println!("[Engine] Initializing GPU accelerator (wgpu)...");
        let solver = py.allow_threads(|| {
            pollster::block_on(GPUEigenSolver::new())
        }).map_err(|e| PyRuntimeError::new_err(format!("GPU init failed: {}. Will use CPU fallback.", e)))?;
        println!("[Engine] GPU initialized successfully");
        self.gpu_solver = Some(solver);
        Ok(())
    }

    #[pyo3(signature = (mol_json, output_dir=None, grid_resolution=50, grid_padding=3.0))]
    fn calculate<'py>(
        &mut self,
        py: Python<'py>,
        mol_json: &str,
        output_dir: Option<&str>,
        grid_resolution: usize,
        grid_padding: f64,
    ) -> PyResult<PyObject> {
        println!("╔══════════════════════════════════════════════╗");
        println!("║    Hückel Molecular Orbital Calculation      ║");
        println!("╚══════════════════════════════════════════════╝");

        let mol: Molecule = serde_json::from_str(mol_json)
            .map_err(|e| PyRuntimeError::new_err(format!("Failed to parse molecule JSON: {}", e)))?;

        let n = mol.num_atoms();
        let num_electrons = mol.num_pi_electrons();
        println!("  Atoms:        {}", n);
        println!("  π-electrons:  {}", num_electrons);
        println!("  Charge:       {:?}", mol.charge.unwrap_or(0));

        let _ = self.initialize_gpu(py);

        println!("\n[Step 1] Building Hamiltonian matrix...");
        let t0 = Instant::now();
        let h = build_hamiltonian(&mol);
        let s = build_overlap_matrix(&mol);
        println!("  Matrix size: {}x{} ({:.2} MB)",
                 n, n, (n * n * 8) as f64 / (1024.0 * 1024.0));
        println!("  Completed in {:.2?}", t0.elapsed());

        println!("\n[Step 2] Solving eigenvalue problem Hc = εSc...");
        let t1 = Instant::now();

        let num_eigen = (num_electrons / 2 + 10).max(n.min(50));

        let result: EigenResult = if let Some(ref gpu) = self.gpu_solver {
            if n > 200 {
                println!("  Using GPU (Lanczos method) for {} eigenpairs", num_eigen);
                py.allow_threads(|| gpu.solve_eigen(&h, Some(num_eigen)))
            } else {
                println!("  Using CPU (Jacobi method) for small matrix");
                py.allow_threads(|| solve_eigen_jacobi(&h, 500, 1e-10))
            }
        } else {
            println!("  Using CPU (Jacobi method)");
            py.allow_threads(|| solve_eigen_jacobi(&h, 500, 1e-10))
        };

        println!("  Found {} eigenvalues", result.eigenvalues.len());
        println!("  Completed in {:.2?}", t1.elapsed());

        let homo_idx = if num_electrons > 0 { (num_electrons - 1) / 2 } else { 0 };
        let lumo_idx = homo_idx + 1;
        let gap = if lumo_idx < result.eigenvalues.len() && homo_idx < result.eigenvalues.len() {
            result.eigenvalues[lumo_idx] - result.eigenvalues[homo_idx]
        } else { 0.0 };

        println!("\n  Key Results:");
        println!("    HOMO (ε_{}) = {:.6} eV", homo_idx + 1, result.eigenvalues[homo_idx.min(result.eigenvalues.len() - 1)]);
        if lumo_idx < result.eigenvalues.len() {
            println!("    LUMO (ε_{}) = {:.6} eV", lumo_idx + 1, result.eigenvalues[lumo_idx]);
        }
        println!("    HOMO-LUMO Gap = {:.6} eV", gap);

        println!("\n[Step 3] Computing electron density grid...");
        let t2 = Instant::now();

        let (bb_min, bb_max) = mol.bounding_box();
        let mut grid_origin = [
            bb_min[0] - grid_padding,
            bb_min[1] - grid_padding,
            bb_min[2] - grid_padding,
        ];
        let grid_size = [
            (bb_max[0] - bb_min[0]) + 2.0 * grid_padding,
            (bb_max[1] - bb_min[1]) + 2.0 * grid_padding,
            (bb_max[2] - bb_min[2]) + 2.0 * grid_padding,
        ];
        let max_size = grid_size[0].max(grid_size[1]).max(grid_size[2]);
        let grid_spacing = max_size / grid_resolution as f64;
        let grid_dims = [
            (grid_size[0] / grid_spacing).ceil() as usize + 1,
            (grid_size[1] / grid_spacing).ceil() as usize + 1,
            (grid_size[2] / grid_spacing).ceil() as usize + 1,
        ];
        let total_voxels = grid_dims[0] * grid_dims[1] * grid_dims[2];

        println!("  Grid dimensions: {} x {} x {} = {} voxels",
                 grid_dims[0], grid_dims[1], grid_dims[2], total_voxels);
        println!("  Grid spacing:    {:.4} Å", grid_spacing);
        println!("  Memory estimate: {:.2} MB", (total_voxels * 4) as f64 / (1024.0 * 1024.0));

        let density = py.allow_threads(|| compute_electron_density(
            &mol,
            &result.eigenvectors,
            num_electrons,
            grid_origin,
            grid_dims,
            grid_spacing,
        ));
        println!("  Completed in {:.2?}", t2.elapsed());

        let calc_output = CalculationOutput {
            eigenvalues: result.eigenvalues.clone(),
            homo_lumo_gap: gap,
            electron_density_grid: density,
            grid_dims,
            grid_origin,
            grid_spacing,
            num_electrons,
        };

        if let Some(dir) = output_dir {
            println!("\n[Step 4] Generating output files...");
            let t3 = Instant::now();
            let path = Path::new(dir);
            fs::create_dir_all(path).ok();

            let svg_path = path.join("energy_levels.svg");
            let svg = generate_energy_level_svg(&result.eigenvalues, num_electrons, 900, 700);
            fs::write(&svg_path, &svg).ok();
            println!("  ✓ SVG energy level diagram: {}", svg_path.display());

            let bin_path = path.join("electron_density.bin");
            write_density_binary(&calc_output, bin_path.to_str().unwrap())
                .map_err(|e| PyRuntimeError::new_err(format!("Write binary failed: {}", e)))?;
            println!("  ✓ Binary density data:      {}", bin_path.display());

            let mol_json_str = molecule_to_json(&mol, &calc_output);
            let density_b64 = density_to_base64(&calc_output);
            let html = get_html_template(&mol_json_str, &density_b64, &svg);
            let html_path = path.join("visualization.html");
            fs::write(&html_path, html).ok();
            println!("  ✓ HTML 3D visualization:    {}", html_path.display());

            println!("  Completed in {:.2?}", t3.elapsed());
        }

        println!("\n╔══════════════════════════════════════════════╗");
        println!("║          Calculation Complete!               ║");
        println!("╚══════════════════════════════════════════════╝");

        Python::with_gil(|py| {
            let dict = pyo3::types::PyDict::new(py);
            dict.set_item("eigenvalues", &result.eigenvalues)?;
            dict.set_item("homo_lumo_gap", gap)?;
            dict.set_item("num_atoms", n)?;
            dict.set_item("num_electrons", num_electrons)?;
            dict.set_item("grid_dims", grid_dims.to_vec())?;
            dict.set_item("grid_origin", grid_origin.to_vec())?;
            dict.set_item("grid_spacing", grid_spacing)?;
            Ok(dict.to_object(py))
        })
    }

    fn __repr__(&self) -> PyResult<String> {
        Ok(format!(
            "HuckelEngine(use_gpu={}, gpu_initialized={})",
            self.use_gpu,
            self.gpu_solver.is_some()
        ))
    }
}

#[pyfunction]
fn load_molecule_json(path: &str) -> PyResult<String> {
    let content = fs::read_to_string(path)
        .map_err(|e| PyRuntimeError::new_err(format!("Cannot read file: {}", e)))?;
    Ok(content)
}

#[pyfunction]
fn run_calculation(
    mol_json: &str,
    output_dir: &str,
    use_gpu: bool,
    grid_resolution: usize,
    grid_padding: f64,
) -> PyResult<PyObject> {
    Python::with_gil(|py| {
        let mut engine = HuckelEngine::new(use_gpu);
        engine.calculate(py, mol_json, Some(output_dir), grid_resolution, grid_padding)
    })
}

#[pymodule]
fn huckel_engine(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<HuckelEngine>()?;
    m.add_function(wrap_pyfunction!(load_molecule_json, m)?)?;
    m.add_function(wrap_pyfunction!(run_calculation, m)?)?;
    Ok(())
}
