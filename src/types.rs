use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Atom {
    pub element: String,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Bond {
    pub atom1: usize,
    pub atom2: usize,
    pub bond_order: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub charge: Option<i32>,
    pub multiplicity: Option<i32>,
}

#[derive(Debug, Clone)]
pub struct EigenResult {
    pub eigenvalues: Vec<f64>,
    pub eigenvectors: Vec<Vec<f64>>,
}

#[derive(Debug, Clone, Serialize)]
pub struct CalculationOutput {
    pub eigenvalues: Vec<f64>,
    pub homo_lumo_gap: f64,
    pub electron_density_grid: Vec<f32>,
    pub grid_dims: [usize; 3],
    pub grid_origin: [f64; 3],
    pub grid_spacing: f64,
    pub num_electrons: usize,
}

impl Molecule {
    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }

    pub fn num_pi_electrons(&self) -> usize {
        let mut count = 0;
        for atom in &self.atoms {
            let pi_e = match atom.element.as_str() {
                "C" => 1,
                "N" => 1,
                "O" => 1,
                "S" => 1,
                "P" => 1,
                "B" => 0,
                _ => 1,
            };
            count += pi_e;
        }
        if let Some(charge) = self.charge {
            count = (count as i32 - charge) as usize;
        }
        count
    }

    pub fn bounding_box(&self) -> ([f64; 3], [f64; 3]) {
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for atom in &self.atoms {
            min[0] = min[0].min(atom.x);
            min[1] = min[1].min(atom.y);
            min[2] = min[2].min(atom.z);
            max[0] = max[0].max(atom.x);
            max[1] = max[1].max(atom.y);
            max[2] = max[2].max(atom.z);
        }
        (min, max)
    }

    pub fn atom_radius(element: &str) -> f64 {
        match element {
            "H" => 0.31,
            "C" => 0.76,
            "N" => 0.71,
            "O" => 0.66,
            "F" => 0.57,
            "S" => 1.05,
            "P" => 1.07,
            "Cl" => 1.02,
            "B" => 0.84,
            _ => 1.0,
        }
    }
}

pub fn gaussian_3d(x: f64, y: f64, z: f64, cx: f64, cy: f64, cz: f64, sigma: f64) -> f64 {
    let dx = x - cx;
    let dy = y - cy;
    let dz = z - cz;
    let r2 = dx * dx + dy * dy + dz * dz;
    let norm = 1.0 / ((2.0 * PI) * sigma * sigma).sqrt();
    norm * (-r2 / (2.0 * sigma * sigma)).exp()
}
