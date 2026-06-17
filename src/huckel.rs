use crate::types::{Atom, Bond, Molecule};

const ALPHA_C: f64 = -11.2;
const ALPHA_N: f64 = -13.9;
const ALPHA_O: f64 = -15.8;
const ALPHA_S: f64 = -11.6;
const ALPHA_P: f64 = -10.0;
const ALPHA_B: f64 = -8.5;

const BETA_CC: f64 = -0.815;
const BETA_CN: f64 = -0.85;
const BETA_CO: f64 = -0.90;
const BETA_CS: f64 = -0.70;
const BETA_CP: f64 = -0.70;
const BETA_NN: f64 = -0.80;
const BETA_NO: f64 = -0.92;
const BETA_OO: f64 = -0.95;

fn get_alpha(atom: &Atom) -> f64 {
    match atom.element.as_str() {
        "C" => ALPHA_C,
        "N" => ALPHA_N,
        "O" => ALPHA_O,
        "S" => ALPHA_S,
        "P" => ALPHA_P,
        "B" => ALPHA_B,
        _ => ALPHA_C,
    }
}

fn get_beta(a1: &Atom, a2: &Atom, bond_order: f64) -> f64 {
    let base = match (a1.element.as_str(), a2.element.as_str()) {
        ("C", "C") => BETA_CC,
        ("C", "N") | ("N", "C") => BETA_CN,
        ("C", "O") | ("O", "C") => BETA_CO,
        ("C", "S") | ("S", "C") => BETA_CS,
        ("C", "P") | ("P", "C") => BETA_CP,
        ("N", "N") => BETA_NN,
        ("N", "O") | ("O", "N") => BETA_NO,
        ("O", "O") => BETA_OO,
        _ => BETA_CC,
    };
    let k = match bond_order {
        b if b < 1.1 => 1.0,
        b if b < 1.6 => 1.1,
        b if b < 2.1 => 1.2,
        _ => 1.3,
    };
    base * k
}

fn bond_distance(a1: &Atom, a2: &Atom) -> f64 {
    let dx = a1.x - a2.x;
    let dy = a1.y - a2.y;
    let dz = a1.z - a2.z;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn infer_bonds(mol: &Molecule) -> Vec<Bond> {
    let mut bonds = mol.bonds.clone();
    if !bonds.is_empty() {
        return bonds;
    }
    let n = mol.atoms.len();
    for i in 0..n {
        for j in (i + 1)..n {
            let r = bond_distance(&mol.atoms[i], &mol.atoms[j]);
            let r_cov = Molecule::atom_radius(&mol.atoms[i].element)
                + Molecule::atom_radius(&mol.atoms[j].element);
            if r < r_cov * 1.2 {
                bonds.push(Bond {
                    atom1: i,
                    atom2: j,
                    bond_order: 1.0,
                });
            }
        }
    }
    bonds
}

pub fn build_hamiltonian(mol: &Molecule) -> Vec<Vec<f64>> {
    let n = mol.num_atoms();
    let mut h = vec![vec![0.0; n]; n];

    for i in 0..n {
        h[i][i] = get_alpha(&mol.atoms[i]);
    }

    let bonds = infer_bonds(mol);
    for bond in &bonds {
        let i = bond.atom1;
        let j = bond.atom2;
        if i < n && j < n {
            let beta = get_beta(&mol.atoms[i], &mol.atoms[j], bond.bond_order);
            h[i][j] = beta;
            h[j][i] = beta;
        }
    }
    h
}

pub fn build_overlap_matrix(mol: &Molecule) -> Vec<Vec<f64>> {
    let n = mol.num_atoms();
    let mut s = vec![vec![0.0; n]; n];
    let bonds = infer_bonds(mol);

    for i in 0..n {
        s[i][i] = 1.0;
    }
    for bond in &bonds {
        let i = bond.atom1;
        let j = bond.atom2;
        if i < n && j < n {
            let r = bond_distance(&mol.atoms[i], &mol.atoms[j]);
            let a = 2.0;
            let overlap = (1.0 + a * r / 2.0) * (-a * r).exp();
            s[i][j] = overlap;
            s[j][i] = overlap;
        }
    }
    s
}

pub fn compute_electron_density(
    mol: &Molecule,
    eigenvectors: &Vec<Vec<f64>>,
    num_electrons: usize,
    grid_origin: [f64; 3],
    grid_dims: [usize; 3],
    grid_spacing: f64,
) -> Vec<f32> {
    let n = mol.num_atoms();
    let num_orbitals = num_electrons / 2;
    let total_voxels = grid_dims[0] * grid_dims[1] * grid_dims[2];
    let mut density = vec![0.0f32; total_voxels];

    for iz in 0..grid_dims[2] {
        let z = grid_origin[2] + iz as f64 * grid_spacing;
        for iy in 0..grid_dims[1] {
            let y = grid_origin[1] + iy as f64 * grid_spacing;
            for ix in 0..grid_dims[0] {
                let x = grid_origin[0] + ix as f64 * grid_spacing;
                let idx = iz * grid_dims[1] * grid_dims[0] + iy * grid_dims[0] + ix;

                let mut ao_values = vec![0.0f64; n];
                for (a_idx, atom) in mol.atoms.iter().enumerate() {
                    let sigma = Molecule::atom_radius(&atom.element) * 0.8;
                    ao_values[a_idx] = gaussian(x, y, z, atom.x, atom.y, atom.z, sigma);
                }

                let mut rho = 0.0;
                for orb in 0..num_orbitals {
                    let mut psi = 0.0;
                    for a in 0..n {
                        psi += eigenvectors[orb][a] * ao_values[a];
                    }
                    rho += 2.0 * psi * psi;
                }
                density[idx] = rho as f32;
            }
        }
    }
    density
}

fn gaussian(x: f64, y: f64, z: f64, cx: f64, cy: f64, cz: f64, sigma: f64) -> f64 {
    let dx = x - cx;
    let dy = y - cy;
    let dz = z - cz;
    let r2 = dx * dx + dy * dy + dz * dz;
    let pi_sq = std::f64::consts::PI.powi(3);
    let norm = 1.0 / (pi_sq * sigma.powi(6)).sqrt();
    norm * (-r2 / (sigma * sigma)).exp()
}
