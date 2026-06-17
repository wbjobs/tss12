use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Algorithm {
    SimpleHuckel,
    ExtendedHuckel,
}

impl Default for Algorithm {
    fn default() -> Self { Algorithm::SimpleHuckel }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Atom {
    pub element: String,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Debug, Clone, Copy)]
pub struct EHMOParams {
    pub atomic_number: u32,
    pub is_transition_metal: bool,
    pub valence_e: u32,
    pub hs_1: f64,
    pub hs_2s: f64,
    pub hs_2p: f64,
    pub hs_3d: f64,
    pub zeta_s: f64,
    pub zeta_p: f64,
    pub zeta_d: f64,
}

impl EHMOParams {
    pub fn num_valence_orbitals(&self) -> usize {
        if self.is_transition_metal { 9 } else if self.atomic_number > 2 { 4 } else { 1 }
    }
}

pub fn get_atomic_number(element: &str) -> u32 {
    match element {
        "H" => 1, "He" => 2, "Li" => 3, "Be" => 4, "B" => 5, "C" => 6,
        "N" => 7, "O" => 8, "F" => 9, "Ne" => 10, "Na" => 11, "Mg" => 12,
        "Al" => 13, "Si" => 14, "P" => 3, "S" => 16, "Cl" => 17, "Ar" => 18,
        "K" => 19, "Ca" => 20, "Sc" => 21, "Ti" => 22, "V" => 23, "Cr" => 24,
        "Mn" => 25, "Fe" => 26, "Co" => 27, "Ni" => 28, "Cu" => 29, "Zn" => 30,
        "Ga" => 31, "Ge" => 32, "As" => 33, "Se" => 34, "Br" => 35, "Kr" => 36,
        "Rb" => 37, "Sr" => 38, "Y" => 39, "Zr" => 40, "Nb" => 41, "Mo" => 42,
        "Tc" => 43, "Ru" => 44, "Rh" => 45, "Pd" => 46, "Ag" => 47, "Cd" => 48,
        "In" => 49, "Sn" => 50, "Sb" => 51, "Te" => 52, "I" => 53, "Xe" => 54,
        "Cs" => 55, "Ba" => 56, "La" => 57,
        "Ce" => 58, "Pr" => 59, "Nd" => 60, "Pm" => 61, "Sm" => 62, "Eu" => 63,
        "Gd" => 64, "Tb" => 65, "Dy" => 66, "Ho" => 67, "Er" => 68, "Tm" => 69,
        "Yb" => 70, "Lu" => 71,
        "Hf" => 72, "Ta" => 73, "W" => 74, "Re" => 75, "Os" => 76, "Ir" => 77,
        "Pt" => 78, "Au" => 79, "Hg" => 80, "Tl" => 81, "Pb" => 82, "Bi" => 83,
        _ => 0,
    }
}

pub fn is_heavy_metal(element: &str) -> bool {
    let z = get_atomic_number(element);
    if z == 0 { return false; }
    matches!(element,
        "Sc"|"Ti"|"V"|"Cr"|"Mn"|"Fe"|"Co"|"Ni"|"Cu"|"Zn"|
        "Y"|"Zr"|"Nb"|"Mo"|"Tc"|"Ru"|"Rh"|"Pd"|"Ag"|"Cd"|
        "Hf"|"Ta"|"W"|"Re"|"Os"|"Ir"|"Pt"|"Au"|"Hg"|
        "La"|"Ce"|"Pr"|"Nd"|"Pm"|"Sm"|"Eu"|"Gd"|"Tb"|"Dy"|
        "Ho"|"Er"|"Tm"|"Yb"|"Lu"
    ) || z >= 78
}

pub fn contains_heavy_metal(mol: &Molecule) -> bool {
    mol.atoms.iter().any(|a| is_heavy_metal(&a.element))
}

pub fn auto_select_algorithm(mol: &Molecule) -> Algorithm {
    if contains_heavy_metal(mol) {
        Algorithm::ExtendedHuckel
    } else {
        Algorithm::SimpleHuckel
    }
}

pub fn get_ehmo_params(element: &str) -> Option<EHMOParams> {
    let z = get_atomic_number(element);
    if z == 0 { return None; }
    let tm = is_heavy_metal(element);
    let (ve, h1s, h2s, h2p, h3d, zs, zp, zd) = match element {
        "H" => (1, -13.6, 0.0, 0.0, 0.0, 1.24, 0.0, 0.0),
        "C" => (4, 0.0, -21.4, -11.4, 0.0, 1.625, 1.625, 0.0),
        "N" => (5, 0.0, -26.0, -13.4, 0.0, 1.95, 1.95, 0.0),
        "O" => (6, 0.0, -32.3, -15.8, 0.0, 2.275, 2.275, 0.0),
        "F" => (7, 0.0, -40.0, -18.1, 0.0, 2.55, 2.55, 0.0),
        "B" => (3, 0.0, -15.2, -8.5, 0.0, 1.30, 1.30, 0.0),
        "S" => (6, 0.0, -20.0, -13.3, -8.0, 2.122, 1.827, 1.383),
        "P" => (5, 0.0, -18.6, -10.0, -7.0, 1.88, 1.55, 1.30),
        "Cl" => (7, 0.0, -24.5, -15.0, -9.0, 2.36, 2.04, 1.42),
        "Si" => (4, 0.0, -17.3, -8.3, -6.0, 1.64, 1.39, 1.30),
        "Fe" => (8, 0.0, 0.0, 0.0, -12.6, 2.20, 2.20, 4.65),
        "Cu" => (11, 0.0, 0.0, 0.0, -14.0, 2.40, 2.40, 4.70),
        "Zn" => (12, 0.0, 0.0, 0.0, -10.8, 2.45, 2.45, 4.80),
        "Pt" => (10, 0.0, 0.0, 0.0, -12.2, 2.50, 2.50, 5.76),
        "Au" => (11, 0.0, 0.0, 0.0, -14.0, 2.60, 2.60, 5.90),
        "Pd" => (10, 0.0, 0.0, 0.0, -12.1, 2.40, 2.40, 5.45),
        "Ag" => (11, 0.0, 0.0, 0.0, -14.0, 2.50, 2.50, 5.50),
        "Rh" => (9, 0.0, 0.0, 0.0, -12.1, 2.35, 2.35, 5.35),
        "Ru" => (8, 0.0, 0.0, 0.0, -12.1, 2.30, 2.30, 5.20),
        "Ir" => (9, 0.0, 0.0, 0.0, -12.4, 2.45, 2.45, 5.80),
        "Os" => (8, 0.0, 0.0, 0.0, -12.4, 2.40, 2.40, 5.70),
        "Ni" => (10, 0.0, 0.0, 0.0, -12.8, 2.25, 2.25, 4.75),
        "Co" => (9, 0.0, 0.0, 0.0, -12.6, 2.20, 2.20, 4.70),
        "Mn" => (7, 0.0, 0.0, 0.0, -11.4, 2.15, 2.15, 4.50),
        "Cr" => (6, 0.0, 0.0, 0.0, -11.2, 2.10, 2.10, 4.45),
        "V" => (5, 0.0, 0.0, 0.0, -10.1, 2.05, 2.05, 4.30),
        "Ti" => (4, 0.0, 0.0, 0.0, -9.2, 2.00, 2.00, 4.20),
        "W" => (6, 0.0, 0.0, 0.0, -11.2, 2.40, 2.40, 5.60),
        "Mo" => (6, 0.0, 0.0, 0.0, -10.5, 2.30, 2.30, 5.10),
        "Ta" => (5, 0.0, 0.0, 0.0, -9.5, 2.35, 2.35, 5.50),
        "Nb" => (5, 0.0, 0.0, 0.0, -9.8, 2.20, 2.20, 4.95),
        "Zr" => (4, 0.0, 0.0, 0.0, -8.8, 2.15, 2.15, 4.85),
        "Hf" => (4, 0.0, 0.0, 0.0, -9.0, 2.40, 2.40, 5.50),
        "Hg" => (12, 0.0, 0.0, 0.0, -10.4, 2.65, 2.65, 6.00),
        "Cd" => (12, 0.0, 0.0, 0.0, -10.0, 2.50, 2.50, 4.95),
        _ => (z.min(8) as u32, 0.0, -20.0, -12.0, -9.0, 2.0, 1.8, 4.5),
    };
    Some(EHMOParams {
        atomic_number: z,
        is_transition_metal: tm,
        valence_e: ve,
        hs_1: h1s, hs_2s: h2s, hs_2p: h2p, hs_3d: h3d,
        zeta_s: zs, zeta_p: zp, zeta_d: zd,
    })
}

pub fn num_valence_orbitals(element: &str) -> usize {
    if let Some(p) = get_ehmo_params(element) {
        p.num_valence_orbitals()
    } else {
        match get_atomic_number(element) {
            1 => 1,
            z if z <= 2 => 1,
            z if z <= 10 => 4,
            _ => 9,
        }
    }
}

pub fn total_ehmo_basis_size(mol: &Molecule) -> usize {
    mol.atoms.iter().map(|a| num_valence_orbitals(&a.element)).sum()
}

pub fn atom_orbital_offset(mol: &Molecule, atom_idx: usize) -> usize {
    mol.atoms.iter()
        .take(atom_idx)
        .map(|a| num_valence_orbitals(&a.element))
        .sum()
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
    pub algorithm: Algorithm,
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
    pub algorithm: String,
}

impl Molecule {
    pub fn num_atoms(&self) -> usize { self.atoms.len() }

    pub fn num_pi_electrons(&self) -> usize {
        let mut count = 0;
        for atom in &self.atoms {
            let pi_e = match atom.element.as_str() {
                "H" => 0,
                "C" => 1, "N" => 1, "O" => 1, "S" => 1, "P" => 1, "B" => 0,
                "F" => 2, "Cl" => 2, "Br" => 2, "I" => 2,
                e if is_heavy_metal(e) => {
                    get_ehmo_params(e).map(|p| p.valence_e as usize).unwrap_or(0)
                },
                _ => 1,
            };
            count += pi_e;
        }
        if let Some(charge) = self.charge {
            count = ((count as i32) - charge).max(0) as usize;
        }
        count
    }

    pub fn total_valence_electrons(&self) -> usize {
        let mut total: usize = self.atoms.iter()
            .map(|a| get_ehmo_params(&a.element).map(|p| p.valence_e as usize).unwrap_or(0))
            .sum();
        if let Some(charge) = self.charge {
            total = ((total as i32) - charge).max(0) as usize;
        }
        total
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
            "H" => 0.31, "He" => 0.28,
            "Li" => 1.28, "Be" => 0.96, "B" => 0.84, "C" => 0.76,
            "N" => 0.71, "O" => 0.66, "F" => 0.57, "Ne" => 0.58,
            "Na" => 1.66, "Mg" => 1.41, "Al" => 1.21, "Si" => 1.11,
            "P" => 1.07, "S" => 1.05, "Cl" => 1.02, "Ar" => 1.06,
            "K" => 2.03, "Ca" => 1.76, "Sc" => 1.70, "Ti" => 1.60,
            "V" => 1.53, "Cr" => 1.39, "Mn" => 1.39, "Fe" => 1.32,
            "Co" => 1.26, "Ni" => 1.24, "Cu" => 1.32, "Zn" => 1.22,
            "Ga" => 1.22, "Ge" => 1.20, "As" => 1.19, "Se" => 1.20,
            "Br" => 1.20, "Kr" => 1.16,
            "Rb" => 2.20, "Sr" => 1.95, "Y" => 1.90, "Zr" => 1.75,
            "Nb" => 1.64, "Mo" => 1.54, "Tc" => 1.47, "Ru" => 1.46,
            "Rh" => 1.42, "Pd" => 1.39, "Ag" => 1.45, "Cd" => 1.44,
            "In" => 1.42, "Sn" => 1.39, "Sb" => 1.39, "Te" => 1.38,
            "I" => 1.39, "Xe" => 1.40,
            "Cs" => 2.44, "Ba" => 2.15, "La" => 2.07, "Hf" => 1.67,
            "Ta" => 1.49, "W" => 1.41, "Re" => 1.37, "Os" => 1.35,
            "Ir" => 1.36, "Pt" => 1.39, "Au" => 1.42, "Hg" => 1.48,
            "Tl" => 1.45, "Pb" => 1.46, "Bi" => 1.48,
            _ => 1.5,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReactionInput {
    pub reactant: Molecule,
    pub product: Molecule,
    pub num_steps: Option<usize>,
    pub interpolation: Option<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct IntermediateState {
    pub step: usize,
    pub t: f64,
    pub molecule: Molecule,
    pub eigenvalues: Vec<f64>,
    pub total_energy: f64,
    pub homo_lumo_gap: f64,
    pub density: Vec<f32>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ReactionPathResult {
    pub steps: Vec<IntermediateState>,
    pub grid_dims: [usize; 3],
    pub grid_origin: [f64; 3],
    pub grid_spacing: f64,
    pub num_electrons: usize,
    pub algorithm: String,
    pub reactant_name: String,
    pub product_name: String,
}

pub fn gaussian_3d(x: f64, y: f64, z: f64, cx: f64, cy: f64, cz: f64, sigma: f64) -> f64 {
    let dx = x - cx; let dy = y - cy; let dz = z - cz;
    let r2 = dx * dx + dy * dy + dz * dz;
    let norm = 1.0 / ((2.0 * PI) * sigma * sigma).sqrt();
    norm * (-r2 / (2.0 * sigma * sigma)).exp()
}

pub fn slater_type_orbital_r(r: f64, zeta: f64, n: u32) -> f64 {
    let norm = ((2.0 * zeta).powi(2 * n as i32 + 1) / (2.0 * (2..=2*n).product::<u32>() as f64)).sqrt();
    norm * r.powi(n as i32 - 1) * (-zeta * r).exp()
}
