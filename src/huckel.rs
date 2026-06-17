use crate::types::{Atom, Bond, Molecule, Algorithm, get_ehmo_params, total_ehmo_basis_size, atom_orbital_offset, num_valence_orbitals, is_heavy_metal};

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

const WOLFSBERG_HELMHOLZ_K: f64 = 1.75;

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum OrbitalType {
    S, Px, Py, Pz, Dxy, Dyz, Dzx, Dx2y2, Dz2,
}

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

fn direction_cosines(a1: &Atom, a2: &Atom) -> (f64, f64, f64, f64) {
    let r = bond_distance(a1, a2);
    if r < 1e-12 { return (0.0, 0.0, 0.0, 1e-12); }
    let l = (a2.x - a1.x) / r;
    let m = (a2.y - a1.y) / r;
    let nn = (a2.z - a1.z) / r;
    (l, m, nn, r)
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

fn slater_s_s_overlap(zeta_a: f64, zeta_b: f64, r: f64) -> f64 {
    if r < 1e-12 { return 1.0; }
    let p = 0.5 * (zeta_a + zeta_b) * r;
    let e = (-2.0 * p).exp();
    let s_s = 1.0 + p + p * p / 3.0;
    s_s * e
}

fn slater_p_sigma_overlap(zeta_a: f64, zeta_b: f64, r: f64) -> f64 {
    if r < 1e-12 { return 0.0; }
    let p = 0.5 * (zeta_a + zeta_b) * r;
    let e = (-2.0 * p).exp();
    let sqrt3 = 3.0_f64.sqrt();
    -(p / sqrt3) * e
}

fn slater_p_pi_overlap(zeta_a: f64, zeta_b: f64, r: f64) -> f64 {
    if r < 1e-12 { return 0.0; }
    let p = 0.5 * (zeta_a + zeta_b) * r;
    let e = (-2.0 * p).exp();
    (p * (1.0 + p * 0.3)) * e * 0.8
}

fn slater_s_d_overlap(zeta_a: f64, zeta_b: f64, r: f64) -> f64 {
    if r < 1e-12 { return 0.0; }
    let p = 0.5 * (zeta_a + zeta_b) * r;
    let e = (-2.0 * p).exp();
    0.1 * p * p * e
}

fn slater_p_d_overlap(zeta_a: f64, zeta_b: f64, r: f64) -> f64 {
    if r < 1e-12 { return 0.0; }
    let p = 0.5 * (zeta_a + zeta_b) * r;
    let e = (-2.0 * p).exp();
    0.15 * p * p * e
}

fn slater_d_sigma_overlap(zeta_a: f64, zeta_b: f64, r: f64) -> f64 {
    if r < 1e-12 { return 1.0; }
    let p = 0.5 * (zeta_a + zeta_b) * r;
    let e = (-2.0 * p).exp();
    (1.0 + p * 0.5) * p * p * e * 0.3
}

fn slater_d_pi_overlap(zeta_a: f64, zeta_b: f64, r: f64) -> f64 {
    if r < 1e-12 { return 0.0; }
    let p = 0.5 * (zeta_a + zeta_b) * r;
    let e = (-2.0 * p).exp();
    0.2 * p * p * e
}

fn slater_d_delta_overlap(zeta_a: f64, zeta_b: f64, r: f64) -> f64 {
    if r < 1e-12 { return 0.0; }
    let p = 0.5 * (zeta_a + zeta_b) * r;
    let e = (-2.0 * p).exp();
    0.1 * p * p * e
}

fn orbital_principal_number(element: &str, orbital: OrbitalType) -> u32 {
    use OrbitalType::*;
    let z = crate::types::get_atomic_number(element);
    let has_3d = z >= 21;
    let has_4d = z >= 39;
    let has_5d = z >= 57;
    match orbital {
        S | Px | Py | Pz => {
            if z <= 2 { 1 }
            else if z <= 10 { 2 }
            else if z <= 18 { 3 }
            else if z <= 36 { 4 }
            else if z <= 54 { 5 }
            else { 6 }
        }
        Dxy | Dyz | Dzx | Dx2y2 | Dz2 => {
            if has_5d { 5 }
            else if has_4d { 4 }
            else if has_3d { 3 }
            else { 3 }
        }
    }
}

fn orbital_hii_for(element: &str, orbital: OrbitalType) -> f64 {
    use OrbitalType::*;
    let params = get_ehmo_params(element);
    match orbital {
        S => params.hii_s,
        Px | Py | Pz => params.hii_p,
        Dxy | Dyz | Dzx | Dx2y2 | Dz2 => params.hii_d,
    }
}

fn orbital_zeta_for(element: &str, orbital: OrbitalType) -> f64 {
    use OrbitalType::*;
    let params = get_ehmo_params(element);
    match orbital {
        S => params.zeta_s,
        Px | Py | Pz => params.zeta_p,
        Dxy | Dyz | Dzx | Dx2y2 | Dz2 => params.zeta_d,
    }
}

fn atom_orbital_type_list(orbitals_per_atom: usize) -> Vec<OrbitalType> {
    use OrbitalType::*;
    let all = vec![S, Px, Py, Pz, Dxy, Dyz, Dzx, Dx2y2, Dz2];
    all.into_iter().take(orbitals_per_atom).collect()
}

fn compute_orbital_overlap(
    a1: &Atom, a2: &Atom,
    orb1: OrbitalType, orb2: OrbitalType,
    zeta1: f64, zeta2: f64,
) -> f64 {
    use OrbitalType::*;
    let (l, m, nn, r) = direction_cosines(a1, a2);

    if r < 1e-12 {
        return if orb1 == orb2 && std::ptr::eq(a1, a2) { 1.0 } else { 0.0 };
    }

    let s_sigma = slater_s_s_overlap(zeta1, zeta2, r);
    let p_sigma_12 = slater_p_sigma_overlap(zeta1, zeta2, r);
    let p_sigma_21 = slater_p_sigma_overlap(zeta2, zeta1, r);
    let p_pi = slater_p_pi_overlap(zeta1, zeta2, r);
    let s_d_sigma = slater_s_d_overlap(zeta1, zeta2, r);
    let p_d_sigma = slater_p_d_overlap(zeta1, zeta2, r);
    let d_sigma = slater_d_sigma_overlap(zeta1, zeta2, r);
    let d_pi = slater_d_pi_overlap(zeta1, zeta2, r);
    let d_delta = slater_d_delta_overlap(zeta1, zeta2, r);

    let pp_sigma = |z1, z2, rr| {
        let ppp = 0.5 * (z1 + z2) * rr;
        (-2.0 * ppp).exp() * (1.0 + ppp + ppp * ppp / 5.0) * 0.9
    };
    let pp_sigma_val = pp_sigma(zeta1, zeta2, r);

    match (orb1, orb2) {
        (S, S) => s_sigma,
        (S, Px) => l * p_sigma_21,
        (S, Py) => m * p_sigma_21,
        (S, Pz) => nn * p_sigma_21,
        (Px, S) => -l * p_sigma_12,
        (Py, S) => -m * p_sigma_12,
        (Pz, S) => -nn * p_sigma_12,

        (Px, Px) => l * l * pp_sigma_val + (1.0 - l * l) * p_pi,
        (Px, Py) => l * m * (pp_sigma_val - p_pi),
        (Px, Pz) => l * nn * (pp_sigma_val - p_pi),
        (Py, Px) => m * l * (pp_sigma_val - p_pi),
        (Py, Py) => m * m * pp_sigma_val + (1.0 - m * m) * p_pi,
        (Py, Pz) => m * nn * (pp_sigma_val - p_pi),
        (Pz, Px) => nn * l * (pp_sigma_val - p_pi),
        (Pz, Py) => nn * m * (pp_sigma_val - p_pi),
        (Pz, Pz) => nn * nn * pp_sigma_val + (1.0 - nn * nn) * p_pi,

        (S, Dx2y2) | (Dx2y2, S) => (l * l - m * m) * s_d_sigma,
        (S, Dz2) | (Dz2, S) => (3.0 * nn * nn - 1.0) * s_d_sigma / 3.0_f64.sqrt(),
        (S, Dxy) | (Dxy, S) => 2.0 * l * m * s_d_sigma,
        (S, Dyz) | (Dyz, S) => 2.0 * m * nn * s_d_sigma,
        (S, Dzx) | (Dzx, S) => 2.0 * nn * l * s_d_sigma,

        (Px, Dx2y2) | (Dx2y2, Px) => l * (l * l - 3.0 * m * m) * p_d_sigma,
        (Px, Dz2) | (Dz2, Px) => l * (3.0 * nn * nn - 1.0) * p_d_sigma / 3.0_f64.sqrt(),
        (Px, Dxy) | (Dxy, Px) => 3.0 * l * l * m * p_d_sigma,
        (Px, Dyz) | (Dyz, Px) => 2.0 * l * m * nn * p_d_sigma,
        (Px, Dzx) | (Dzx, Px) => 2.0 * l * l * nn * p_d_sigma,

        (Py, Dx2y2) | (Dx2y2, Py) => m * (3.0 * l * l - m * m) * p_d_sigma,
        (Py, Dz2) | (Dz2, Py) => m * (3.0 * nn * nn - 1.0) * p_d_sigma / 3.0_f64.sqrt(),
        (Py, Dxy) | (Dxy, Py) => 3.0 * l * m * m * p_d_sigma,
        (Py, Dyz) | (Dyz, Py) => 2.0 * m * m * nn * p_d_sigma,
        (Py, Dzx) | (Dzx, Py) => 2.0 * l * m * nn * p_d_sigma,

        (Pz, Dx2y2) | (Dx2y2, Pz) => nn * (l * l - m * m) * p_d_sigma,
        (Pz, Dz2) | (Dz2, Pz) => nn * (3.0 * nn * nn - 5.0) * p_d_sigma / 3.0_f64.sqrt(),
        (Pz, Dxy) | (Dxy, Pz) => 2.0 * l * m * nn * p_d_sigma,
        (Pz, Dyz) | (Dyz, Pz) => m * (3.0 * nn * nn - 1.0) * p_d_sigma,
        (Pz, Dzx) | (Dzx, Pz) => l * (3.0 * nn * nn - 1.0) * p_d_sigma,

        (Dx2y2, Dx2y2) => d_d_mat(l, m, nn, d_sigma, d_pi, d_delta, "xxyy_xxyy"),
        (Dz2, Dz2) => d_d_mat(l, m, nn, d_sigma, d_pi, d_delta, "zz_zz"),
        (Dxy, Dxy) => d_d_mat(l, m, nn, d_sigma, d_pi, d_delta, "xy_xy"),
        (Dyz, Dyz) => d_d_mat(l, m, nn, d_sigma, d_pi, d_delta, "yz_yz"),
        (Dzx, Dzx) => d_d_mat(l, m, nn, d_sigma, d_pi, d_delta, "zx_zx"),

        (Dx2y2, Dz2) | (Dz2, Dx2y2) => 3.0_f64.sqrt() * nn * nn * (l * l - m * m) * (d_sigma - d_pi) * 0.5,
        (Dx2y2, Dxy) | (Dxy, Dx2y2) => l * m * (l * l - m * m) * (d_sigma - d_delta),
        (Dx2y2, Dyz) | (Dyz, Dx2y2) => m * nn * (l * l - m * m) * d_pi,
        (Dx2y2, Dzx) | (Dzx, Dx2y2) => l * nn * (l * l - m * m) * d_pi,

        (Dz2, Dxy) | (Dxy, Dz2) => 3.0_f64.sqrt() * l * m * nn * nn * (d_sigma - d_pi),
        (Dz2, Dyz) | (Dyz, Dz2) => 3.0_f64.sqrt() * m * nn * (3.0 * nn * nn - 1.0) * (d_sigma - d_pi) / 3.0,
        (Dz2, Dzx) | (Dzx, Dz2) => 3.0_f64.sqrt() * l * nn * (3.0 * nn * nn - 1.0) * (d_sigma - d_pi) / 3.0,

        (Dxy, Dyz) | (Dyz, Dxy) => 2.0 * l * m * m * nn * d_delta,
        (Dxy, Dzx) | (Dzx, Dxy) => 2.0 * l * l * m * nn * d_delta,
        (Dyz, Dzx) | (Dzx, Dyz) => l * m * nn * nn * d_delta,

        _ => 0.0,
    }
}

fn d_d_mat(l: f64, m: f64, nn: f64, s: f64, p: f64, d: f64, mode: &str) -> f64 {
    let l2 = l * l; let m2 = m * m; let n2 = nn * nn;
    let l4 = l2 * l2; let m4 = m2 * m2; let n4 = n2 * n2;
    let lm = l * m;
    match mode {
        "xxyy_xxyy" => {
            (l4 - 6.0 * l2 * m2 + m4) * s +
            4.0 * (l2 * n2 + m2 * n2 - lm * lm) * p +
            (12.0 * lm * lm - 4.0 * l2 * m2) * d
        }
        "zz_zz" => {
            (3.0 * n2 - 1.0).powi(2) * s / 3.0 +
            4.0 * n2 * (1.0 - n2).powi(2) * p +
            (1.0 - n2).powi(2) * d
        }
        "xy_xy" => {
            4.0 * l2 * m2 * s +
            (l2 + m2 - 4.0 * l2 * m2) * n2 * p +
            (l4 + m4 - 4.0 * lm * lm) * d
        }
        "yz_yz" => {
            4.0 * m2 * n2 * s +
            (m2 + n2 - 4.0 * m2 * n2) * l2 * p +
            (m4 + n4 - 4.0 * m2 * n2) * d
        }
        "zx_zx" => {
            4.0 * n2 * l2 * s +
            (n2 + l2 - 4.0 * n2 * l2) * m2 * p +
            (n4 + l4 - 4.0 * n2 * l2) * d
        }
        _ => 0.0,
    }
}

fn wolfsberg_helmholz(hii: f64, hjj: f64, sij: f64) -> f64 {
    if sij.abs() < 1e-10 { return 0.0; }
    WOLFSBERG_HELMHOLZ_K * sij * (hii + hjj) * 0.5
}

pub fn build_ehmo_hamiltonian(mol: &Molecule) -> Vec<Vec<f64>> {
    let total = total_ehmo_basis_size(mol);
    let mut h = vec![vec![0.0f64; total]; total];

    let atoms = &mol.atoms;
    for (a1_idx, atom1) in atoms.iter().enumerate() {
        let n1 = num_valence_orbitals(&atom1.element);
        let offset1 = atom_orbital_offset(mol, a1_idx);
        let orb_list1 = atom_orbital_type_list(n1);
        for (local1, orb1) in orb_list1.iter().enumerate() {
            let i = offset1 + local1;
            let hii = orbital_hii_for(&atom1.element, *orb1);
            h[i][i] = hii;

            for (a2_idx, atom2) in atoms.iter().enumerate() {
                if a2_idx == a1_idx { continue; }
                let n2 = num_valence_orbitals(&atom2.element);
                let offset2 = atom_orbital_offset(mol, a2_idx);
                let orb_list2 = atom_orbital_type_list(n2);
                for (local2, orb2) in orb_list2.iter().enumerate() {
                    let j = offset2 + local2;
                    if j < total && i < total {
                        let zeta1 = orbital_zeta_for(&atom1.element, *orb1);
                        let zeta2 = orbital_zeta_for(&atom2.element, *orb2);
                        let sij = compute_orbital_overlap(atom1, atom2, *orb1, *orb2, zeta1, zeta2);
                        let hjj = orbital_hii_for(&atom2.element, *orb2);
                        h[i][j] = wolfsberg_helmholz(hii, hjj, sij);
                    }
                }
            }
        }
    }
    h
}

pub fn build_ehmo_overlap(mol: &Molecule) -> Vec<Vec<f64>> {
    let total = total_ehmo_basis_size(mol);
    let mut s = vec![vec![0.0f64; total]; total];

    let atoms = &mol.atoms;
    for (a1_idx, atom1) in atoms.iter().enumerate() {
        let n1 = num_valence_orbitals(&atom1.element);
        let offset1 = atom_orbital_offset(mol, a1_idx);
        let orb_list1 = atom_orbital_type_list(n1);
        for (local1, orb1) in orb_list1.iter().enumerate() {
            let i = offset1 + local1;
            if i < total { s[i][i] = 1.0; }

            for (a2_idx, atom2) in atoms.iter().enumerate().skip(a1_idx + 1) {
                let n2 = num_valence_orbitals(&atom2.element);
                let offset2 = atom_orbital_offset(mol, a2_idx);
                let orb_list2 = atom_orbital_type_list(n2);
                for (local2, orb2) in orb_list2.iter().enumerate() {
                    let j = offset2 + local2;
                    if i < total && j < total {
                        let zeta1 = orbital_zeta_for(&atom1.element, *orb1);
                        let zeta2 = orbital_zeta_for(&atom2.element, *orb2);
                        let sij = compute_orbital_overlap(atom1, atom2, *orb1, *orb2, zeta1, zeta2);
                        s[i][j] = sij;
                        s[j][i] = sij;
                    }
                }
            }
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
    algorithm: Algorithm,
) -> Vec<f32> {
    let (n_basis, total_voxels, is_ehmo) = match algorithm {
        Algorithm::ExtendedHuckel => {
            let t = total_ehmo_basis_size(mol);
            (t, grid_dims[0] * grid_dims[1] * grid_dims[2], true)
        }
        _ => (mol.num_atoms(), grid_dims[0] * grid_dims[1] * grid_dims[2], false),
    };
    let num_orbitals = (num_electrons / 2).max(1).min(eigenvectors.len());
    let mut density = vec![0.0f32; total_voxels];

    for iz in 0..grid_dims[2] {
        let z = grid_origin[2] + iz as f64 * grid_spacing;
        for iy in 0..grid_dims[1] {
            let y = grid_origin[1] + iy as f64 * grid_spacing;
            for ix in 0..grid_dims[0] {
                let x = grid_origin[0] + ix as f64 * grid_spacing;
                let idx = iz * grid_dims[1] * grid_dims[0] + iy * grid_dims[0] + ix;

                let ao_values = if is_ehmo {
                    eval_ehmo_ao_grid(mol, x, y, z, n_basis)
                } else {
                    let mut v = vec![0.0f64; n_basis];
                    for (a_idx, atom) in mol.atoms.iter().enumerate() {
                        let sigma = Molecule::atom_radius(&atom.element) * 0.8;
                        v[a_idx] = gaussian(x, y, z, atom.x, atom.y, atom.z, sigma);
                    }
                    v
                };

                let mut rho = 0.0;
                for orb in 0..num_orbitals {
                    let mut psi = 0.0;
                    for a in 0..n_basis {
                        if a < eigenvectors[orb].len() && a < ao_values.len() {
                            psi += eigenvectors[orb][a] * ao_values[a];
                        }
                    }
                    rho += 2.0 * psi * psi;
                }
                density[idx] = rho as f32;
            }
        }
    }
    density
}

fn eval_ehmo_ao_grid(mol: &Molecule, x: f64, y: f64, z: f64, n_basis: usize) -> Vec<f64> {
    use OrbitalType::*;
    let mut vals = vec![0.0f64; n_basis];

    for (a_idx, atom) in mol.atoms.iter().enumerate() {
        let n_orb = num_valence_orbitals(&atom.element);
        let offset = atom_orbital_offset(mol, a_idx);
        let orbs = atom_orbital_type_list(n_orb);

        let dx = x - atom.x;
        let dy = y - atom.y;
        let dz = z - atom.z;
        let r = (dx * dx + dy * dy + dz * dz).sqrt();

        for (local, orb) in orbs.iter().enumerate() {
            let idx = offset + local;
            if idx >= n_basis { continue; }
            let zeta = orbital_zeta_for(&atom.element, *orb);
            let n_prin = orbital_principal_number(&atom.element, *orb) as f64;

            let radial = if r < 20.0 / zeta.max(0.1) {
                let z2z = 2.0 * zeta;
                let norm = z2z.powf(n_prin + 0.5) / (2.0 * (2.0 * n_prin).sqrt() * n_prin.powf(n_prin as i32));
                r.powf(n_prin - 1.0) * (-zeta * r).exp() * norm
            } else {
                0.0
            };

            let angular: f64 = match orb {
                S => 1.0 / std::f64::consts::PI.sqrt(),
                Px => dx,
                Py => dy,
                Pz => dz,
                Dxy => dx * dy,
                Dyz => dy * dz,
                Dzx => dz * dx,
                Dx2y2 => dx * dx - dy * dy,
                Dz2 => 2.0 * dz * dz - dx * dx - dy * dy,
            };

            vals[idx] = radial * angular;
        }
    }
    vals
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
