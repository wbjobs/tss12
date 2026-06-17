use crate::types::{EigenResult, Algorithm};

pub fn symmetrize(matrix: &mut Vec<Vec<f64>>) {
    let n = matrix.len();
    if n == 0 { return; }
    for i in 0..n {
        let m = matrix[i].len();
        for j in (i + 1)..m {
            if j < n && i < matrix[j].len() {
                let avg = 0.5 * (matrix[i][j] + matrix[j][i]);
                if avg.is_finite() {
                    matrix[i][j] = avg;
                    matrix[j][i] = avg;
                } else {
                    matrix[i][j] = 0.0;
                    matrix[j][i] = 0.0;
                }
            }
        }
    }
}

pub fn regularize_positive_definite(s: &mut Vec<Vec<f64>>, min_diag: f64) {
    let n = s.len();
    if n == 0 { return; }
    for i in 0..n {
        if i < s[i].len() {
            if s[i][i] <= 0.0 || !s[i][i].is_finite() {
                s[i][i] = min_diag;
            } else if s[i][i] < min_diag {
                s[i][i] = min_diag;
            }
        }
    }
    let eps = 1e-8;
    for i in 0..n {
        if i < s[i].len() {
            s[i][i] += eps;
        }
    }
}

fn ensure_square(matrix: &Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    let n = matrix.len();
    let mut m = matrix.clone();
    for row in &mut m {
        row.resize(n, 0.0);
    }
    m
}

fn sanitize_matrix(matrix: &mut Vec<Vec<f64>>) {
    let n = matrix.len();
    for i in 0..n {
        for j in 0..n {
            if !matrix[i][j].is_finite() {
                matrix[i][j] = 0.0;
            }
            let abs = matrix[i][j].abs();
            if abs > 1e10 {
                matrix[i][j] = matrix[i][j] * 1e-10;
            }
        }
    }
}

pub fn solve_eigen_jacobi(
    matrix: &Vec<Vec<f64>>,
    max_iter: usize,
    tolerance: f64,
) -> EigenResult {
    let mut a = ensure_square(matrix);
    sanitize_matrix(&mut a);
    symmetrize(&mut a);
    let n = a.len();

    let mut v = vec![vec![0.0f64; n]; n];
    for i in 0..n { v[i][i] = 1.0; }

    let mut iterations = 0;
    loop {
        let mut max_off = 0.0;
        let mut p = 0usize;
        let mut q = 1usize;
        for i in 0..n {
            for j in (i + 1)..n {
                let abs_val = a[i][j].abs();
                if abs_val > max_off {
                    max_off = abs_val;
                    p = i; q = j;
                }
            }
        }
        if max_off < tolerance || iterations >= max_iter { break; }

        let app = a[p][p];
        let aqq = a[q][q];
        let apq = a[p][q];
        if apq.abs() < 1e-20 { iterations += 1; continue; }

        let theta = (aqq - app) / (2.0 * apq);
        let t = if theta >= 0.0 {
            1.0 / (theta + (1.0 + theta * theta).sqrt())
        } else {
            1.0 / (theta - (1.0 + theta * theta).sqrt())
        };
        let c = 1.0 / (1.0 + t * t).sqrt();
        let s = t * c;

        a[p][p] = app - t * apq;
        a[q][q] = aqq + t * apq;
        a[p][q] = 0.0;
        a[q][p] = 0.0;

        for i in 0..n {
            if i != p && i != q {
                let aip = a[i][p];
                let aiq = a[i][q];
                let nip = c * aip - s * aiq;
                let niq = s * aip + c * aiq;
                a[i][p] = nip; a[p][i] = nip;
                a[i][q] = niq; a[q][i] = niq;
            }
        }
        for i in 0..n {
            let vip = v[i][p];
            let viq = v[i][q];
            v[i][p] = c * vip - s * viq;
            v[i][q] = s * vip + c * viq;
        }
        iterations += 1;
    }

    let mut eigenvalues: Vec<f64> = (0..n).map(|i| a[i][i]).collect();
    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by(|&i, &j| eigenvalues[i].partial_cmp(&eigenvalues[j]).unwrap());

    let sorted_eigenvalues: Vec<f64> = indices.iter().map(|&i| eigenvalues[i]).collect();
    let sorted_eigenvectors: Vec<Vec<f64>> = indices.iter().map(|&i| {
        (0..n).map(|j| v[j][i]).collect()
    }).collect();

    EigenResult {
        eigenvalues: sorted_eigenvalues,
        eigenvectors: sorted_eigenvectors,
        algorithm: Algorithm::SimpleHuckel,
    }
}

pub fn solve_generalized_eigen(
    h: &Vec<Vec<f64>>,
    s: &Vec<Vec<f64>>,
    max_iter: usize,
    tolerance: f64,
) -> EigenResult {
    let mut h = ensure_square(h);
    let mut s = ensure_square(s);
    sanitize_matrix(&mut h);
    sanitize_matrix(&mut s);
    symmetrize(&mut h);
    symmetrize(&mut s);

    let n = h.len();
    if n == 0 {
        return EigenResult {
            eigenvalues: vec![],
            eigenvectors: vec![],
            algorithm: Algorithm::ExtendedHuckel,
        };
    }

    regularize_positive_definite(&mut s, 1e-4);

    let s_inv_sqrt = match matrix_inverse_sqrt_safe(&s) {
        Some(m) => m,
        None => {
            eprintln!("[WARN] S^(-1/2) failed, using diagonal approximation");
            let mut diag = vec![vec![0.0f64; n]; n];
            for i in 0..n { diag[i][i] = 1.0 / s[i][i].sqrt().max(1e-6); }
            diag
        }
    };

    let mut h_prime = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let mut sum = 0.0;
            for k in 0..n {
                for l in 0..n {
                    sum += s_inv_sqrt[i][k] * h[k][l] * s_inv_sqrt[l][j];
                }
            }
            h_prime[i][j] = sum;
            h_prime[j][i] = sum;
        }
    }
    sanitize_matrix(&mut h_prime);
    symmetrize(&mut h_prime);

    let intermediate = solve_eigen_jacobi(&h_prime, max_iter, tolerance);

    let num_eig = intermediate.eigenvalues.len();
    let mut eigenvectors = vec![vec![0.0f64; n]; num_eig];
    for (eig_idx, eig_vec) in intermediate.eigenvectors.iter().enumerate() {
        for i in 0..n {
            let mut sum = 0.0;
            for k in 0..n {
                sum += s_inv_sqrt[i][k] * eig_vec[k];
            }
            eigenvectors[eig_idx][i] = sum;
        }
        let norm_sq: f64 = eigenvectors[eig_idx].iter()
            .enumerate()
            .map(|(i, vi)| vi * (0..n).map(|j| s[i][j] * eigenvectors[eig_idx][j]).sum::<f64>())
            .sum();
        if norm_sq > 1e-20 {
            let inv_norm = 1.0 / norm_sq.sqrt();
            for v in eigenvectors[eig_idx].iter_mut() { *v *= inv_norm; }
        }
    }

    EigenResult {
        eigenvalues: intermediate.eigenvalues,
        eigenvectors,
        algorithm: Algorithm::ExtendedHuckel,
    }
}

fn matrix_cholesky_safe(matrix: &Vec<Vec<f64>>) -> Option<Vec<Vec<f64>>> {
    let n = matrix.len();
    let mut l = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let mut sum = 0.0;
            for k in 0..j { sum += l[i][k] * l[j][k]; }
            if i == j {
                let diag = matrix[i][i] - sum;
                if diag <= 0.0 || !diag.is_finite() {
                    return None;
                }
                l[i][j] = diag.sqrt();
            } else {
                if l[j][j] < 1e-20 { return None; }
                l[i][j] = (matrix[i][j] - sum) / l[j][j];
            }
        }
    }
    Some(l)
}

fn matrix_inverse_sqrt_safe(matrix: &Vec<Vec<f64>>) -> Option<Vec<Vec<f64>>> {
    let n = matrix.len();
    let l = matrix_cholesky_safe(matrix)?;
    let l_inv = matrix_inverse_triangular_safe(&l)?;
    let mut result = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let mut sum = 0.0;
            for k in j..n { sum += l_inv[k][i] * l_inv[k][j]; }
            if sum.is_finite() {
                result[i][j] = sum;
                result[j][i] = sum;
            } else {
                return None;
            }
        }
    }
    Some(result)
}

fn matrix_inverse_triangular_safe(l: &Vec<Vec<f64>>) -> Option<Vec<Vec<f64>>> {
    let n = l.len();
    let mut inv = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        if l[i][i].abs() < 1e-20 { return None; }
        inv[i][i] = 1.0 / l[i][i];
        for j in (i + 1)..n {
            let mut sum = 0.0;
            for k in i..j { sum += l[j][k] * inv[k][i]; }
            inv[j][i] = -sum / l[j][j];
            if !inv[j][i].is_finite() { return None; }
        }
    }
    Some(inv)
}
