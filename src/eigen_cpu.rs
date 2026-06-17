use crate::types::EigenResult;

pub fn solve_eigen_jacobi(matrix: &Vec<Vec<f64>>, max_iter: usize, tolerance: f64) -> EigenResult {
    let n = matrix.len();
    let mut a = matrix.clone();
    let mut v = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        v[i][i] = 1.0;
    }

    let mut iterations = 0;
    loop {
        let mut max_off = 0.0;
        let mut p = 0;
        let mut q = 1;
        for i in 0..n {
            for j in (i + 1)..n {
                let abs_val = a[i][j].abs();
                if abs_val > max_off {
                    max_off = abs_val;
                    p = i;
                    q = j;
                }
            }
        }

        if max_off < tolerance || iterations >= max_iter {
            break;
        }

        let app = a[p][p];
        let aqq = a[q][q];
        let apq = a[p][q];
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
                a[i][p] = c * aip - s * aiq;
                a[p][i] = a[i][p];
                a[i][q] = s * aip + c * aiq;
                a[q][i] = a[i][q];
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
    let sorted_eigenvectors: Vec<Vec<f64>> = indices
        .iter()
        .map(|&i| {
            let mut col = vec![0.0f64; n];
            for j in 0..n {
                col[j] = v[j][i];
            }
            col
        })
        .collect();

    EigenResult {
        eigenvalues: sorted_eigenvalues,
        eigenvectors: sorted_eigenvectors,
    }
}

pub fn solve_generalized_eigen(
    h: &Vec<Vec<f64>>,
    s: &Vec<Vec<f64>>,
    max_iter: usize,
    tolerance: f64,
) -> EigenResult {
    let n = h.len();
    let s_inv_sqrt = matrix_inverse_sqrt(s);

    let mut h_prime = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in 0..n {
            let mut sum = 0.0;
            for k in 0..n {
                for l in 0..n {
                    sum += s_inv_sqrt[i][k] * h[k][l] * s_inv_sqrt[l][j];
                }
            }
            h_prime[i][j] = sum;
        }
    }

    let intermediate = solve_eigen_jacobi(&h_prime, max_iter, tolerance);

    let mut eigenvectors = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in 0..n {
            let mut sum = 0.0;
            for k in 0..n {
                sum += s_inv_sqrt[i][k] * intermediate.eigenvectors[j][k];
            }
            eigenvectors[j][i] = sum;
        }
    }

    EigenResult {
        eigenvalues: intermediate.eigenvalues,
        eigenvectors,
    }
}

fn matrix_cholesky(matrix: &Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    let n = matrix.len();
    let mut l = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let mut sum = 0.0;
            for k in 0..j {
                sum += l[i][k] * l[j][k];
            }
            if i == j {
                l[i][j] = (matrix[i][i] - sum).sqrt();
            } else {
                l[i][j] = (matrix[i][j] - sum) / l[j][j];
            }
        }
    }
    l
}

fn matrix_inverse_sqrt(matrix: &Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    let n = matrix.len();
    let l = matrix_cholesky(matrix);
    let l_inv = matrix_inverse_triangular(&l);
    let mut result = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in 0..n {
            let mut sum = 0.0;
            for k in j..n {
                sum += l_inv[k][i] * l_inv[k][j];
            }
            result[i][j] = sum;
            result[j][i] = sum;
        }
    }
    result
}

fn matrix_inverse_triangular(l: &Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    let n = l.len();
    let mut inv = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        inv[i][i] = 1.0 / l[i][i];
        for j in (i + 1)..n {
            let mut sum = 0.0;
            for k in i..j {
                sum += l[j][k] * inv[k][i];
            }
            inv[j][i] = -sum / l[j][j];
        }
    }
    inv
}
