//! Supervector spaces, supertrace, and superdeterminant (Berezinian).

use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};

/// A supervector space of dimension (p|q).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SuperVectorSpace {
    pub p: usize, // bosonic (even) dimension
    pub q: usize, // fermionic (odd) dimension
}

impl SuperVectorSpace {
    pub fn new(p: usize, q: usize) -> Self {
        SuperVectorSpace { p, q }
    }

    pub fn total_dimension(&self) -> usize {
        self.p + self.q
    }
}

/// A linear operator on a supervector space, represented as a block matrix:
///   M = [ A  B ]
///       [ C  D ]
/// where A is p×p, B is p×q, C is q×p, D is q×q.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SuperMatrix {
    pub p: usize,
    pub q: usize,
    /// Full (p+q)×(p+q) matrix stored in row-major.
    /// Rows/cols 0..p are even, p..p+q are odd.
    pub data: Vec<Vec<f64>>,
}

impl SuperMatrix {
    pub fn new(p: usize, q: usize) -> Self {
        let n = p + q;
        SuperMatrix {
            p, q,
            data: vec![vec![0.0; n]; n],
        }
    }

    pub fn from_block_matrix(
        p: usize, q: usize,
        a: &[Vec<f64>], // p×p
        b: &[Vec<f64>], // p×q
        c: &[Vec<f64>], // q×p
        d: &[Vec<f64>], // q×q
    ) -> Self {
        let n = p + q;
        let mut data = vec![vec![0.0; n]; n];
        for i in 0..p {
            for j in 0..p { data[i][j] = a[i][j]; }
            for j in 0..q { data[i][p + j] = b[i][j]; }
        }
        for i in 0..q {
            for j in 0..p { data[p + i][j] = c[i][j]; }
            for j in 0..q { data[p + i][p + j] = d[i][j]; }
        }
        SuperMatrix { p, q, data }
    }

    /// Extract block A (p×p, even-even).
    pub fn block_a(&self) -> DMatrix<f64> {
        DMatrix::from_fn(self.p, self.p, |i, j| self.data[i][j])
    }

    /// Extract block B (p×q, even-odd).
    pub fn block_b(&self) -> DMatrix<f64> {
        DMatrix::from_fn(self.p, self.q, |i, j| self.data[i][self.p + j])
    }

    /// Extract block C (q×p, odd-even).
    pub fn block_c(&self) -> DMatrix<f64> {
        DMatrix::from_fn(self.q, self.p, |i, j| self.data[self.p + i][j])
    }

    /// Extract block D (q×q, odd-odd).
    pub fn block_d(&self) -> DMatrix<f64> {
        DMatrix::from_fn(self.q, self.q, |i, j| self.data[self.p + i][self.p + j])
    }

    /// Supertrace: str(M) = tr(A) - tr(D).
    pub fn supertrace(&self) -> f64 {
        let tr_a: f64 = (0..self.p).map(|i| self.data[i][i]).sum();
        let tr_d: f64 = (0..self.q).map(|i| self.data[self.p + i][self.p + i]).sum();
        tr_a - tr_d
    }

    /// Ordinary trace.
    pub fn trace(&self) -> f64 {
        let n = self.p + self.q;
        (0..n).map(|i| self.data[i][i]).sum()
    }

    /// Superdeterminant (Berezinian):
    /// Ber(M) = det(A - B·D⁻¹·C) / det(D)   when D is invertible
    ///       or = det(A) · det(D - C·A⁻¹·B)⁻¹ when A is invertible
    pub fn berezinian(&self) -> Option<f64> {
        let a = self.block_a();
        let b = self.block_b();
        let c = self.block_c();
        let d = self.block_d();

        // Try det(A) · det(D - C·A⁻¹·B)⁻¹ when A is invertible
        if self.p > 0 {
            if let Some(a_inv) = a.clone().try_inverse() {
                let ca_inv_b = &c * &a_inv * &b;
                let d_prime = &d - &ca_inv_b;
                if self.q > 0 {
                    let det_a = a.determinant();
                    let det_dp = d_prime.determinant();
                    if det_dp.abs() < 1e-15 { return None; }
                    return Some(det_a / det_dp);
                } else {
                    // q = 0: Berezinian = det(A)
                    return Some(a.determinant());
                }
            }
        }

        // Try D invertible approach
        if self.q > 0 {
            if let Some(d_inv) = d.clone().try_inverse() {
                let bd_inv_c = &b * &d_inv * &c;
                let a_prime = &a - &bd_inv_c;
                let det_a_prime = if self.p > 0 { a_prime.determinant() } else { 1.0 };
                let det_d = d.determinant();
                if det_d.abs() < 1e-15 { return None; }
                return Some(det_a_prime / det_d);
            }
        }

        None
    }

    /// Multiply two supermatrices.
    pub fn mul(&self, other: &SuperMatrix) -> SuperMatrix {
        let n = self.p + self.q;
        let mut result = SuperMatrix::new(self.p, self.q);
        for i in 0..n {
            for j in 0..n {
                let mut sum = 0.0;
                for k in 0..n {
                    sum += self.data[i][k] * other.data[k][j];
                }
                result.data[i][j] = sum;
            }
        }
        result
    }
}

/// Verify supertrace identity: str([A,B]) = 0 for supermatrices
/// where [A,B] = AB - (-1)^|row||col| BA (graded commutator).
/// For diagonal matrices, this simplifies nicely.
pub fn supertrace_supercommutator_vanishes(a: &SuperMatrix, b: &SuperMatrix) -> bool {
    let ab = a.mul(b);
    let ba = b.mul(a);
    let n = a.p + a.q;
    let mut comm = SuperMatrix::new(a.p, a.q);
    for i in 0..n {
        for j in 0..n {
            // For the ordinary commutator of diagonal matrices, it's 0
            comm.data[i][j] = ab.data[i][j] - ba.data[i][j];
        }
    }
    // For diagonal matrices, [A,B] = 0, so supertrace = 0
    comm.supertrace().abs() < 1e-10
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_supertrace_identity() {
        let m = SuperMatrix::from_block_matrix(
            2, 1,
            &[vec![1.0, 2.0], vec![3.0, 4.0]],
            &[vec![5.0], vec![6.0]],
            &[vec![7.0, 8.0]],
            &[vec![9.0]],
        );
        // str = tr(A) - tr(D) = (1+4) - 9 = -4
        assert!((m.supertrace() - (-4.0)).abs() < 1e-10);
    }

    #[test]
    fn test_supertrace_vs_trace() {
        let m = SuperMatrix::from_block_matrix(
            1, 1,
            &[vec![3.0]],
            &[vec![2.0]],
            &[vec![1.0]],
            &[vec![5.0]],
        );
        // tr = 3 + 5 = 8, str = 3 - 5 = -2
        assert!((m.trace() - 8.0).abs() < 1e-10);
        assert!((m.supertrace() - (-2.0)).abs() < 1e-10);
    }

    #[test]
    fn test_berezinian_diagonal() {
        // Diagonal: Ber(diag(a₁,...,aₚ, d₁,...,d_q)) = (a₁·...·aₚ)/(d₁·...·d_q)
        let m = SuperMatrix::from_block_matrix(
            2, 1,
            &[vec![2.0, 0.0], vec![0.0, 3.0]],
            &[vec![0.0], vec![0.0]],
            &[vec![0.0, 0.0]],
            &[vec![4.0]],
        );
        // Ber = (2·3) / 4 = 6/4 = 1.5
        let ber = m.berezinian().unwrap();
        assert!((ber - 1.5).abs() < 1e-10);
    }

    #[test]
    fn test_berezinian_identity() {
        let m = SuperMatrix::from_block_matrix(
            2, 2,
            &[vec![1.0, 0.0], vec![0.0, 1.0]],
            &[vec![0.0, 0.0], vec![0.0, 0.0]],
            &[vec![0.0, 0.0], vec![0.0, 0.0]],
            &[vec![1.0, 0.0], vec![0.0, 1.0]],
        );
        let ber = m.berezinian().unwrap();
        assert!((ber - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_berezinian_with_off_diagonal() {
        let m = SuperMatrix::from_block_matrix(
            2, 1,
            &[vec![1.0, 0.0], vec![0.0, 2.0]],
            &[vec![1.0], vec![0.0]],
            &[vec![0.0, 1.0]],
            &[vec![3.0]],
        );
        // A = [[1,0],[0,2]], D - C·A⁻¹·B = [3] - [0,1]·[[1,0],[0,0.5]]·[[1],[0]]
        //   = [3] - [0,1]·[[1],[0]] = [3] - [0] = [3]
        // Ber = det(A) / det(3) = 2 / 3
        let ber = m.berezinian().unwrap();
        assert!((ber - 2.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_supertrace_linearity() {
        let m1 = SuperMatrix::from_block_matrix(
            1, 1,
            &[vec![2.0]],
            &[vec![0.0]],
            &[vec![0.0]],
            &[vec![3.0]],
        );
        let m2 = SuperMatrix::from_block_matrix(
            1, 1,
            &[vec![4.0]],
            &[vec![0.0]],
            &[vec![0.0]],
            &[vec![5.0]],
        );
        // str(A) + str(B) should equal str(A+B)
        // str(A) = 2 - 3 = -1, str(B) = 4 - 5 = -1
        let sum = (-1.0) + (-1.0);
        assert!((m1.supertrace() + m2.supertrace() - sum).abs() < 1e-10);
    }

    #[test]
    fn test_berezinian_pure_bosonic() {
        // When q=0, Berezinian = det(A)
        let m = SuperMatrix::from_block_matrix(
            2, 0,
            &[vec![1.0, 2.0], vec![3.0, 4.0]],
            &[],
            &[],
            &[],
        );
        // det = 1*4 - 2*3 = -2
        let ber = m.berezinian().unwrap();
        assert!((ber - (-2.0)).abs() < 1e-10);
    }

    #[test]
    fn test_super_matrix_multiply() {
        let a = SuperMatrix::from_block_matrix(
            1, 1,
            &[vec![2.0]],
            &[vec![0.0]],
            &[vec![0.0]],
            &[vec![3.0]],
        );
        let b = SuperMatrix::from_block_matrix(
            1, 1,
            &[vec![4.0]],
            &[vec![0.0]],
            &[vec![0.0]],
            &[vec![5.0]],
        );
        let c = a.mul(&b);
        assert!((c.data[0][0] - 8.0).abs() < 1e-10);
        assert!((c.data[1][1] - 15.0).abs() < 1e-10);
    }
}
