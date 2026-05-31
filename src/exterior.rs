//! Exterior algebra Λ(V) as a superalgebra.

use crate::graded::{GradedElement, Parity};
use serde::{Deserialize, Serialize};

/// Exterior algebra on n generators.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExteriorAlgebra {
    pub n_generators: usize,
    pub generator_labels: Vec<String>,
}

impl ExteriorAlgebra {
    pub fn new(n: usize, labels: Vec<String>) -> Self {
        let labels = if labels.len() == n { labels } else {
            (0..n).map(|i| format!("e{}", i)).collect()
        };
        ExteriorAlgebra { n_generators: n, generator_labels: labels }
    }

    /// Dimension of the exterior algebra: 2^n.
    pub fn dimension(&self) -> usize {
        1 << self.n_generators
    }

    /// Dimension of the k-forms subspace: C(n, k).
    pub fn form_dimension(&self, k: usize) -> usize {
        if k > self.n_generators { 0 } else {
            let mut binom: usize = 1;
            for i in 0..k.min(self.n_generators - k) {
                binom = binom * (self.n_generators - i) / (i + 1);
            }
            binom
        }
    }

    /// Wedge product of two monomials represented as bitmask.
    /// Returns None if they share a generator (result is 0).
    pub fn wedge(&self, a_mask: u64, b_mask: u64) -> Option<u64> {
        if a_mask & b_mask != 0 {
            None // anticommute => 0 if same generator appears
        } else {
            Some(a_mask | b_mask)
        }
    }

    /// Sign from reordering: counts transpositions to merge two sorted sequences.
    /// If a has bits set at positions a₁ < a₂ < ... and b at b₁ < b₂ < ...,
    /// the Koszul sign is (-1)^(number of pairs (i,j) with aᵢ > bⱼ).
    pub fn koszul_sign(&self, a_mask: u64, b_mask: u64) -> i64 {
        let mut sign: i64 = 1;
        let mut a_bits: Vec<u64> = vec![];
        let mut b_bits: Vec<u64> = vec![];
        for i in 0..self.n_generators {
            if a_mask & (1 << i) != 0 { a_bits.push(i as u64); }
            if b_mask & (1 << i) != 0 { b_bits.push(i as u64); }
        }
        for &ai in &a_bits {
            for &bj in &b_bits {
                if ai > bj { sign *= -1; }
            }
        }
        sign
    }

    /// Full wedge product with Koszul sign.
    pub fn signed_wedge(&self, a_mask: u64, b_mask: u64) -> (i64, u64) {
        if a_mask & b_mask != 0 {
            (0, 0) // zero
        } else {
            (self.koszul_sign(a_mask, b_mask), a_mask | b_mask)
        }
    }

    /// The exterior algebra is a superalgebra with even/odd grading by form degree.
    /// All generators are odd.
    pub fn generators_as_graded(&self) -> Vec<GradedElement> {
        self.generator_labels.iter().enumerate().map(|(i, l)| {
            GradedElement::odd(1.0, &format!("{}(d{})", l, i))
        }).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_exterior_dimension() {
        let e = ExteriorAlgebra::new(3, vec!["e0".into(), "e1".into(), "e2".into()]);
        assert_eq!(e.dimension(), 8);
    }

    #[test]
    fn test_form_dimension() {
        let e = ExteriorAlgebra::new(4, vec![]);
        assert_eq!(e.form_dimension(0), 1);
        assert_eq!(e.form_dimension(1), 4);
        assert_eq!(e.form_dimension(2), 6);
        assert_eq!(e.form_dimension(3), 4);
        assert_eq!(e.form_dimension(4), 1);
    }

    #[test]
    fn test_wedge_product() {
        let e = ExteriorAlgebra::new(3, vec!["e0".into(), "e1".into(), "e2".into()]);
        assert_eq!(e.wedge(0b001, 0b010), Some(0b011));
        assert_eq!(e.wedge(0b001, 0b001), None); // e0 ∧ e0 = 0
    }

    #[test]
    fn test_exterior_is_superalgebra() {
        let e = ExteriorAlgebra::new(2, vec!["θ".into(), "η".into()]);
        let gens = e.generators_as_graded();
        assert!(gens[0].is_odd());
        assert!(gens[1].is_odd());
    }

    #[test]
    fn test_koszul_sign() {
        let e = ExteriorAlgebra::new(3, vec![]);
        // e0 ∧ e1: no transpositions needed => +1
        assert_eq!(e.koszul_sign(0b001, 0b010), 1);
        // e1 ∧ e0: transposition needed => -1
        assert_eq!(e.koszul_sign(0b010, 0b001), -1);
    }
}
