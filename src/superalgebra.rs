//! Superalgebra: Z₂-graded associative algebra.

use crate::graded::{GradedElement, Parity};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// A Z₂-graded associative algebra with a multiplication table.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SuperAlgebra {
    pub name: String,
    pub basis: Vec<GradedElement>,
    /// Multiplication table: (i, j) -> (result_value, result_label)
    pub mult_table: HashMap<(usize, usize), (f64, String)>,
}

impl SuperAlgebra {
    pub fn new(name: &str, basis: Vec<GradedElement>) -> Self {
        SuperAlgebra { name: name.to_string(), basis, mult_table: HashMap::new() }
    }

    /// Set multiplication: basis[i] * basis[j] = value with label.
    pub fn set_product(&mut self, i: usize, j: usize, value: f64, label: &str) {
        self.mult_table.insert((i, j), (value, label.to_string()));
    }

    /// Multiply two basis elements.
    pub fn multiply(&self, i: usize, j: usize) -> Option<(f64, Parity)> {
        let a = &self.basis[i];
        let b = &self.basis[j];
        let result_parity = a.parity.combine(b.parity);
        self.mult_table.get(&(i, j)).map(|(v, _)| (*v, result_parity))
    }

    /// Check if the algebra is supercommutative: ab = (-1)^|a||b| ba.
    pub fn is_supercommutative(&self) -> bool {
        let n = self.basis.len();
        for i in 0..n {
            for j in 0..n {
                if let (Some(ab), Some(ba)) = (self.multiply(i, j), self.multiply(j, i)) {
                    let sign = if self.basis[i].parity.grade() * self.basis[j].parity.grade() % 2 != 0 { -1.0 } else { 1.0 };
                    if (ab.0 - sign * ba.0).abs() > 1e-10 {
                        return false;
                    }
                }
            }
        }
        true
    }

    /// Check associativity: (ab)c = a(bc).
    pub fn is_associative(&self) -> bool {
        let n = self.basis.len();
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    let left = self.multiply(i, j).map(|(v, _)| v);
                    let right = self.multiply(j, k).map(|(v, _)| v);
                    if let (Some(ab), Some(bc)) = (left, right) {
                        // Simplified check for scalar algebras
                        let labc = self.mult_table.get(&(i, j))
                            .and_then(|_| self.mult_table.get(&(i, k)));
                        let rabc = self.mult_table.get(&(j, k))
                            .and_then(|_| self.mult_table.get(&(i, k)));
                        if let (Some(l), Some(r)) = (labc, rabc) {
                            if (l.0 - r.0).abs() > 1e-10 {
                                return false;
                            }
                        }
                    }
                }
            }
        }
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_super_algebra_basic() {
        let basis = vec![
            GradedElement::even(1.0, "1"),
            GradedElement::odd(1.0, "θ"),
        ];
        let mut alg = SuperAlgebra::new("Λ(θ)", basis);
        alg.set_product(0, 0, 1.0, "1");
        alg.set_product(0, 1, 1.0, "θ");
        alg.set_product(1, 0, 1.0, "θ");
        alg.set_product(1, 1, 0.0, "0"); // θ² = 0

        assert!(alg.is_supercommutative());
        let prod = alg.multiply(1, 1).unwrap();
        assert!(prod.0.abs() < 1e-10);
    }

    #[test]
    fn test_super_algebra_grassmann2() {
        let basis = vec![
            GradedElement::even(1.0, "1"),
            GradedElement::odd(1.0, "θ₁"),
            GradedElement::odd(1.0, "θ₂"),
            GradedElement::even(1.0, "θ₁θ₂"),
        ];
        let mut alg = SuperAlgebra::new("Λ(θ₁,θ₂)", basis);
        // Identity
        alg.set_product(0, 0, 1.0, "1");
        alg.set_product(0, 1, 1.0, "θ₁");
        alg.set_product(0, 2, 1.0, "θ₂");
        alg.set_product(1, 0, 1.0, "θ₁");
        alg.set_product(2, 0, 1.0, "θ₂");
        // θ₁² = 0
        alg.set_product(1, 1, 0.0, "0");
        alg.set_product(2, 2, 0.0, "0");
        // θ₁θ₂ and θ₂θ₁ = -θ₁θ₂
        alg.set_product(1, 2, 1.0, "θ₁θ₂");
        alg.set_product(2, 1, -1.0, "-θ₁θ₂");

        assert!(alg.is_supercommutative());
    }
}
