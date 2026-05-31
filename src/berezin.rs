//! Berezin integral: integration on Grassmann algebra.
//!
//! The Berezin integral of a Grassmann variable θ is defined as:
//!   ∫ θ dθ = 1  (and ∫ 1 dθ = 0)
//!
//! For n variables: ∫ θ₁θ₂...θₙ dθₙ...dθ₁ = 1

use crate::exterior::ExteriorAlgebra;
use serde::{Deserialize, Serialize};

/// Berezin integral over a Grassmann algebra with n generators.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BerezinIntegrator {
    pub n_variables: usize,
    pub labels: Vec<String>,
}

impl BerezinIntegrator {
    pub fn new(n: usize, labels: Vec<String>) -> Self {
        let labels = if labels.len() == n { labels } else {
            (0..n).map(|i| format!("θ{}", i)).collect()
        };
        BerezinIntegrator { n_variables: n, labels }
    }

    /// Integrate over a single Grassmann variable.
    /// A function of θ is f = a + bθ.
    /// ∫(a + bθ) dθ = b (coefficient of θ).
    pub fn integrate_single(coeff_constant: f64, coeff_theta: f64) -> f64 {
        coeff_theta
    }

    /// Integrate over all n Grassmann variables.
    /// The function is represented as coefficients indexed by θ-monomial bitmask.
    /// ∫ f(θ₁,...,θₙ) dθₙ...dθ₁ = coefficient of θ₁θ₂...θₙ.
    pub fn integrate(&self, coefficients: &std::collections::HashMap<u64, f64>) -> f64 {
        let top_mask = (1u64 << self.n_variables) - 1;
        coefficients.get(&top_mask).copied().unwrap_or(0.0)
    }

    /// Iterated Berezin integral: integrate over variable `var_index`.
    /// For a function f(θ₀,...,θ_{n-1}), integrating over θ_k:
    /// Extract the coefficient of θ_k and drop θ_k from all monomials.
    pub fn integrate_over(&self, coefficients: &std::collections::HashMap<u64, f64>, var_index: usize) -> std::collections::HashMap<u64, f64> {
        let bit = 1u64 << var_index;
        let mut result = std::collections::HashMap::new();
        for (&mask, &val) in coefficients {
            if mask & bit != 0 {
                // This monomial contains θ_{var_index}, so it contributes
                let reduced_mask = mask ^ bit; // remove the bit
                // Sign from moving θ_{var_index} to the front
                let sign = if Self::count_bits_below(mask, var_index) % 2 == 0 { 1.0 } else { -1.0 };
                *result.entry(reduced_mask).or_insert(0.0) += sign * val;
            }
            // Monomials without θ_{var_index} contribute 0 (their derivative)
        }
        result
    }

    fn count_bits_below(mask: u64, bit: usize) -> u32 {
        let below_mask = (1u64 << bit) - 1;
        (mask & below_mask).count_ones()
    }

    /// Full Berezin integral by iterated integration.
    pub fn integrate_iterated(&self, coefficients: &std::collections::HashMap<u64, f64>) -> f64 {
        let mut coeffs = coefficients.clone();
        // Integrate θ₀, then θ₁, ..., then θ_{n-1}
        for i in 0..self.n_variables {
            coeffs = self.integrate_over(&coeffs, i);
        }
        // Final result should be a single number (mask = 0)
        coeffs.get(&0).copied().unwrap_or(0.0)
    }

    /// Change of variables formula for Berezin integral:
    /// ∫ f(θ) dθ = ∫ f(θ') Ber(J) dθ'
    /// where Ber(J) is the Berezinian of the Jacobian.
    pub fn change_of_variables_sign(&self, is_orientation_preserving: bool) -> f64 {
        if is_orientation_preserving { 1.0 } else { -1.0 }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    #[test]
    fn test_berezin_single_theta() {
        // ∫ θ dθ = 1
        assert!((BerezinIntegrator::integrate_single(0.0, 1.0) - 1.0).abs() < 1e-10);
        // ∫ 1 dθ = 0
        assert!(BerezinIntegrator::integrate_single(1.0, 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_berezin_two_variables() {
        // ∫ θ₁θ₂ dθ₂ dθ₁ = 1
        let b = BerezinIntegrator::new(2, vec![]);
        let mut coeffs = HashMap::new();
        coeffs.insert(0b11, 1.0); // θ₁θ₂
        assert!((b.integrate(&coeffs) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_berezin_iterated() {
        let b = BerezinIntegrator::new(3, vec![]);
        let mut coeffs = HashMap::new();
        coeffs.insert(0b111, 2.0); // 2·θ₁θ₂θ₃
        assert!((b.integrate_iterated(&coeffs) - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_berezin_no_top_form() {
        let b = BerezinIntegrator::new(2, vec![]);
        let mut coeffs = HashMap::new();
        coeffs.insert(0b01, 1.0); // just θ₁, not θ₁θ₂
        assert!(b.integrate(&coeffs).abs() < 1e-10);
    }

    #[test]
    fn test_integrate_over_single() {
        let b = BerezinIntegrator::new(2, vec![]);
        let mut coeffs = HashMap::new();
        coeffs.insert(0b11, 1.0); // θ₁θ₂
        let result = b.integrate_over(&coeffs, 0); // integrate θ₁
        assert!((result.get(&0b10).unwrap() - 1.0).abs() < 1e-10); // left with θ₂
    }

    #[test]
    fn test_berezin_linear_combination() {
        // ∫ (3θ₁ + 2θ₂) dθ₂dθ₁ = 0 (no θ₁θ₂ term)
        let b = BerezinIntegrator::new(2, vec![]);
        let mut coeffs = HashMap::new();
        coeffs.insert(0b01, 3.0);
        coeffs.insert(0b10, 2.0);
        assert!(b.integrate(&coeffs).abs() < 1e-10);
    }

    #[test]
    fn test_berezin_1var() {
        let b = BerezinIntegrator::new(1, vec![]);
        let mut coeffs = HashMap::new();
        coeffs.insert(0b1, 7.0);
        assert!((b.integrate(&coeffs) - 7.0).abs() < 1e-10);
    }

    #[test]
    fn test_berezin_4var() {
        let b = BerezinIntegrator::new(4, vec![]);
        let mut coeffs = HashMap::new();
        coeffs.insert(0b1111, 1.0); // θ₁θ₂θ₃θ₄
        assert!((b.integrate(&coeffs) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_berezin_general_function() {
        // f(θ₁,θ₂) = 3 + 2θ₁ + θ₂ + 5θ₁θ₂
        // ∫ f dθ₂dθ₁ = 5
        let b = BerezinIntegrator::new(2, vec![]);
        let mut coeffs = HashMap::new();
        coeffs.insert(0b00, 3.0);
        coeffs.insert(0b01, 2.0);
        coeffs.insert(0b10, 1.0);
        coeffs.insert(0b11, 5.0);
        assert!((b.integrate(&coeffs) - 5.0).abs() < 1e-10);
    }
}
