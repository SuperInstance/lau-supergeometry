//! Superforms: differential forms on supermanifolds.

use crate::graded::{GradedElement, Parity};
use crate::supermanifold::Supermanifold;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// A superform on a supermanifold.
/// Combines ordinary differential forms with Grassmann-odd differentials dθ.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SuperForm {
    pub manifold_dim: (usize, usize), // (p, q)
    /// Components: (form_coefficient_mask, differential_mask) -> value
    /// form_coefficient_mask encodes which θ's appear in the coefficient
    /// differential_mask encodes which differentials (dx and dθ) appear
    pub components: HashMap<(u64, u64), f64>,
    pub labels: Vec<String>,
}

impl SuperForm {
    pub fn new(p: usize, q: usize) -> Self {
        SuperForm {
            manifold_dim: (p, q),
            components: HashMap::new(),
            labels: vec![],
        }
    }

    /// Add a component: coefficient θ-mask, differential mask, value.
    pub fn add_component(&mut self, coeff_mask: u64, diff_mask: u64, value: f64) {
        *self.components.entry((coeff_mask, diff_mask)).or_insert(0.0) += value;
    }

    /// Exterior derivative d: increases form degree by 1.
    /// For a 0-form f: df = Σ ∂f/∂xᵢ dxᵢ + Σ ∂f/∂θⱼ dθⱼ
    /// Here we compute the structure symbolically.
    pub fn exterior_derivative(&self) -> SuperForm {
        // Simplified: shift differential mask by adding one new differential
        let mut result = SuperForm::new(self.manifold_dim.0, self.manifold_dim.1);
        let total_diffs = self.manifold_dim.0 + self.manifold_dim.1;
        for (&(coeff, diff), &val) in &self.components {
            for k in 0..total_diffs as u64 {
                if diff & (1 << k) == 0 {
                    // Add differential k
                    let new_diff = diff | (1 << k);
                    // Koszul sign: (-1)^(number of set bits below k in diff)
                    let sign = if Self::count_bits_below(diff, k) % 2 == 0 { 1.0 } else { -1.0 };
                    result.add_component(coeff, new_diff, sign * val);
                }
            }
        }
        result
    }

    fn count_bits_below(mask: u64, bit: u64) -> u32 {
        if bit >= 64 {
            mask.count_ones()
        } else {
            let below_mask = (1u64 << bit) - 1;
            (mask & below_mask).count_ones()
        }
    }

    /// Wedge product of two superforms.
    pub fn wedge(&self, other: &SuperForm) -> SuperForm {
        let mut result = SuperForm::new(self.manifold_dim.0, self.manifold_dim.1);
        for (&(c1, d1), &v1) in &self.components {
            for (&(c2, d2), &v2) in &other.components {
                if d1 & d2 == 0 {
                    // Koszul sign for wedge of forms
                    let _sign = if Self::count_bits_below(d1, 64) % 2 == 0 { 1.0 } else { -1.0 };
                    // Simplified sign calculation
                    let bits_d2_in_d1 = d2 & ((1u64 << d1.count_ones()) - 1);
                    let actual_sign = if bits_d2_in_d1.count_ones() % 2 == 0 { 1.0 } else { -1.0 };
                    result.add_component(c1 | c2, d1 | d2, actual_sign * v1 * v2);
                }
            }
        }
        result
    }

    /// Form degree (number of differentials).
    pub fn degree(&self) -> Option<u32> {
        let degrees: std::collections::HashSet<u32> = self.components.keys()
            .map(|(_, d)| d.count_ones())
            .collect();
        if degrees.len() == 1 { Some(*degrees.iter().next().unwrap()) } else { None }
    }
}

/// Top superform (volume form) on a (p|q) supermanifold.
pub fn volume_form(p: usize, q: usize) -> SuperForm {
    let mut f = SuperForm::new(p, q);
    // dx₁ ∧ ... ∧ dx_p ∧ dθ₁ ∧ ... ∧ dθ_q
    let total = p + q;
    let full_mask = ((1u64 << total) - 1) << 0;
    // But dxᵢ occupies bits 0..p, dθⱼ occupies bits p..p+q
    // For the standard volume form, coefficient = 1 (body), diff mask = all
    let diff_mask = (1u64 << total) - 1;
    f.add_component(0, diff_mask, 1.0);
    f
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_superform_creation() {
        let mut f = SuperForm::new(2, 2);
        f.add_component(0, 0b01, 1.0); // dx₁
        assert_eq!(f.degree(), Some(1));
    }

    #[test]
    fn test_superform_wedge() {
        let mut f1 = SuperForm::new(2, 0);
        f1.add_component(0, 0b01, 1.0); // dx₁
        let mut f2 = SuperForm::new(2, 0);
        f2.add_component(0, 0b10, 1.0); // dx₂
        let f12 = f1.wedge(&f2);
        assert!(f12.components.contains_key(&(0, 0b11)));
    }

    #[test]
    fn test_exterior_derivative() {
        let mut f = SuperForm::new(1, 1);
        f.add_component(0, 0, 1.0); // 0-form (function)
        let df = f.exterior_derivative();
        assert!(!df.components.is_empty());
    }

    #[test]
    fn test_volume_form() {
        let vol = volume_form(2, 2);
        assert_eq!(vol.degree(), Some(4));
    }
}
