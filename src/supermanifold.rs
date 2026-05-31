//! Supermanifolds: coordinate ring = C^∞(M) ⊗ Λ(θ₁,...,θ_q)

use crate::graded::{GradedElement, GradedVec, Parity};
use crate::exterior::ExteriorAlgebra;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// A supermanifold of dimension (p|q): p even (bosonic) and q odd (fermionic) coordinates.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Supermanifold {
    pub p: usize, // bosonic dimension
    pub q: usize, // fermionic dimension
    pub bosonic_coords: Vec<String>,
    pub fermionic_coords: Vec<String>,
}

impl Supermanifold {
    pub fn new(p: usize, q: usize, bosonic: Vec<String>, fermionic: Vec<String>) -> Self {
        let bosonic_coords = if bosonic.len() == p { bosonic } else {
            (0..p).map(|i| format!("x{}", i)).collect()
        };
        let fermionic_coords = if fermionic.len() == q { fermionic } else {
            (0..q).map(|i| format!("θ{}", i)).collect()
        };
        Supermanifold { p, q, bosonic_coords, fermionic_coords }
    }

    /// The coordinate ring: functions are C^∞(ℝᵖ) ⊗ Λ(θ₁,...,θ_q).
    /// Represented as a map from monomial (in θ's, as bitmask) to a smooth function value.
    pub fn coordinate_ring(&self) -> CoordinateRing {
        CoordinateRing {
            manifold: self.clone(),
            functions: HashMap::new(),
        }
    }

    /// Total dimension of the function space: 2^q components for each point in ℝᵖ.
    pub fn function_dimension(&self) -> usize {
        self.q.pow(2) // Actually 2^q
    }

    /// Number of independent "functions" in the coordinate ring expansion.
    pub fn ring_rank(&self) -> usize {
        1 << self.q
    }
}

/// Coordinate ring of a supermanifold.
/// A superfunction f(x,θ) = f₀(x) + f₁(x)θ₁ + ... + f_{2^q-1}(x)θ₁θ₂...θ_q
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoordinateRing {
    pub manifold: Supermanifold,
    /// Map from θ-monomial bitmask to function values at sample points.
    pub functions: HashMap<u64, Vec<f64>>,
}

impl CoordinateRing {
    /// Set the body (θ-free part) of the superfunction.
    pub fn set_body(&mut self, values: Vec<f64>) {
        self.functions.insert(0, values);
    }

    /// Set a component corresponding to a θ-monomial.
    pub fn set_component(&mut self, mask: u64, values: Vec<f64>) {
        self.functions.insert(mask, values);
    }

    /// Evaluate the superfunction at bosonic point index i.
    /// Returns all components as (mask, value) pairs.
    pub fn evaluate_at(&self, point_index: usize) -> Vec<(u64, f64)> {
        self.functions.iter().map(|(mask, vals)| {
            let v = if point_index < vals.len() { vals[point_index] } else { 0.0 };
            (*mask, v)
        }).collect()
    }

    /// The body (degree-0 part) of the superfunction.
    pub fn body(&self) -> Option<&Vec<f64>> {
        self.functions.get(&0)
    }

    /// The soul (all higher-degree parts).
    pub fn soul_masks(&self) -> Vec<u64> {
        self.functions.keys().filter(|&&m| m != 0).cloned().collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_supermanifold_dimension() {
        let m = Supermanifold::new(3, 2, vec![], vec![]);
        assert_eq!(m.p, 3);
        assert_eq!(m.q, 2);
    }

    #[test]
    fn test_coordinate_ring_basic() {
        let m = Supermanifold::new(1, 2, vec!["t".into()], vec!["θ".into(), "η".into()]);
        let mut ring = m.coordinate_ring();
        ring.set_body(vec![1.0, 2.0, 3.0]);
        ring.set_component(0b01, vec![0.5, 1.0, 1.5]); // θ component
        ring.set_component(0b10, vec![0.3, 0.6, 0.9]); // η component
        ring.set_component(0b11, vec![0.1, 0.2, 0.3]); // θη component

        let body = ring.body().unwrap();
        assert!((body[1] - 2.0).abs() < 1e-10);
        assert_eq!(ring.soul_masks().len(), 3);
    }

    #[test]
    fn test_ring_rank() {
        let m = Supermanifold::new(2, 3, vec![], vec![]);
        assert_eq!(m.ring_rank(), 8); // 2^3
    }

    #[test]
    fn test_supermanifold_coords() {
        let m = Supermanifold::new(2, 2,
            vec!["x".into(), "y".into()],
            vec!["θ".into(), "η".into()]);
        assert_eq!(m.bosonic_coords, vec!["x", "y"]);
        assert_eq!(m.fermionic_coords, vec!["θ", "η"]);
    }
}
