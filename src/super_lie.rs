//! Super Lie algebra: supercommutator satisfying super Jacobi identity.

use crate::graded::{GradedElement, Parity};
use crate::supercommutator::supercommutator;

/// Super Jacobi identity check:
/// (-1)^|a||c| [a,[b,c]] + (-1)^|b||a| [b,[c,a]] + (-1)^|c||b| [c,[a,b]] = 0
pub fn super_jacobi_check(
    a: &GradedElement,
    b: &GradedElement,
    c: &GradedElement,
) -> bool {
    let bc = supercommutator(b, c);
    let a_bc = supercommutator(a, &bc);
    let ca = supercommutator(c, a);
    let b_ca = supercommutator(b, &ca);
    let ab = supercommutator(a, b);
    let c_ab = supercommutator(c, &ab);

    let sign1 = if a.parity.grade() * c.parity.grade() % 2 != 0 { -1.0 } else { 1.0 };
    let sign2 = if b.parity.grade() * a.parity.grade() % 2 != 0 { -1.0 } else { 1.0 };
    let sign3 = if c.parity.grade() * b.parity.grade() % 2 != 0 { -1.0 } else { 1.0 };

    let total = sign1 * a_bc.value + sign2 * b_ca.value + sign3 * c_ab.value;
    total.abs() < 1e-10
}

/// A super Lie algebra defined by a set of generators with structure constants.
#[derive(Debug, Clone)]
pub struct SuperLieAlgebra {
    pub name: String,
    pub generators: Vec<GradedElement>,
}

impl SuperLieAlgebra {
    pub fn new(name: &str, generators: Vec<GradedElement>) -> Self {
        SuperLieAlgebra { name: name.to_string(), generators }
    }

    /// Compute supercommutator of generators i and j.
    pub fn bracket(&self, i: usize, j: usize) -> GradedElement {
        supercommutator(&self.generators[i], &self.generators[j])
    }

    /// Verify super Jacobi for all triples of generators.
    pub fn verify_jacobi(&self) -> bool {
        let n = self.generators.len();
        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    if !super_jacobi_check(&self.generators[i], &self.generators[j], &self.generators[k]) {
                        return false;
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
    fn test_jacobi_even_elements() {
        let a = GradedElement::even(1.0, "H");
        let b = GradedElement::even(2.0, "X");
        let c = GradedElement::even(3.0, "Y");
        assert!(super_jacobi_check(&a, &b, &c));
    }

    #[test]
    fn test_jacobi_odd_elements() {
        let a = GradedElement::odd(1.0, "Q₁");
        let b = GradedElement::odd(1.0, "Q₂");
        let c = GradedElement::even(1.0, "H");
        assert!(super_jacobi_check(&a, &b, &c));
    }

    #[test]
    fn test_jacobi_mixed() {
        let a = GradedElement::odd(2.0, "Q");
        let b = GradedElement::even(3.0, "P");
        let c = GradedElement::odd(1.0, "R");
        assert!(super_jacobi_check(&a, &b, &c));
    }

    #[test]
    fn test_super_lie_algebra_jacobi() {
        let gens = vec![
            GradedElement::even(1.0, "H"),
            GradedElement::odd(1.0, "Q"),
            GradedElement::odd(1.0, "Q̄"),
        ];
        let sla = SuperLieAlgebra::new("sqm", gens);
        assert!(sla.verify_jacobi());
    }
}
