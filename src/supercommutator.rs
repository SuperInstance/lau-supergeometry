//! Supercommutator: [a, b] = ab - (-1)^|a||b| ba

use crate::graded::{GradedElement, Parity};

/// Sign factor (-1)^(|a|*|b|) using grades directly.
fn super_sign(a: Parity, b: Parity) -> f64 {
    if a.grade() * b.grade() == 1 { -1.0 } else { 1.0 }
}

/// Compute the supercommutator [a, b] = ab - (-1)^|a||b| ba.
pub fn supercommutator(a: &GradedElement, b: &GradedElement) -> GradedElement {
    let ab = a.clone() * b.clone();
    let sign = super_sign(a.parity, b.parity);
    let ba = b.clone() * a.clone();
    GradedElement {
        value: ab.value - sign * ba.value,
        parity: ab.parity,
        label: format!("[{},{}]", a.label, b.label),
    }
}

/// Extension for matrix-like supercommutator on pairs (value, parity).
/// Returns (result_value, result_parity).
pub fn supercomm_scalar(
    (a_val, a_par): (f64, Parity),
    (b_val, b_par): (f64, Parity),
) -> (f64, Parity) {
    let ab = a_val * b_val;
    let ba = b_val * a_val;
    let sign = super_sign(a_par, b_par);
    let result = ab - sign * ba;
    let result_parity = a_par.combine(b_par);
    (result, result_parity)
}

/// Verify graded symmetry: [a,b] = -(-1)^|a||b| [b,a]
pub fn check_graded_symmetry(a: &GradedElement, b: &GradedElement) -> bool {
    let ab = supercommutator(a, b);
    let ba = supercommutator(b, a);
    let sign = super_sign(a.parity, b.parity);
    (ab.value - (-sign * ba.value)).abs() < 1e-10
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_supercommutator_even_even() {
        // [even, even] = ab - ba (ordinary commutator)
        let a = GradedElement::even(3.0, "x");
        let b = GradedElement::even(5.0, "y");
        let result = supercommutator(&a, &b);
        // For scalars: ab = ba, so commutator = 0
        assert!(result.value.abs() < 1e-10);
    }

    #[test]
    fn test_supercommutator_odd_odd() {
        // [θ, η] = θη - (-1)^(1*1) ηθ = θη + ηθ = 2θη (anticommutator)
        let a = GradedElement::odd(2.0, "θ");
        let b = GradedElement::odd(3.0, "η");
        let result = supercommutator(&a, &b);
        // θη = 6, ηθ = 6, sign = -1 => 6 - (-1)*6 = 12
        assert!((result.value - 12.0).abs() < 1e-10);
        assert_eq!(result.parity, Parity::Even);
    }

    #[test]
    fn test_supercommutator_even_odd() {
        // [x, θ] = xθ - (-1)^(0*1) θx = xθ - θx = 0 (scalars commute)
        let a = GradedElement::even(4.0, "x");
        let b = GradedElement::odd(2.0, "θ");
        let result = supercommutator(&a, &b);
        assert!(result.value.abs() < 1e-10);
    }

    #[test]
    fn test_graded_symmetry() {
        let a = GradedElement::odd(1.0, "θ");
        let b = GradedElement::odd(1.0, "η");
        assert!(check_graded_symmetry(&a, &b));
    }

    #[test]
    fn test_supercomm_scalar() {
        let r = supercomm_scalar((3.0, Parity::Odd), (2.0, Parity::Odd));
        assert!((r.0 - 12.0).abs() < 1e-10);
        assert_eq!(r.1, Parity::Even);
    }
}
