//! Z₂-graded elements and parity.

use serde::{Deserialize, Serialize};
use std::fmt;
use std::ops::{Add, Mul, Neg, Sub};

/// Parity (Z₂ grade): Even or Odd.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Parity {
    Even,
    Odd,
}

impl Parity {
    /// Z₂ addition (XOR).
    pub fn combine(self, other: Parity) -> Parity {
        match (self, other) {
            (Parity::Even, Parity::Even) => Parity::Even,
            (Parity::Even, Parity::Odd) => Parity::Odd,
            (Parity::Odd, Parity::Even) => Parity::Odd,
            (Parity::Odd, Parity::Odd) => Parity::Even,
        }
    }

    /// (-1)^|p| sign factor.
    pub fn sign(self) -> i64 {
        match self {
            Parity::Even => 1,
            Parity::Odd => -1,
        }
    }

    /// Grade (0 for even, 1 for odd).
    pub fn grade(self) -> u8 {
        match self {
            Parity::Even => 0,
            Parity::Odd => 1,
        }
    }
}

/// A Z₂-graded element with a scalar value and parity.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GradedElement {
    pub value: f64,
    pub parity: Parity,
    pub label: String,
}

impl GradedElement {
    pub fn even(value: f64, label: &str) -> Self {
        GradedElement { value, parity: Parity::Even, label: label.to_string() }
    }

    pub fn odd(value: f64, label: &str) -> Self {
        GradedElement { value, parity: Parity::Odd, label: label.to_string() }
    }

    pub fn is_even(&self) -> bool { self.parity == Parity::Even }
    pub fn is_odd(&self) -> bool { self.parity == Parity::Odd }

    /// The sign factor (-1)^|a|.
    pub fn sign(&self) -> i64 {
        self.parity.sign()
    }
}

impl fmt::Display for GradedElement {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let g = if self.is_even() { "₀" } else { "₁" };
        write!(f, "{}{}({})", self.label, g, self.value)
    }
}

impl Add for GradedElement {
    type Output = GradedElement;
    fn add(self, rhs: GradedElement) -> GradedElement {
        assert_eq!(self.parity, rhs.parity, "Cannot add elements of different parity");
        GradedElement {
            value: self.value + rhs.value,
            parity: self.parity,
            label: format!("({}+{})", self.label, rhs.label),
        }
    }
}

impl Sub for GradedElement {
    type Output = GradedElement;
    fn sub(self, rhs: GradedElement) -> GradedElement {
        assert_eq!(self.parity, rhs.parity, "Cannot subtract elements of different parity");
        GradedElement {
            value: self.value - rhs.value,
            parity: self.parity,
            label: format!("({}-{})", self.label, rhs.label),
        }
    }
}

impl Neg for GradedElement {
    type Output = GradedElement;
    fn neg(self) -> GradedElement {
        GradedElement { value: -self.value, parity: self.parity, label: format!("-{}", self.label) }
    }
}

impl Mul for GradedElement {
    type Output = GradedElement;
    fn mul(self, rhs: GradedElement) -> GradedElement {
        GradedElement {
            value: self.value * rhs.value,
            parity: self.parity.combine(rhs.parity),
            label: format!("{}·{}", self.label, rhs.label),
        }
    }
}

/// A general Z₂-graded vector space element: sum of even and odd components.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GradedVec {
    pub even_part: Vec<(f64, String)>,
    pub odd_part: Vec<(f64, String)>,
}

impl GradedVec {
    pub fn new() -> Self {
        GradedVec { even_part: vec![], odd_part: vec![] }
    }

    pub fn add_even(&mut self, coeff: f64, label: &str) {
        self.even_part.push((coeff, label.to_string()));
    }

    pub fn add_odd(&mut self, coeff: f64, label: &str) {
        self.odd_part.push((coeff, label.to_string()));
    }

    pub fn parity(&self) -> Option<Parity> {
        let has_even = !self.even_part.is_empty();
        let has_odd = !self.odd_part.is_empty();
        match (has_even, has_odd) {
            (true, false) => Some(Parity::Even),
            (false, true) => Some(Parity::Odd),
            (true, true) | (false, false) => None, // mixed or zero
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parity_combine() {
        assert_eq!(Parity::Even.combine(Parity::Even), Parity::Even);
        assert_eq!(Parity::Even.combine(Parity::Odd), Parity::Odd);
        assert_eq!(Parity::Odd.combine(Parity::Even), Parity::Odd);
        assert_eq!(Parity::Odd.combine(Parity::Odd), Parity::Even);
    }

    #[test]
    fn test_parity_sign() {
        assert_eq!(Parity::Even.sign(), 1);
        assert_eq!(Parity::Odd.sign(), -1);
    }

    #[test]
    fn test_graded_element_mul_parity() {
        let a = GradedElement::odd(2.0, "θ");
        let b = GradedElement::odd(3.0, "η");
        let c = a * b;
        assert_eq!(c.parity, Parity::Even);
        assert!((c.value - 6.0).abs() < 1e-10);
    }

    #[test]
    fn test_graded_element_add_same_parity() {
        let a = GradedElement::even(1.0, "x");
        let b = GradedElement::even(2.0, "y");
        let c = a + b;
        assert_eq!(c.parity, Parity::Even);
        assert!((c.value - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_graded_element_neg() {
        let a = GradedElement::odd(3.0, "θ");
        let neg_a = -a;
        assert!((neg_a.value - (-3.0)).abs() < 1e-10);
        assert_eq!(neg_a.parity, Parity::Odd);
    }

    #[test]
    fn test_grades() {
        assert_eq!(Parity::Even.grade(), 0);
        assert_eq!(Parity::Odd.grade(), 1);
    }

    #[test]
    fn test_graded_element_display() {
        let a = GradedElement::even(42.0, "x");
        let s = format!("{}", a);
        assert!(s.contains("x"));
    }
    #[test]
    #[should_panic]
    fn test_graded_element_add_different_parity_panics() {
        let a = GradedElement::even(1.0, "x");
        let b = GradedElement::odd(2.0, "θ");
        let _ = a + b;
    }
}
