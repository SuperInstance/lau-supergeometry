# lau-supergeometry

Supergeometry extends geometry to include anti-commuting (Grassmann) coordinates. The result: a mathematical framework where bosonic (inclusive) and fermionic (exclusive) degrees of freedom live in the same space. It's the language of supersymmetry, and it naturally models systems with both shared and exclusive resources.

## The math in 60 seconds

A **supermanifold** has both even (commuting) and odd (anti-commuting) coordinates. The coordinate ring is C∞(ℝᵖ) ⊗ Λ(θ₁…θq) — smooth functions tensored with an exterior algebra. Key structures:

- **Z₂-grading:** every element is either Even or Odd
- **Supercommutator:** [a,b] = ab - (-1)^|a||b|ba — graded antisymmetry
- **Super Lie algebras:** satisfy the super Jacobi identity
- **Berezin integral:** ∫θ dθ = 1, ∫1 dθ = 0 — integration on Grassmann variables
- **Supertrace:** str(A) = tr(A) - tr(D) for a block matrix
- **Berezinian:** the superdeterminant det(A-D·C⁻¹·B)·det(C) — generalizes the Jacobian

References: Manin, *Gauge Field Theory and Complex Geometry* (1988); DeWitt, *Supermanifolds* (1984)

## Quick start

```rust
use lau_supergeometry::{GradedElement, SuperLieAlgebra, BerezinIntegral, SuperMatrix};

// Create a Z₂-graded element
let even = GradedElement::even(vec![1.0, 2.0]);
let odd = GradedElement::odd(vec![0.5]);

// Supercommutator: [even, odd] = even·odd - (-1)^(0·1)·odd·even = 0
let bracket = even.supercommutator(&odd);

// Super Lie algebra with generators
let sla = SuperLieAlgebra::from_generators(&["H", "Q"]);
assert!(sla.verify_super_jacobi());

// Berezin integral: ∫θ dθ = 1
let integral = BerezinIntegral::of_monomial(1); // θ¹ → 1
assert_eq!(integral, 1.0);

// Supertrace and Berezinian (superdeterminant)
let matrix = SuperMatrix::from_blocks(&a, &b, &c, &d);
let strace = matrix.supertrace();
let ber = matrix.berezinian(); // superdeterminant
```

## Key types

| Type | What it is |
|------|-----------|
| `GradedElement` | A Z₂-graded vector — Even (bosonic) or Odd (fermionic) |
| `Supercommutator` | The graded bracket [a,b] = ab - (-1)^|a||b|ba |
| `SuperLieAlgebra` | A Lie algebra satisfying the super Jacobi identity |
| `Supermanifold` | (p|q)-dimensional space with even + odd coordinates |
| `BerezinIntegral` | Integration over Grassmann (anti-commuting) variables |
| `SuperMatrix` | Block matrix with supertrace and Berezinian |

## Contributing

[Open an issue](https://github.com/SuperInstance/lau-supergeometry/issues) or PR.
