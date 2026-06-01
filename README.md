# lau-supergeometry

**Supergeometry in pure Rust** — Z₂-graded algebras, supercommutators, Berezin integration, supermanifolds, supervector spaces with supertrace/superdeterminant (Berezinian), and an application to agent state spaces with fermionic/bosonic degrees of freedom.

---

## What This Does

This crate implements the algebraic foundations of **supergeometry** — the mathematical language behind supersymmetry, superstring theory, and graded algebra. It provides:

- **Z₂-graded elements** with even/odd parity and sign rules
- **Supercommutators** and graded symmetry checks
- **Super Lie algebras** with Jacobi identity verification
- **Super algebras** (associative, with multiplication tables)
- **Exterior algebras** Λ(V) with wedge products and Koszul signs
- **Symmetric algebras** S(V) (super-symmetric version)
- **Supermanifolds** (p|q) with coordinate rings C^∞(ℝ^p) ⊗ Λ(θ₁,...,θ_q)
- **Superforms** — differential forms on supermanifolds with exterior derivatives
- **Berezin integration** — integration over Grassmann variables (∫θ dθ = 1, ∫1 dθ = 0)
- **Supervector spaces** with block matrices, supertrace, and the **Berezinian** (superdeterminant)
- **Agent state spaces** — modeling multi-agent systems with fermionic exclusion constraints

All built on `nalgebra` for linear algebra and `serde` for serialization.

---

## Key Idea

In supergeometry, every object carries a **Z₂ grading** (even or odd). The fundamental sign rule is:

```
ab = (-1)^|a||b| ba
```

where |a| is the parity of a. Even elements commute normally; odd elements anticommute. This single rule generates:
- The supercommutator: `[a, b] = ab - (-1)^|a||b| ba`
- Exterior algebras (where all generators are odd → anticommuting)
- The Berezin integral (which is really differentiation on Grassmann variables)
- The supertrace: `str(M) = tr(A) - tr(D)` and superdeterminant (Berezinian): `Ber(M) = det(A) · det(D - CA⁻¹B)⁻¹`

---

## Install

Add to your `Cargo.toml`:

```toml
[dependencies]
lau-supergeometry = "0.1.0"
```

Requires Rust 2021 edition. Dependencies: `nalgebra`, `serde`, `num-traits`.

---

## Quick Start

```rust
use lau_supergeometry::*;

// Z₂-graded elements
let x = GradedElement::even(3.0, "x");
let theta = GradedElement::odd(1.0, "θ");
let eta = GradedElement::odd(2.0, "η");

// Odd × Odd = Even (with sign flip in commutator)
let product = theta.clone() * eta.clone();
assert_eq!(product.parity, Parity::Even);

// Supercommutator: [θ, η] = θη + ηθ = 2θη (anticommutator for odd elements)
let bracket = supercommutator(&theta, &eta);

// Berezin integral over 2 Grassmann variables
let integrator = BerezinIntegrator::new(2, vec![]);
let mut coeffs = std::collections::HashMap::new();
coeffs.insert(0b11, 5.0); // 5·θ₁θ₂
assert_eq!(integrator.integrate(&coeffs), 5.0); // ∫ 5·θ₁θ₂ dθ₂dθ₁ = 5

// Supermatrix with supertrace and Berezinian
let m = SuperMatrix::from_block_matrix(
    2, 1,
    &[vec![2.0, 0.0], vec![0.0, 3.0]],  // A (2×2)
    &[vec![0.0], vec![0.0]],              // B (2×1)
    &[vec![0.0, 0.0]],                    // C (1×2)
    &[vec![4.0]],                         // D (1×1)
);
assert_eq!(m.supertrace(), (2.0 + 3.0) - 4.0);  // tr(A) - tr(D) = 1
assert_eq!(m.berezinian(), Some(6.0 / 4.0));     // det(A)/det(D) = 1.5

// Agent state with fermionic exclusion
let mut sys = MultiAgentSystem::new();
let mut a1 = AgentState::new("agent-1");
a1.add_bosonic("position", 3.14);
a1.add_fermionic("leader", true);
let mut a2 = AgentState::new("agent-2");
a2.add_fermionic("leader", true);  // conflict!
sys.add_agent(a1);
sys.add_agent(a2);
assert!(sys.check_exclusion().is_err());  // Two agents can't both be leader
```

---

## API Reference

### `graded` — Z₂-Graded Elements

| Type / Function | Description |
|---|---|
| `Parity` | `Even` or `Odd` — the Z₂ grade |
| `Parity::combine(a, b)` | Z₂ addition (XOR) |
| `Parity::sign()` | Returns `(-1)^|p|`: `+1` for even, `-1` for odd |
| `GradedElement` | A labeled scalar with parity. Supports `+`, `-`, `×` (parity-aware) |
| `GradedVec` | General graded vector with separate even/odd components |

### `supercommutator` — Graded Brackets

| Function | Description |
|---|---|
| `supercommutator(a, b)` | `[a, b] = ab - (-1)^|a||b| ba` |
| `supercomm_scalar((val, par), ...)` | Scalar-only version |
| `check_graded_symmetry(a, b)` | Verifies `[a,b] = -(-1)^|a||b| [b,a]` |

### `super_lie` — Super Lie Algebras

| Function | Description |
|---|---|
| `super_jacobi_check(a, b, c)` | Checks the graded Jacobi identity |
| `SuperLieAlgebra` | Collection of generators with bracket and Jacobi verification |

### `superalgebra` — Z₂-Graded Associative Algebras

| Type / Method | Description |
|---|---|
| `SuperAlgebra::new(name, basis)` | Create algebra with labeled basis elements |
| `set_product(i, j, val, label)` | Define multiplication table |
| `multiply(i, j)` | Multiply basis elements |
| `is_supercommutative()` | Check `ab = (-1)^|a||b| ba` |
| `is_associative()` | Check `(ab)c = a(bc)` |

### `exterior` — Exterior Algebra Λ(V)

| Type / Method | Description |
|---|---|
| `ExteriorAlgebra::new(n, labels)` | Λ(V) on n generators |
| `dimension()` | Returns `2^n` |
| `form_dimension(k)` | Returns `C(n, k)` |
| `wedge(a_mask, b_mask)` | Wedge product (returns `None` if generators overlap → 0) |
| `koszul_sign(a_mask, b_mask)` | Sign from reordering |
| `signed_wedge(a, b)` | Full wedge with Koszul sign |

### `symmetric` — Super-Symmetric Algebra S(V)

| Method | Description |
|---|---|
| `sym_product(i, j)` | Super-symmetric product (nilpotent for odd self-products) |
| `is_supercommutative()` | Check graded commutativity |

### `supermanifold` — Supermanifolds (p|q)

| Type / Method | Description |
|---|---|
| `Supermanifold::new(p, q, ...)` | Dimension (p\|q) with bosonic/fermionic coordinates |
| `coordinate_ring()` | Returns `CoordinateRing` — functions f(x,θ) |
| `CoordinateRing::set_body(values)` | Set the θ-free part |
| `CoordinateRing::set_component(mask, values)` | Set a θ-monomial component |
| `ring_rank()` | Returns `2^q` |

### `superform` — Differential Forms on Supermanifolds

| Type / Method | Description |
|---|---|
| `SuperForm::new(p, q)` | Create form on (p\|q) manifold |
| `add_component(coeff_mask, diff_mask, value)` | Add a term |
| `exterior_derivative()` | Compute dω |
| `wedge(other)` | Wedge product of superforms |
| `volume_form(p, q)` | Standard volume form `dx₁∧...∧dx_p∧dθ₁∧...∧dθ_q` |

### `berezin` — Berezin Integration

| Type / Method | Description |
|---|---|
| `BerezinIntegrator::new(n, labels)` | Integrator over n Grassmann variables |
| `integrate(coefficients)` | Extracts coefficient of top form θ₁...θ_n |
| `integrate_over(coeffs, var_index)` | Iterated integration over one variable |
| `integrate_iterated(coeffs)` | Full iterated Berezin integral |

### `supervector` — Supertrace & Berezinian

| Type / Method | Description |
|---|---|
| `SuperVectorSpace::new(p, q)` | (p\|q) supervector space |
| `SuperMatrix::new(p, q)` | Block matrix `[A B; C D]` on a (p\|q) space |
| `SuperMatrix::from_block_matrix(...)` | Construct from A, B, C, D blocks |
| `supertrace()` | `str(M) = tr(A) - tr(D)` |
| `berezinian()` | `Ber(M) = det(A) / det(D - CA⁻¹B)` |
| `block_a/b/c/d()` | Extract individual blocks |

### `agent_state` — Multi-Agent Systems

| Type / Method | Description |
|---|---|
| `AgentState::new(id)` | Agent with bosonic/fermionic DOFs |
| `add_bosonic(name, value)` | Add continuous DOF (e.g., position) |
| `add_fermionic(name, occupied)` | Add binary DOF (e.g., role assignment) |
| `MultiAgentSystem::check_exclusion()` | Verify no two agents share a fermionic DOF |
| `MultiAgentSystem::berezinian()` | Berezinian of combined state matrix |

---

## How It Works

The crate is structured in layers:

1. **Foundation** (`graded`): `Parity` enum and `GradedElement` with full arithmetic (`Add`, `Sub`, `Mul`, `Neg`). Multiplication combines parities via XOR.

2. **Algebra** (`supercommutator`, `superalgebra`, `super_lie`): Graded brackets, associative algebras with multiplication tables, and super Lie algebras with Jacobi identity checks.

3. **Classical algebras** (`exterior`, `symmetric`): Exterior (Grassmann) and symmetric algebras as special cases of Z₂-graded algebras. Elements represented as bitmasks for efficiency.

4. **Geometry** (`supermanifold`, `superform`): Supermanifolds with coordinate rings, superforms with exterior derivatives and wedge products.

5. **Integration** (`berezin`): Berezin integral — picks out the coefficient of the top Grassmann monomial. Supports iterated integration and change of variables.

6. **Linear algebra** (`supervector`): Block matrix representation, supertrace, and the Berezinian (superdeterminant) computed via the quotient formula `det(A - BD⁻¹C) / det(D)`.

7. **Application** (`agent_state`): Models multi-agent systems where bosonic DOFs are shareable (position, velocity) and fermionic DOFs are exclusive (role assignment, task ownership). Exclusion constraints are checked via fermionic occupancy.

---

## The Math

### Z₂ Grading
Every element has parity |a| ∈ {0, 1}. The sign rule (-1)^|a||b| governs all exchanges.

### Supercommutator
```
[a, b] = ab - (-1)^|a||b| ba
```
- Two even elements: ordinary commutator
- Two odd elements: anticommutator
- Mixed: one commutes past the other

### Berezin Integral
For Grassmann variable θ with f = a + bθ:
```
∫ f dθ = b      (∫ θ dθ = 1, ∫ 1 dθ = 0)
```
This is formally "differentiation" — the Berezin integral coincides with the Grassmann derivative.

### Supertrace and Berezinian
For a block matrix M = `[A B; C D]`:
```
str(M) = tr(A) - tr(D)
Ber(M) = det(A) · det(D - CA⁻¹B)⁻¹
```
The Berezinian is the correct notion of determinant for supermatrices — it reduces to `det(A)/det(D)` for diagonal matrices and satisfies `Ber(MN) = Ber(M)·Ber(N)`.

### Supermanifolds
A supermanifold of dimension (p|q) has coordinate ring C^∞(ℝ^p) ⊗ Λ(θ₁,...,θ_q). Functions expand as:
```
f(x, θ) = f₀(x) + f₁(x)θ₁ + ... + f_{2^q-1}(x)θ₁θ₂...θ_q
```

---

## Tests

**58 tests** covering:
- Parity arithmetic and sign rules
- Supercommutator for all parity combinations
- Jacobi identity verification
- Exterior algebra wedge products and Koszul signs
- Berezin integration (single, iterated, multi-variable)
- Supertrace linearity and Berezinian for diagonal, identity, and off-diagonal matrices
- Agent exclusion constraints
- Roundtrip raise/lower index operations

Run with:
```bash
cargo test
```

---

## License

MIT
