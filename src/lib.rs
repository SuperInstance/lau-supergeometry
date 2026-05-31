#![deny(unsafe_code)]
//! # lau-supergeometry
//!
//! Supergeometry — the mathematical framework behind supersymmetry and graded algebra.
//!
//! Provides Z₂-graded (super) algebras, supercommutators, super Lie algebras,
//! supermanifolds, superforms, Berezin integration, supervector spaces with
//! supertrace/superdeterminant (Berezinian), and an application to agent state
//! spaces with fermionic (exclusive) and bosonic (inclusive) degrees of freedom.

pub mod graded;
pub mod supercommutator;
pub mod super_lie;
pub mod superalgebra;
pub mod exterior;
pub mod symmetric;
pub mod supermanifold;
pub mod superform;
pub mod berezin;
pub mod supervector;
pub mod agent_state;

pub use graded::*;
pub use supercommutator::*;
pub use super_lie::*;
pub use superalgebra::*;
pub use exterior::*;
pub use symmetric::*;
pub use supermanifold::*;
pub use superform::*;
pub use berezin::*;
pub use supervector::*;
pub use agent_state::*;
