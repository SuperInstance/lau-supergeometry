//! Application: Agent state spaces with fermionic (exclusive) and bosonic (inclusive)
//! degrees of freedom.

use crate::graded::{GradedElement, GradedVec, Parity};
use crate::supervector::{SuperMatrix, SuperVectorSpace};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// A degree of freedom in an agent's state space.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum DegreeOfFreedom {
    /// Bosonic: can take any real value, multiple agents can share the same value.
    Bosonic { name: String, value: f64 },
    /// Fermionic: binary/exclusive state, only one agent can occupy.
    Fermionic { name: String, occupied: bool },
}

/// An agent's state in a super-vector space.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AgentState {
    pub agent_id: String,
    pub bosonic_dofs: Vec<(String, f64)>,
    pub fermionic_dofs: Vec<(String, bool)>,
}

impl AgentState {
    pub fn new(agent_id: &str) -> Self {
        AgentState {
            agent_id: agent_id.to_string(),
            bosonic_dofs: vec![],
            fermionic_dofs: vec![],
        }
    }

    /// Add a bosonic DOF (e.g., position, velocity).
    pub fn add_bosonic(&mut self, name: &str, value: f64) {
        self.bosonic_dofs.push((name.to_string(), value));
    }

    /// Add a fermionic DOF (e.g., task assigned, role occupied).
    pub fn add_fermionic(&mut self, name: &str, occupied: bool) {
        self.fermionic_dofs.push((name.to_string(), occupied));
    }

    /// Dimension of the agent's state space: (bosonic | fermionic).
    pub fn dimension(&self) -> (usize, usize) {
        (self.bosonic_dofs.len(), self.fermionic_dofs.len())
    }

    /// Get the super-vector space for this agent.
    pub fn super_vector_space(&self) -> SuperVectorSpace {
        let (p, q) = self.dimension();
        SuperVectorSpace::new(p, q)
    }

    /// Encode the state as a vector (bosonic values + fermionic indicators).
    pub fn to_vector(&self) -> Vec<f64> {
        let mut v: Vec<f64> = self.bosonic_dofs.iter().map(|(_, v)| *v).collect();
        v.extend(self.fermionic_dofs.iter().map(|(_, occ)| if *occ { 1.0 } else { 0.0 }));
        v
    }
}

/// A multi-agent system with fermionic exclusion constraints.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MultiAgentSystem {
    pub agents: Vec<AgentState>,
}

impl MultiAgentSystem {
    pub fn new() -> Self {
        MultiAgentSystem { agents: vec![] }
    }

    pub fn add_agent(&mut self, agent: AgentState) {
        self.agents.push(agent);
    }

    /// Check fermionic exclusion: no two agents can occupy the same fermionic DOF.
    pub fn check_exclusion(&self) -> Result<(), Vec<(String, String, String)>> {
        let mut violations = vec![];
        let mut occupied: HashMap<String, (String, bool)> = HashMap::new();

        for agent in &self.agents {
            for (name, occ) in &agent.fermionic_dofs {
                if *occ {
                    if let Some((other_id, _)) = occupied.get(name) {
                        violations.push((
                            agent.agent_id.clone(),
                            other_id.clone(),
                            name.clone(),
                        ));
                    } else {
                        occupied.insert(name.clone(), (agent.agent_id.clone(), true));
                    }
                }
            }
        }

        if violations.is_empty() {
            Ok(())
        } else {
            Err(violations)
        }
    }

    /// Compute the total state dimension across all agents.
    pub fn total_dimension(&self) -> (usize, usize) {
        let p: usize = self.agents.iter().map(|a| a.bosonic_dofs.len()).sum();
        let q: usize = self.agents.iter().map(|a| a.fermionic_dofs.len()).sum();
        (p, q)
    }

    /// Construct the combined super-matrix representation.
    pub fn to_super_matrix(&self) -> SuperMatrix {
        let (p, q) = self.total_dimension();
        let n = p + q;
        // Diagonal matrix with bosonic values and fermionic indicators
        let mut m = SuperMatrix::new(p, q);
        let mut row = 0;
        for agent in &self.agents {
            for (_, v) in &agent.bosonic_dofs {
                m.data[row][row] = *v;
                row += 1;
            }
        }
        for agent in &self.agents {
            for (_, occ) in &agent.fermionic_dofs {
                m.data[row][row] = if *occ { 1.0 } else { 0.0 };
                row += 1;
            }
        }
        m
    }

    /// Compute the Berezinian of the combined state matrix.
    pub fn berezinian(&self) -> Option<f64> {
        self.to_super_matrix().berezinian()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_agent_state_creation() {
        let mut s = AgentState::new("agent-1");
        s.add_bosonic("position", 3.14);
        s.add_fermionic("leader", true);
        assert_eq!(s.dimension(), (1, 1));
        assert_eq!(s.to_vector(), vec![3.14, 1.0]);
    }

    #[test]
    fn test_multi_agent_exclusion_ok() {
        let mut sys = MultiAgentSystem::new();
        let mut a1 = AgentState::new("a1");
        a1.add_fermionic("leader", true);
        let mut a2 = AgentState::new("a2");
        a2.add_fermionic("leader", false);
        sys.add_agent(a1);
        sys.add_agent(a2);
        assert!(sys.check_exclusion().is_ok());
    }

    #[test]
    fn test_multi_agent_exclusion_violation() {
        let mut sys = MultiAgentSystem::new();
        let mut a1 = AgentState::new("a1");
        a1.add_fermionic("leader", true);
        let mut a2 = AgentState::new("a2");
        a2.add_fermionic("leader", true);
        sys.add_agent(a1);
        sys.add_agent(a2);
        let result = sys.check_exclusion();
        assert!(result.is_err());
        let violations = result.unwrap_err();
        assert_eq!(violations.len(), 1);
        assert_eq!(violations[0].2, "leader");
    }

    #[test]
    fn test_multi_agent_total_dimension() {
        let mut sys = MultiAgentSystem::new();
        let mut a1 = AgentState::new("a1");
        a1.add_bosonic("x", 1.0);
        a1.add_bosonic("y", 2.0);
        a1.add_fermionic("role", true);
        let mut a2 = AgentState::new("a2");
        a2.add_bosonic("x", 3.0);
        a2.add_fermionic("role", false);
        sys.add_agent(a1);
        sys.add_agent(a2);
        assert_eq!(sys.total_dimension(), (3, 2));
    }

    #[test]
    fn test_agent_berezinian() {
        let mut sys = MultiAgentSystem::new();
        let mut a1 = AgentState::new("a1");
        a1.add_bosonic("x", 2.0);
        a1.add_fermionic("role", true);
        sys.add_agent(a1);
        // Ber = 2.0 / 1.0 = 2.0
        let ber = sys.berezinian().unwrap();
        assert!((ber - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_agent_super_vector_space() {
        let mut s = AgentState::new("a");
        s.add_bosonic("x", 1.0);
        s.add_bosonic("y", 2.0);
        s.add_fermionic("f1", true);
        s.add_fermionic("f2", false);
        let sv = s.super_vector_space();
        assert_eq!(sv.p, 2);
        assert_eq!(sv.q, 2);
        assert_eq!(sv.total_dimension(), 4);
    }
}
