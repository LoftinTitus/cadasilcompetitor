"""Binding simulation entry points that accept shared simulation configuration.

The caller is responsible for loading the YAML once at the top level and then
passing shared config and per-candidate computed properties into simulation
functions:

    config = load_simulation_config("config/physiology.yaml")
    candidate_properties = CandidateProperties(...)
    run_binding_simulation(candidate, config, candidate_properties)
"""

from __future__ import annotations

from typing import Any

from core.config_loader import CandidateProperties, PhysiologyConfig, SimulationConfig
from simulation.environment import build_candidate_property_state, build_environment_state
from simulation.occupancy_model import prepare_occupancy_model_inputs


def run_binding_simulation(
    candidate: Any,
    simulation_config: SimulationConfig | PhysiologyConfig,
    candidate_properties: CandidateProperties | None = None,
) -> dict[str, Any]:
    return {
        "candidate": candidate,
        "environment": build_environment_state(simulation_config),
        "candidate_properties": build_candidate_property_state(candidate_properties),
        "occupancy_model_inputs": prepare_occupancy_model_inputs(
            candidate,
            simulation_config,
            candidate_properties,
        ),
    }
