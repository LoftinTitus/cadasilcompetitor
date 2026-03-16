"""Occupancy model inputs derived from the shared simulation configuration."""

from __future__ import annotations

from typing import Any

from core.config_loader import CandidateProperties, PhysiologyConfig, SimulationConfig
from simulation.environment import build_candidate_property_state, build_environment_state


def prepare_occupancy_model_inputs(
    candidate: Any,
    simulation_config: SimulationConfig | PhysiologyConfig,
    candidate_properties: CandidateProperties | None = None,
) -> dict[str, Any]:
    physiology = (
        simulation_config.physiology
        if isinstance(simulation_config, SimulationConfig)
        else simulation_config
    )
    return {
        "candidate": candidate,
        "environment": build_environment_state(simulation_config),
        "candidate_properties": build_candidate_property_state(candidate_properties),
        "exposure_duration_s": physiology.exposure_duration_s,
    }
