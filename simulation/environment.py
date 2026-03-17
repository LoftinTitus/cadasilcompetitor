"""Environment helpers backed by the shared simulation configuration."""

from __future__ import annotations

from dataclasses import asdict
from typing import Any

from core.config_loader import (
    CandidateProperties,
    PhysiologyConfig,
    SimulationConfig,
    TargetConfig,
    validate_candidate_properties,
)


def build_environment_state(
    simulation_config: SimulationConfig | PhysiologyConfig,
) -> dict[str, Any]:
    if isinstance(simulation_config, SimulationConfig):
        return asdict(simulation_config)

    return {"physiology": asdict(simulation_config)}


def build_candidate_property_state(
    candidate_properties: CandidateProperties | None,
) -> dict[str, float] | None:
    if candidate_properties is None:
        return None

    return asdict(validate_candidate_properties(candidate_properties))


def build_target_context(target_config: TargetConfig) -> dict[str, str]:
    screening_focus_by_compartment = {
        "glycocalyx": "surface_occupancy_under_flow",
        "basement_membrane": "subendothelial_retention",
        "perivascular_ecm": "ecm_penetration_and_retention",
    }

    return {
        "compartment": target_config.compartment,
        "exposure_side": target_config.exposure_side,
        "biological_rationale_label": target_config.biological_rationale_label,
        "screening_focus": screening_focus_by_compartment[target_config.compartment],
    }
