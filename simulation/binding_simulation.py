"""Binding simulation entry points that accept shared physiology configuration.

The caller is responsible for loading the YAML once at the top level and then
passing the resulting ``PhysiologyConfig`` into simulation functions:

    physiology = load_physiology_config("config/physiology.yaml")
    run_binding_simulation(candidate, physiology)
"""

from __future__ import annotations

from typing import Any

from core.config_loader import PhysiologyConfig
from simulation.environment import build_environment_state
from simulation.occupancy_model import prepare_occupancy_model_inputs


def run_binding_simulation(
    candidate: Any,
    physiology_config: PhysiologyConfig,
) -> dict[str, Any]:
    return {
        "candidate": candidate,
        "environment": build_environment_state(physiology_config),
        "occupancy_model_inputs": prepare_occupancy_model_inputs(
            candidate,
            physiology_config,
        ),
    }
