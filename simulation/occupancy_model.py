"""Occupancy model inputs derived from the shared physiology configuration."""

from __future__ import annotations

from typing import Any

from core.config_loader import PhysiologyConfig
from simulation.environment import build_environment_state


def prepare_occupancy_model_inputs(
    candidate: Any,
    physiology_config: PhysiologyConfig,
) -> dict[str, Any]:
    environment = build_environment_state(physiology_config)
    return {
        "candidate": candidate,
        "physiology": environment,
        "exposure_duration_s": physiology_config.exposure_duration_s,
    }
