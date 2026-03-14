"""Environment helpers backed by the shared physiology configuration."""

from __future__ import annotations

from core.config_loader import PhysiologyConfig

# This will be a lot more complex in the coming pushes
def build_environment_state(physiology_config: PhysiologyConfig) -> dict[str, float]:
    return {
        "pH": physiology_config.pH,
        "salt_concentration_mM": physiology_config.salt_concentration_mM,
        "temperature_K": physiology_config.temperature_K,
        "exposure_duration_s": physiology_config.exposure_duration_s,
    }
