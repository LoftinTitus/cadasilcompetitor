"""Environment helpers backed by the shared simulation configuration."""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
from typing import Any

from core.config_loader import (
    CandidateProperties,
    HSVariantConfig,
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


def build_candidate_context(candidate: Any) -> dict[str, Any]:
    sequence = _read_candidate_value(candidate, "sequence")
    candidate_id = _read_candidate_value(candidate, "candidate_id") or _read_candidate_value(
        candidate,
        "id",
    )

    context: dict[str, Any] = {
        "candidate_id": candidate_id,
        "sequence": sequence,
        "length": len(sequence) if isinstance(sequence, str) else _read_candidate_value(candidate, "length"),
        "label": candidate_id or sequence or type(candidate).__name__,
    }

    approx_net_charge = _read_candidate_value(candidate, "approx_net_charge")
    if approx_net_charge is not None:
        context["approx_net_charge"] = approx_net_charge

    return {key: value for key, value in context.items() if value is not None}


def build_hs_variant_context(hs_variant: HSVariantConfig | None) -> dict[str, Any] | None:
    if hs_variant is None:
        return None

    sulfation_pattern = dict(hs_variant.sulfation_pattern)
    chain_length = dict(hs_variant.chain_length)
    simulation_annotations = dict(hs_variant.simulation_annotations)

    return {
        "variant_id": hs_variant.variant_id,
        "display_name": hs_variant.display_name,
        "gag_class": hs_variant.gag_class,
        "panel_role": hs_variant.panel_role,
        "degree_of_polymerization": chain_length.get("degree_of_polymerization"),
        "total_sulfates": sulfation_pattern.get("total_sulfates"),
        "contains_3_o_sulfation": sulfation_pattern.get("contains_3_o_sulfation"),
        "expected_relative_affinity": simulation_annotations.get("expected_relative_affinity"),
    }


def _read_candidate_value(candidate: Any, field_name: str) -> Any:
    if candidate is None:
        return None
    if isinstance(candidate, dict):
        return candidate.get(field_name)
    if is_dataclass(candidate):
        return getattr(candidate, field_name, None)
    return getattr(candidate, field_name, None)
