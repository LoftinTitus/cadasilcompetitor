"""Heuristic candidate-property estimation for intro screening.

This module converts a peptide candidate into ``CandidateProperties`` so the
transport models can run without manually specifying per-candidate kinetics and
transport constants. The estimates are intentionally simple and should be
treated as screening priors, not calibrated biophysical truth.
"""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any

from core.config_loader import CandidateProperties
from models.property_surrogate import (
    DEFAULT_PROPERTY_ESTIMATOR_PATH,
    PropertySurrogate,
    try_load_property_surrogate,
)
from peptide.metadata import compute_sequence_metadata

_PROPERTY_SURROGATE_CACHE: dict[Path, PropertySurrogate | None] = {}


def estimate_candidate_properties(
    candidate: Any,
    *,
    learned_estimator: PropertySurrogate | None = None,
    property_estimator_path: str | Path | None = DEFAULT_PROPERTY_ESTIMATOR_PATH,
) -> CandidateProperties:
    """Estimate transport and kinetic properties from peptide sequence features."""

    heuristic_properties = _estimate_candidate_properties_heuristic(candidate)
    resolved_estimator = learned_estimator or _load_cached_property_estimator(property_estimator_path)
    if resolved_estimator is None:
        return heuristic_properties

    learned_properties = resolved_estimator.predict_candidate_properties(
        _candidate_to_mapping(candidate)
    )
    return _blend_candidate_properties(
        heuristic_properties=heuristic_properties,
        learned_properties=learned_properties,
        learned_weight=0.7,
    )


def _estimate_candidate_properties_heuristic(candidate: Any) -> CandidateProperties:
    """Estimate transport and kinetic properties from peptide sequence features."""

    metadata = _extract_candidate_metadata(candidate)
    warning_flags = _extract_flag_list(candidate, "warning_flags")
    filter_flags = _extract_flag_list(candidate, "filter_flags")

    length = int(metadata["length"])
    approx_net_charge = float(metadata["approx_net_charge"])
    basic_fraction = float(metadata["basic_fraction"])
    hydrophobic_fraction = float(metadata["hydrophobic_fraction"])
    acidic_fraction = float(metadata["acidic_fraction"])
    polar_fraction = float(metadata["polar_fraction"])
    longest_basic_run = int(metadata["longest_basic_run"])
    longest_hydrophobic_run = int(metadata["longest_hydrophobic_run"])
    dominant_residue_fraction = float(metadata["dominant_residue_fraction"])

    cpp_like_warning = "cpp_like_uptake_risk" in warning_flags
    aggregation_warning = "aggregation_risk" in warning_flags
    low_complexity_flag = "low_complexity" in filter_flags or bool(
        metadata["is_low_complexity"]
    )

    diffusion_coefficient_um2_s = _clamp(
        125.0
        - (3.0 * max(length - 8, 0))
        - (18.0 * hydrophobic_fraction)
        - (8.0 * basic_fraction)
        - (6.0 * max(approx_net_charge - 6.0, 0.0)),
        20.0,
        140.0,
    )

    barrier_permeability_cm_s = _clamp(
        2.5e-7
        + (3.0e-7 * basic_fraction)
        + (2.0e-7 * hydrophobic_fraction)
        + (1.5e-7 if cpp_like_warning else 0.0)
        - (1.0e-7 * acidic_fraction)
        - (2.0e-8 * max(length - 12, 0)),
        1.0e-8,
        5.0e-6,
    )

    clearance_rate_per_s = _clamp(
        4.0e-4
        + (9.0e-5 * max(approx_net_charge - 3.0, 0.0))
        + (6.0e-4 * barrier_permeability_cm_s / 1.0e-6)
        + (2.0e-4 if cpp_like_warning else 0.0)
        + (1.0e-4 if aggregation_warning else 0.0),
        1.0e-4,
        0.01,
    )

    association_rate_M_inv_s = _clamp(
        5.0e4
        + (1.8e5 * basic_fraction)
        + (1.5e4 * approx_net_charge)
        + (2.5e4 * max(longest_basic_run - 2, 0))
        - (4.0e4 * acidic_fraction),
        2.0e4,
        8.0e5,
    )

    dissociation_rate_per_s = _clamp(
        0.7
        - (0.05 * approx_net_charge)
        - (0.55 * basic_fraction)
        + (0.25 * acidic_fraction)
        + (0.12 * hydrophobic_fraction)
        + (0.08 if aggregation_warning else 0.0),
        0.03,
        1.5,
    )

    half_life_s = _clamp(
        1800.0
        + (350.0 * length)
        + (1200.0 * hydrophobic_fraction)
        - (450.0 * basic_fraction)
        - (600.0 if low_complexity_flag else 0.0)
        - (450.0 if cpp_like_warning else 0.0),
        900.0,
        43200.0,
    )

    enzymatic_degradation_rate_per_s = _clamp(
        5.0e-5
        + (4.0e-5 * polar_fraction)
        + (3.0e-5 * basic_fraction)
        + (2.0e-5 * max(length - 10, 0) / 10.0)
        + (4.0e-5 if low_complexity_flag else 0.0),
        1.0e-5,
        0.003,
    )

    spontaneous_degradation_rate_per_s = _clamp(
        1.0e-5
        + (2.0e-5 * acidic_fraction)
        + (1.0e-5 * hydrophobic_fraction)
        + (5.0e-6 * max(longest_hydrophobic_run - 2, 0)),
        1.0e-6,
        5.0e-4,
    )

    return CandidateProperties(
        diffusion_coefficient_um2_s=diffusion_coefficient_um2_s,
        clearance_rate_per_s=clearance_rate_per_s,
        association_rate_M_inv_s=association_rate_M_inv_s,
        dissociation_rate_per_s=dissociation_rate_per_s,
        barrier_permeability_cm_s=barrier_permeability_cm_s,
        half_life_s=half_life_s,
        enzymatic_degradation_rate_per_s=enzymatic_degradation_rate_per_s,
        spontaneous_degradation_rate_per_s=spontaneous_degradation_rate_per_s,
    )


def estimate_candidate_properties_with_diagnostics(candidate: Any) -> dict[str, Any]:
    """Return estimated properties plus the features used to derive them."""

    metadata = _extract_candidate_metadata(candidate)
    heuristic_properties = _estimate_candidate_properties_heuristic(candidate)
    learned_estimator = _load_cached_property_estimator(DEFAULT_PROPERTY_ESTIMATOR_PATH)
    properties = estimate_candidate_properties(candidate, learned_estimator=learned_estimator)
    return {
        "candidate_properties": properties,
        "candidate_properties_dict": asdict(properties),
        "heuristic_candidate_properties_dict": asdict(heuristic_properties),
        "property_estimation_mode": "learned_blend" if learned_estimator is not None else "heuristic",
        "sequence_metadata": metadata,
        "warning_flags": _extract_flag_list(candidate, "warning_flags"),
        "filter_flags": _extract_flag_list(candidate, "filter_flags"),
    }


def _extract_candidate_metadata(candidate: Any) -> dict[str, object]:
    sequence = _read_candidate_value(candidate, "sequence")
    if not isinstance(sequence, str) or not sequence.strip():
        raise ValueError(
            "Candidate property estimation requires a candidate with a non-empty 'sequence'."
        )

    metadata = compute_sequence_metadata(sequence)
    for field_name in (
        "approx_net_charge",
        "basic_fraction",
        "hydrophobic_fraction",
        "acidic_fraction",
        "polar_fraction",
        "longest_basic_run",
        "longest_hydrophobic_run",
    ):
        value = _read_candidate_value(candidate, field_name)
        if value is not None:
            metadata[field_name] = value
    return metadata


def _extract_flag_list(candidate: Any, field_name: str) -> list[str]:
    value = _read_candidate_value(candidate, field_name)
    if value is None:
        return []
    if isinstance(value, list):
        return [str(item) for item in value]
    if isinstance(value, str):
        if not value.strip():
            return []
        return [item for item in value.split(";") if item]
    return [str(value)]


def _read_candidate_value(candidate: Any, field_name: str) -> Any:
    if candidate is None:
        return None
    if isinstance(candidate, dict):
        return candidate.get(field_name)
    if is_dataclass(candidate):
        return getattr(candidate, field_name, None)
    return getattr(candidate, field_name, None)


def _clamp(value: float, minimum: float, maximum: float) -> float:
    return max(minimum, min(value, maximum))


def _candidate_to_mapping(candidate: Any) -> dict[str, Any]:
    if isinstance(candidate, dict):
        return dict(candidate)
    metadata = _extract_candidate_metadata(candidate)
    return {
        "sequence": _read_candidate_value(candidate, "sequence"),
        **metadata,
        "warning_flags": _extract_flag_list(candidate, "warning_flags"),
        "filter_flags": _extract_flag_list(candidate, "filter_flags"),
    }


def _blend_candidate_properties(
    *,
    heuristic_properties: CandidateProperties,
    learned_properties: CandidateProperties,
    learned_weight: float,
) -> CandidateProperties:
    heuristic_dict = asdict(heuristic_properties)
    learned_dict = asdict(learned_properties)
    blended = {
        field_name: ((learned_weight * float(learned_dict[field_name])) + ((1.0 - learned_weight) * float(heuristic_dict[field_name])))
        for field_name in heuristic_dict
    }
    return CandidateProperties(**blended)


def _load_cached_property_estimator(
    property_estimator_path: str | Path | None,
) -> PropertySurrogate | None:
    if property_estimator_path is None:
        return None
    estimator_path = Path(property_estimator_path)
    if estimator_path in _PROPERTY_SURROGATE_CACHE:
        return _PROPERTY_SURROGATE_CACHE[estimator_path]
    estimator = try_load_property_surrogate(estimator_path)
    _PROPERTY_SURROGATE_CACHE[estimator_path] = estimator
    return estimator
