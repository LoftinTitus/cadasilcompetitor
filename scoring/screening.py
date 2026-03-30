"""Candidate screening, scoring, and ranked CSV reporting."""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any, Iterable

from core.config_loader import (
    HSVariantConfig,
    HSVariantPanelConfig,
    SimulationConfig,
    load_hs_variant_panel,
    load_simulation_config,
)
from simulation.binding_simulation import run_binding_simulation, run_binding_simulation_panel
from simulation.candidate_property_estimation import estimate_candidate_properties

TIER_B_VARIANT_IDS: tuple[str, ...] = ("HS-dp6-mixed", "HS-dp5-AT-motif")
NON_SENTINEL_SULFATED_VARIANTS: tuple[str, ...] = (
    "HS-dp4-NS",
    "HS-dp4-NS-2S",
    "HS-dp4-NS-6S",
    "HS-dp6-mixed",
)
UNSULFATED_VARIANT_ID = "HS-dp4-NAc"
ANTICOAGULANT_SENTINEL_ID = "HS-dp5-AT-motif"

RANKED_CSV_FIELDS = [
    "rank",
    "candidate_id",
    "sequence",
    "length",
    "approx_net_charge",
    "basic_fraction",
    "hydrophobic_fraction",
    "warning_flags",
    "filter_flags",
    "screening_status",
    "composite_screen_score",
    "hs_affinity_reward",
    "hs_selectivity_reward",
    "transport_reward",
    "residence_reward",
    "anticoagulant_penalty",
    "developability_penalty",
    "tier_a_sulfated_peak_mean",
    "tier_a_unsulfated_peak",
    "tier_a_anticoagulant_peak",
    "tier_a_best_variant_id",
    "tier_a_best_variant_peak",
    "tier_b_mixed_peak_mean_occupancy",
    "tier_b_mixed_penetration_depth_um",
    "tier_b_anticoagulant_peak_mean_occupancy",
    "tier_b_anticoagulant_penetration_depth_um",
    "estimated_diffusion_coefficient_um2_s",
    "estimated_clearance_rate_per_s",
    "estimated_association_rate_M_inv_s",
    "estimated_dissociation_rate_per_s",
    "estimated_barrier_permeability_cm_s",
    "estimated_half_life_s",
    "estimated_enzymatic_degradation_rate_per_s",
    "estimated_spontaneous_degradation_rate_per_s",
    "risk_flags",
]


def screen_candidates_from_csv(
    candidate_csv_path: str | Path,
    *,
    output_csv_path: str | Path | None = None,
    simulation_config: SimulationConfig | None = None,
    hs_variant_panel: HSVariantPanelConfig | None = None,
    limit: int | None = None,
) -> list[dict[str, object]]:
    candidates = load_candidates_from_csv(candidate_csv_path, limit=limit)
    ranked_rows = rank_candidates(
        candidates,
        simulation_config=simulation_config,
        hs_variant_panel=hs_variant_panel,
    )
    if output_csv_path is not None:
        write_ranked_candidates_csv(output_csv_path, ranked_rows)
    return ranked_rows


def rank_candidates(
    candidates: Iterable[dict[str, Any]],
    *,
    simulation_config: SimulationConfig | None = None,
    hs_variant_panel: HSVariantPanelConfig | None = None,
) -> list[dict[str, object]]:
    resolved_simulation_config = simulation_config or load_simulation_config(
        "config/physiology.yaml"
    )
    resolved_panel = hs_variant_panel or load_hs_variant_panel("config/hs_variant_panel.yaml")

    scored_rows = [
        screen_candidate(
            candidate,
            simulation_config=resolved_simulation_config,
            hs_variant_panel=resolved_panel,
        )
        for candidate in candidates
    ]
    scored_rows.sort(
        key=lambda row: (
            _screening_status_priority(str(row["screening_status"])),
            float(row["composite_screen_score"]),
            float(row["hs_affinity_reward"]),
            float(row["transport_reward"]),
            -float(row["anticoagulant_penalty"]),
        ),
        reverse=True,
    )
    for rank, row in enumerate(scored_rows, start=1):
        row["rank"] = rank
    return scored_rows


def screen_candidate(
    candidate: dict[str, Any],
    *,
    simulation_config: SimulationConfig,
    hs_variant_panel: HSVariantPanelConfig,
) -> dict[str, object]:
    candidate_properties = estimate_candidate_properties(candidate)
    tier_a_panel = run_binding_simulation_panel(
        candidate,
        simulation_config,
        candidate_properties=candidate_properties,
        hs_variant_panel=hs_variant_panel,
        model_tier="A",
    )
    tier_a_result_by_variant = _map_results_by_variant_id(tier_a_panel["results"])

    tier_b_result_by_variant = {
        variant_id: run_binding_simulation(
            candidate,
            simulation_config,
            candidate_properties=candidate_properties,
            hs_variant=_variant_by_id(hs_variant_panel, variant_id),
            model_tier="B",
        )
        for variant_id in TIER_B_VARIANT_IDS
    }

    warning_flags = _normalize_flag_list(candidate.get("warning_flags"))
    filter_flags = _normalize_flag_list(candidate.get("filter_flags"))
    risk_flags = _build_risk_flags(
        warning_flags=warning_flags,
        filter_flags=filter_flags,
        tier_a_result_by_variant=tier_a_result_by_variant,
        tier_b_result_by_variant=tier_b_result_by_variant,
    )

    tier_a_sulfated_peak_mean = _mean(
        _tier_peak(tier_a_result_by_variant[variant_id])
        for variant_id in NON_SENTINEL_SULFATED_VARIANTS
    )
    tier_a_unsulfated_peak = _tier_peak(tier_a_result_by_variant[UNSULFATED_VARIANT_ID])
    tier_a_anticoagulant_peak = _tier_peak(
        tier_a_result_by_variant[ANTICOAGULANT_SENTINEL_ID]
    )
    tier_b_mixed_result = tier_b_result_by_variant["HS-dp6-mixed"]
    tier_b_anticoagulant_result = tier_b_result_by_variant[ANTICOAGULANT_SENTINEL_ID]

    hs_affinity_reward = _clamp01(tier_a_sulfated_peak_mean)
    hs_selectivity_reward = _clamp01(tier_a_sulfated_peak_mean - tier_a_unsulfated_peak)
    transport_reward = _clamp01(
        0.65 * float(tier_b_mixed_result["summary"]["best_case_max_mean_occupancy_fraction"])
        + 0.35
        * _safe_fraction(
            float(tier_b_mixed_result["summary"]["best_case_max_penetration_depth_um"]),
            float(tier_b_mixed_result["operating_point"]["summary"]["domain_length_um"]),
        )
    )
    residence_reward = _clamp01(
        math.log10(
            1.0
            + float(
                tier_a_result_by_variant["HS-dp6-mixed"]["model_parameters"]["residence_time_s"]
            )
        )
        / 4.0
    )
    anticoagulant_penalty = _clamp01(
        0.55 * tier_a_anticoagulant_peak
        + 0.45
        * float(tier_b_anticoagulant_result["summary"]["best_case_max_mean_occupancy_fraction"])
        + (
            0.20
            if tier_a_anticoagulant_peak >= _tier_peak(tier_a_result_by_variant["HS-dp6-mixed"])
            else 0.0
        )
    )
    developability_penalty = _compute_developability_penalty(
        warning_flags=warning_flags,
        filter_flags=filter_flags,
        candidate=candidate,
    )

    composite_screen_score = 100.0 * (
        (0.35 * hs_affinity_reward)
        + (0.20 * hs_selectivity_reward)
        + (0.20 * transport_reward)
        + (0.10 * residence_reward)
        - (0.10 * anticoagulant_penalty)
        - (0.05 * developability_penalty)
    )
    screening_status = _determine_screening_status(
        filter_flags=filter_flags,
        risk_flags=risk_flags,
        anticoagulant_penalty=anticoagulant_penalty,
    )

    tier_a_best_variant_id, tier_a_best_variant_peak = _best_variant_peak(
        tier_a_result_by_variant
    )
    candidate_properties_dict = tier_a_panel["results"][0]["candidate_properties"]

    return {
        "rank": 0,
        "candidate_id": candidate.get("candidate_id", ""),
        "sequence": candidate.get("sequence", ""),
        "length": int(float(candidate.get("length", 0) or 0)),
        "approx_net_charge": float(candidate.get("approx_net_charge", 0) or 0),
        "basic_fraction": float(candidate.get("basic_fraction", 0) or 0),
        "hydrophobic_fraction": float(candidate.get("hydrophobic_fraction", 0) or 0),
        "warning_flags": ";".join(warning_flags),
        "filter_flags": ";".join(filter_flags),
        "screening_status": screening_status,
        "composite_screen_score": round(composite_screen_score, 6),
        "hs_affinity_reward": round(hs_affinity_reward, 6),
        "hs_selectivity_reward": round(hs_selectivity_reward, 6),
        "transport_reward": round(transport_reward, 6),
        "residence_reward": round(residence_reward, 6),
        "anticoagulant_penalty": round(anticoagulant_penalty, 6),
        "developability_penalty": round(developability_penalty, 6),
        "tier_a_sulfated_peak_mean": round(tier_a_sulfated_peak_mean, 6),
        "tier_a_unsulfated_peak": round(tier_a_unsulfated_peak, 6),
        "tier_a_anticoagulant_peak": round(tier_a_anticoagulant_peak, 6),
        "tier_a_best_variant_id": tier_a_best_variant_id,
        "tier_a_best_variant_peak": round(tier_a_best_variant_peak, 6),
        "tier_b_mixed_peak_mean_occupancy": round(
            float(tier_b_mixed_result["summary"]["best_case_max_mean_occupancy_fraction"]),
            6,
        ),
        "tier_b_mixed_penetration_depth_um": round(
            float(tier_b_mixed_result["summary"]["best_case_max_penetration_depth_um"]),
            6,
        ),
        "tier_b_anticoagulant_peak_mean_occupancy": round(
            float(
                tier_b_anticoagulant_result["summary"][
                    "best_case_max_mean_occupancy_fraction"
                ]
            ),
            6,
        ),
        "tier_b_anticoagulant_penetration_depth_um": round(
            float(
                tier_b_anticoagulant_result["summary"][
                    "best_case_max_penetration_depth_um"
                ]
            ),
            6,
        ),
        "estimated_diffusion_coefficient_um2_s": round(
            float(candidate_properties_dict["diffusion_coefficient_um2_s"]),
            6,
        ),
        "estimated_clearance_rate_per_s": round(
            float(candidate_properties_dict["clearance_rate_per_s"]),
            6,
        ),
        "estimated_association_rate_M_inv_s": round(
            float(candidate_properties_dict["association_rate_M_inv_s"]),
            6,
        ),
        "estimated_dissociation_rate_per_s": round(
            float(candidate_properties_dict["dissociation_rate_per_s"]),
            6,
        ),
        "estimated_barrier_permeability_cm_s": round(
            float(candidate_properties_dict["barrier_permeability_cm_s"]),
            12,
        ),
        "estimated_half_life_s": round(float(candidate_properties_dict["half_life_s"]), 6),
        "estimated_enzymatic_degradation_rate_per_s": round(
            float(candidate_properties_dict["enzymatic_degradation_rate_per_s"]),
            9,
        ),
        "estimated_spontaneous_degradation_rate_per_s": round(
            float(candidate_properties_dict["spontaneous_degradation_rate_per_s"]),
            9,
        ),
        "risk_flags": ";".join(risk_flags),
    }


def load_candidates_from_csv(
    path: str | Path,
    *,
    limit: int | None = None,
) -> list[dict[str, Any]]:
    csv_path = Path(path)
    with csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows: list[dict[str, Any]] = []
        for index, row in enumerate(reader, start=1):
            rows.append(_normalize_candidate_row(row))
            if limit is not None and index >= limit:
                break
        return rows


def write_ranked_candidates_csv(
    output_path: str | Path,
    ranked_rows: Iterable[dict[str, object]],
) -> None:
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=RANKED_CSV_FIELDS)
        writer.writeheader()
        for row in ranked_rows:
            writer.writerow({field: row.get(field, "") for field in RANKED_CSV_FIELDS})


def _normalize_candidate_row(row: dict[str, str]) -> dict[str, Any]:
    normalized = dict(row)
    for field_name in (
        "length",
        "approx_net_charge",
        "basic_fraction",
        "hydrophobic_fraction",
        "acidic_fraction",
        "polar_fraction",
        "longest_basic_run",
        "longest_hydrophobic_run",
    ):
        value = normalized.get(field_name, "")
        if value == "":
            continue
        if field_name in {"length", "longest_basic_run", "longest_hydrophobic_run"}:
            normalized[field_name] = int(float(value))
        else:
            normalized[field_name] = float(value)

    for field_name in ("warning_flags", "filter_flags"):
        normalized[field_name] = _normalize_flag_list(normalized.get(field_name))
    return normalized


def _normalize_flag_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(item) for item in value if str(item)]
    if isinstance(value, str):
        if not value.strip():
            return []
        return [item for item in value.split(";") if item]
    return [str(value)]


def _map_results_by_variant_id(
    panel_results: list[dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    return {
        result["hs_variant"]["variant_id"]: result
        for result in panel_results
        if result.get("hs_variant") is not None
    }


def _variant_by_id(
    panel: HSVariantPanelConfig,
    variant_id: str,
) -> HSVariantConfig:
    for variant in panel.variants:
        if variant.variant_id == variant_id:
            return variant
    raise KeyError(f"Could not find HS variant '{variant_id}' in panel '{panel.panel_id}'.")


def _tier_peak(result: dict[str, Any]) -> float:
    summary = result["summary"]
    if "selection_peak_occupancy_fraction" in summary:
        return float(summary["selection_peak_occupancy_fraction"])
    if "best_case_max_mean_occupancy_fraction" in summary:
        return float(summary["best_case_max_mean_occupancy_fraction"])
    return float(summary.get("best_case_max_occupancy_fraction", 0.0))


def _best_variant_peak(
    tier_a_result_by_variant: dict[str, dict[str, Any]]
) -> tuple[str, float]:
    best_variant_id = ""
    best_peak = -1.0
    for variant_id, result in tier_a_result_by_variant.items():
        peak = _tier_peak(result)
        if peak > best_peak:
            best_peak = peak
            best_variant_id = variant_id
    return best_variant_id, best_peak


def _compute_developability_penalty(
    *,
    warning_flags: list[str],
    filter_flags: list[str],
    candidate: dict[str, Any],
) -> float:
    penalty = 0.0
    if "aggregation_risk" in warning_flags:
        penalty += 0.35
    if "cpp_like_uptake_risk" in warning_flags:
        penalty += 0.25
    if filter_flags:
        penalty += 0.6
    if float(candidate.get("hydrophobic_fraction", 0.0) or 0.0) > 0.35:
        penalty += 0.15
    if float(candidate.get("basic_fraction", 0.0) or 0.0) > 0.55:
        penalty += 0.15
    return _clamp01(penalty)


def _build_risk_flags(
    *,
    warning_flags: list[str],
    filter_flags: list[str],
    tier_a_result_by_variant: dict[str, dict[str, Any]],
    tier_b_result_by_variant: dict[str, dict[str, Any]],
) -> list[str]:
    risk_flags: list[str] = []
    if filter_flags:
        risk_flags.append("filter_fail")
    if "aggregation_risk" in warning_flags:
        risk_flags.append("aggregation_risk")
    if "cpp_like_uptake_risk" in warning_flags:
        risk_flags.append("uptake_risk")

    sentinel_peak = _tier_peak(tier_a_result_by_variant[ANTICOAGULANT_SENTINEL_ID])
    mixed_peak = _tier_peak(tier_a_result_by_variant["HS-dp6-mixed"])
    if sentinel_peak >= mixed_peak:
        risk_flags.append("anticoagulant_bias")

    tier_b_mixed_penetration = float(
        tier_b_result_by_variant["HS-dp6-mixed"]["summary"]["best_case_max_penetration_depth_um"]
    )
    tier_b_mixed_domain = float(
        tier_b_result_by_variant["HS-dp6-mixed"]["operating_point"]["summary"]["domain_length_um"]
    )
    if _safe_fraction(tier_b_mixed_penetration, tier_b_mixed_domain) < 0.25:
        risk_flags.append("poor_penetration")

    return risk_flags


def _mean(values: Iterable[float]) -> float:
    values_list = list(values)
    if not values_list:
        return 0.0
    return sum(values_list) / len(values_list)


def _safe_fraction(numerator: float, denominator: float) -> float:
    if denominator <= 0.0:
        return 0.0
    return numerator / denominator


def _clamp01(value: float) -> float:
    return max(0.0, min(value, 1.0))


def _determine_screening_status(
    *,
    filter_flags: list[str],
    risk_flags: list[str],
    anticoagulant_penalty: float,
) -> str:
    if filter_flags:
        return "reject"
    if risk_flags or anticoagulant_penalty >= 0.9:
        return "review"
    return "pass"


def _screening_status_priority(status: str) -> int:
    if status == "pass":
        return 2
    if status == "review":
        return 1
    return 0
