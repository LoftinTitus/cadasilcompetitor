"""Multiobjective proposal ranking helpers."""

from __future__ import annotations

import math
from typing import Any, Iterable

from .features import FEATURE_NAMES, extract_feature_vector
from .surrogate import BootstrapSurrogate

MAXIMIZE_FIELDS: tuple[str, ...] = (
    "hs_affinity_reward",
    "hs_selectivity_reward",
    "transport_reward",
    "residence_reward",
)
MINIMIZE_FIELDS: tuple[str, ...] = (
    "anticoagulant_penalty",
    "off_target_penalty",
    "heparin_penalty",
    "developability_penalty",
)
DEFAULT_TARGET_FIELDS: tuple[str, ...] = MAXIMIZE_FIELDS + MINIMIZE_FIELDS


def pareto_frontier(
    rows: Iterable[dict[str, Any]],
    *,
    maximize_fields: tuple[str, ...] = MAXIMIZE_FIELDS,
    minimize_fields: tuple[str, ...] = MINIMIZE_FIELDS,
    tolerance: float = 0.0,
    max_points: int | None = None,
) -> list[dict[str, Any]]:
    rows_list = list(rows)
    frontier: list[dict[str, Any]] = []
    for candidate in rows_list:
        dominated = False
        for challenger in rows_list:
            if challenger is candidate:
                continue
            if _dominates(
                challenger,
                candidate,
                maximize_fields=maximize_fields,
                minimize_fields=minimize_fields,
                tolerance=tolerance,
            ):
                dominated = True
                break
        if not dominated:
            frontier.append(candidate)

    frontier.sort(
        key=lambda row: (
            float(row.get("composite_screen_score", row.get("predicted_composite_score", 0.0))),
            float(row.get("hs_selectivity_reward", row.get("pred_hs_selectivity_reward", 0.0))),
            -float(row.get("anticoagulant_penalty", row.get("pred_anticoagulant_penalty", 0.0))),
        ),
        reverse=True,
    )
    if max_points is not None and len(frontier) > max_points:
        frontier = frontier[:max_points]
    return frontier


def propose_candidates(
    candidate_pool_rows: list[dict[str, Any]],
    *,
    trained_surrogate: BootstrapSurrogate,
    reference_rows: list[dict[str, Any]],
    top_k: int = 50,
) -> list[dict[str, Any]]:
    seen_sequences = {str(row.get("sequence", "") or "") for row in reference_rows}
    proposals: list[dict[str, Any]] = []
    for candidate in candidate_pool_rows:
        sequence = str(candidate.get("sequence", "") or "")
        if not sequence or sequence in seen_sequences:
            continue
        warning_flags = _normalize_flag_list(candidate.get("warning_flags"))
        filter_flags = _normalize_flag_list(candidate.get("filter_flags"))
        if filter_flags:
            continue

        predictions = _clamp_prediction_dict(trained_surrogate.predict(candidate))
        novelty_score = _compute_novelty_score(candidate, reference_rows)
        acquisition_score = _compute_acquisition_score(
            predictions,
            novelty_score,
            warning_flags=warning_flags,
        )
        predicted_composite_score = 100.0 * (
            (0.35 * predictions["hs_affinity_reward"]["mean"])
            + (0.20 * predictions["hs_selectivity_reward"]["mean"])
            + (0.20 * predictions["transport_reward"]["mean"])
            + (0.10 * predictions["residence_reward"]["mean"])
            - (0.10 * predictions["anticoagulant_penalty"]["mean"])
            - (0.07 * predictions["off_target_penalty"]["mean"])
            - (0.03 * predictions["heparin_penalty"]["mean"])
            - (0.05 * predictions["developability_penalty"]["mean"])
        )
        proposals.append(
            {
                "candidate_id": candidate.get("candidate_id", ""),
                "sequence": sequence,
                "length": candidate.get("length", len(sequence)),
                "approx_net_charge": candidate.get("approx_net_charge", ""),
                "warning_flags": warning_flags,
                "filter_flags": filter_flags,
                "acquisition_score": round(acquisition_score, 6),
                "novelty_score": round(novelty_score, 6),
                "uncertainty_score": round(_mean_uncertainty(predictions), 6),
                "predicted_composite_score": round(predicted_composite_score, 6),
                "pred_hs_affinity_reward": round(predictions["hs_affinity_reward"]["mean"], 6),
                "pred_hs_selectivity_reward": round(
                    predictions["hs_selectivity_reward"]["mean"],
                    6,
                ),
                "pred_transport_reward": round(predictions["transport_reward"]["mean"], 6),
                "pred_residence_reward": round(predictions["residence_reward"]["mean"], 6),
                "pred_anticoagulant_penalty": round(
                    predictions["anticoagulant_penalty"]["mean"],
                    6,
                ),
                "pred_off_target_penalty": round(predictions["off_target_penalty"]["mean"], 6),
                "pred_heparin_penalty": round(predictions["heparin_penalty"]["mean"], 6),
                "pred_developability_penalty": round(
                    predictions["developability_penalty"]["mean"],
                    6,
                ),
                "std_hs_selectivity_reward": round(
                    predictions["hs_selectivity_reward"]["std"],
                    6,
                ),
                "std_transport_reward": round(predictions["transport_reward"]["std"], 6),
                "std_anticoagulant_penalty": round(
                    predictions["anticoagulant_penalty"]["std"],
                    6,
                ),
            }
        )

    proposals.sort(
        key=lambda row: (
            float(row["acquisition_score"]),
            float(row["predicted_composite_score"]),
            float(row["novelty_score"]),
        ),
        reverse=True,
    )
    for rank, row in enumerate(proposals[:top_k], start=1):
        row["proposal_rank"] = rank
        row["warning_flags"] = ";".join(row["warning_flags"])
        row["filter_flags"] = ";".join(row["filter_flags"])
    return proposals[:top_k]


def _compute_acquisition_score(
    predictions: dict[str, dict[str, float]],
    novelty_score: float,
    *,
    warning_flags: list[str],
) -> float:
    score = (
        0.28 * _optimistic_value(predictions["hs_affinity_reward"])
        + 0.27 * _optimistic_value(predictions["hs_selectivity_reward"])
        + 0.20 * _optimistic_value(predictions["transport_reward"])
        + 0.10 * _optimistic_value(predictions["residence_reward"])
        - 0.08 * _conservative_penalty(predictions["anticoagulant_penalty"])
        - 0.04 * _conservative_penalty(predictions["off_target_penalty"])
        - 0.03 * _conservative_penalty(predictions["heparin_penalty"])
        - 0.05 * _conservative_penalty(predictions["developability_penalty"])
        + 0.08 * novelty_score
    )
    if "cpp_like_uptake_risk" in warning_flags:
        score -= 0.18
    if "aggregation_risk" in warning_flags:
        score -= 0.15
    return score


def _optimistic_value(summary: dict[str, float]) -> float:
    return float(summary["mean"]) + (0.5 * float(summary["std"]))


def _conservative_penalty(summary: dict[str, float]) -> float:
    return max(0.0, float(summary["mean"]) - (0.25 * float(summary["std"])))


def _compute_novelty_score(
    candidate: dict[str, Any],
    reference_rows: list[dict[str, Any]],
) -> float:
    if not reference_rows:
        return 1.0

    candidate_vector = extract_feature_vector(candidate)
    reference_vectors = [extract_feature_vector(row) for row in reference_rows]

    feature_scales: list[float] = []
    for feature_index in range(len(FEATURE_NAMES)):
        values = [vector[feature_index] for vector in reference_vectors]
        mean_value = sum(values) / len(values)
        variance = sum((value - mean_value) ** 2 for value in values) / len(values)
        feature_scales.append(math.sqrt(variance) or 1.0)

    min_distance = min(
        math.sqrt(
            sum(
                (
                    (candidate_value - reference_value) / feature_scale
                ) ** 2
                for candidate_value, reference_value, feature_scale in zip(
                    candidate_vector,
                    reference_vector,
                    feature_scales,
                )
            )
        )
        for reference_vector in reference_vectors
    )
    return 1.0 - math.exp(-0.35 * min_distance)


def _mean_uncertainty(predictions: dict[str, dict[str, float]]) -> float:
    return sum(float(summary["std"]) for summary in predictions.values()) / len(predictions)


def _clamp_prediction_dict(
    predictions: dict[str, dict[str, float]],
) -> dict[str, dict[str, float]]:
    bounded: dict[str, dict[str, float]] = {}
    for field_name, summary in predictions.items():
        if field_name in MAXIMIZE_FIELDS or field_name in MINIMIZE_FIELDS:
            bounded[field_name] = {
                key: _clamp01(float(value))
                for key, value in summary.items()
            }
        else:
            bounded[field_name] = dict(summary)
    return bounded


def _dominates(
    challenger: dict[str, Any],
    candidate: dict[str, Any],
    *,
    maximize_fields: tuple[str, ...],
    minimize_fields: tuple[str, ...],
    tolerance: float,
) -> bool:
    at_least_one_strict = False
    for field_name in maximize_fields:
        challenger_value = float(challenger[field_name])
        candidate_value = float(candidate[field_name])
        if challenger_value < (candidate_value - tolerance):
            return False
        if challenger_value > (candidate_value + tolerance):
            at_least_one_strict = True

    for field_name in minimize_fields:
        challenger_value = float(challenger[field_name])
        candidate_value = float(candidate[field_name])
        if challenger_value > (candidate_value + tolerance):
            return False
        if challenger_value < (candidate_value - tolerance):
            at_least_one_strict = True

    return at_least_one_strict


def _clamp01(value: float) -> float:
    return max(0.0, min(value, 1.0))


def _normalize_flag_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(item) for item in value if str(item)]
    if isinstance(value, str):
        return [item for item in value.split(";") if item]
    return [str(value)]
