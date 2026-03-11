"""Structured sequence filters for peptide candidates."""

from __future__ import annotations

from dataclasses import dataclass

# Importing configs for min/max props, what aa does what, etc.
from .config import (
    MAX_BASIC_FRACTION,
    MAX_CONSECUTIVE_BASIC_RESIDUES,
    MAX_CONSECUTIVE_HYDROPHOBIC_RESIDUES,
    MAX_HYDROPHOBIC_FRACTION,
    MAX_REPEATED_IDENTICAL_RESIDUE_RUN,
    MAX_ALLOWED_NET_CHARGE,
    MIN_ALLOWED_NET_CHARGE,
    PEPTIDE_MAX_LENGTH,
    PEPTIDE_MIN_LENGTH,
)
from .metadata import compute_sequence_metadata


@dataclass(slots=True)
class FilterResult:
    filter_name: str
    passed: bool
    severity: str
    message: str


def _make_result(filter_name: str, passed: bool, severity: str, message: str) -> FilterResult:
    return FilterResult(
        filter_name=filter_name,
        passed=passed,
        severity=severity,
        message=message,
    )


def _range_filter(
    metadata: dict[str, object],
    *,
    field_name: str,
    metric_name: str,
    min_value: int | float | None = None,
    max_value: int | float | None = None,
    severity: str = "error",
) -> FilterResult:
    value = float(metadata[field_name])
    if min_value is not None and value < min_value:
        return _make_result(
            f"{field_name}_out_of_range",
            False,
            severity,
            f"{metric_name} {value:g} is below minimum {min_value:g}.",
        )
    if max_value is not None and value > max_value:
        return _make_result(
            f"{field_name}_out_of_range",
            False,
            severity,
            f"{metric_name} {value:g} exceeds maximum {max_value:g}.",
        )
    return _make_result(
        f"{field_name}_out_of_range",
        True,
        severity,
        f"{metric_name} {value:g} is within the allowed range.",
    )


def _max_filter(
    metadata: dict[str, object],
    *,
    field_name: str,
    metric_name: str,
    max_value: int | float,
    severity: str = "error",
) -> FilterResult:
    return _range_filter(
        metadata,
        field_name=field_name,
        metric_name=metric_name,
        min_value=None,
        max_value=max_value,
        severity=severity,
    )


def filter_low_complexity(metadata: dict[str, object]) -> FilterResult:
    passed = not bool(metadata["is_low_complexity"])
    return _make_result(
        "low_complexity",
        passed,
        "error",
        "Sequence triggers the low-complexity heuristic.",
    )


def filter_aggregation_risk(metadata: dict[str, object]) -> FilterResult:
    hydrophobic_fraction = float(metadata["hydrophobic_fraction"])
    dominant_fraction = float(metadata["dominant_residue_fraction"])
    longest_hydrophobic_run = int(metadata["longest_hydrophobic_run"])
    triggered = (
        hydrophobic_fraction >= 0.35 and longest_hydrophobic_run >= 3
    ) or (hydrophobic_fraction >= 0.30 and dominant_fraction >= 0.40)
    return _make_result(
        "aggregation_risk",
        not triggered,
        "warning",
        "Sequence shows a simple aggregation-risk composition pattern.",
    )


def filter_cpp_like_warning(metadata: dict[str, object]) -> FilterResult:
    basic_fraction = float(metadata["basic_fraction"])
    charge = float(metadata["approx_net_charge"])
    longest_basic_run = int(metadata["longest_basic_run"])
    triggered = (basic_fraction >= 0.45 and charge >= 6) or longest_basic_run >= 4
    return _make_result(
        "cpp_like_uptake_risk",
        not triggered,
        "warning",
        "Sequence matches a simple CPP-like uptake-risk heuristic.",
    )


def apply_all_filters(
    sequence: str,
    metadata: dict[str, object] | None = None,
) -> dict[str, object]:
    computed_metadata = metadata or compute_sequence_metadata(sequence)
    results = [
        _range_filter(
            computed_metadata,
            field_name="length",
            metric_name="Length",
            min_value=PEPTIDE_MIN_LENGTH,
            max_value=PEPTIDE_MAX_LENGTH,
        ),
        _range_filter(
            computed_metadata,
            field_name="approx_net_charge",
            metric_name="Approximate net charge",
            min_value=MIN_ALLOWED_NET_CHARGE,
            max_value=MAX_ALLOWED_NET_CHARGE,
        ),
        _max_filter(
            computed_metadata,
            field_name="basic_fraction",
            metric_name="Basic fraction",
            max_value=MAX_BASIC_FRACTION,
        ),
        _max_filter(
            computed_metadata,
            field_name="hydrophobic_fraction",
            metric_name="Hydrophobic fraction",
            max_value=MAX_HYDROPHOBIC_FRACTION,
        ),
        _max_filter(
            computed_metadata,
            field_name="longest_basic_run",
            metric_name="Longest basic run",
            max_value=MAX_CONSECUTIVE_BASIC_RESIDUES,
        ),
        _max_filter(
            computed_metadata,
            field_name="longest_hydrophobic_run",
            metric_name="Longest hydrophobic run",
            max_value=MAX_CONSECUTIVE_HYDROPHOBIC_RESIDUES,
        ),
        _max_filter(
            computed_metadata,
            field_name="longest_identical_run",
            metric_name="Longest identical residue run",
            max_value=MAX_REPEATED_IDENTICAL_RESIDUE_RUN,
        ),
        filter_low_complexity(computed_metadata),
        filter_aggregation_risk(computed_metadata),
        filter_cpp_like_warning(computed_metadata),
    ]

    filter_flags = [
        result.filter_name
        for result in results
        if result.severity == "error" and not result.passed
    ]
    warning_flags = [
        result.filter_name
        for result in results
        if result.severity == "warning" and not result.passed
    ]

    return {
        "passed_filters": not filter_flags,
        "filter_flags": filter_flags,
        "warning_flags": warning_flags,
        "filter_results": results,
    }
