"""Structured sequence filters for peptide candidates."""

from __future__ import annotations

from dataclasses import dataclass

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


def filter_sequence_too_short(metadata: dict[str, object]) -> FilterResult:
    length = int(metadata["length"])
    passed = length >= PEPTIDE_MIN_LENGTH
    return _make_result(
        "sequence_too_short",
        passed,
        "error",
        f"Length {length} is below minimum {PEPTIDE_MIN_LENGTH}.",
    )


def filter_sequence_too_long(metadata: dict[str, object]) -> FilterResult:
    length = int(metadata["length"])
    passed = length <= PEPTIDE_MAX_LENGTH
    return _make_result(
        "sequence_too_long",
        passed,
        "error",
        f"Length {length} exceeds maximum {PEPTIDE_MAX_LENGTH}.",
    )


def filter_net_charge_too_high(metadata: dict[str, object]) -> FilterResult:
    charge = float(metadata["approx_net_charge"])
    passed = charge <= MAX_ALLOWED_NET_CHARGE
    return _make_result(
        "net_charge_too_high",
        passed,
        "error",
        f"Approximate net charge {charge:g} exceeds maximum {MAX_ALLOWED_NET_CHARGE}.",
    )


def filter_net_charge_too_low(metadata: dict[str, object]) -> FilterResult:
    charge = float(metadata["approx_net_charge"])
    passed = charge >= MIN_ALLOWED_NET_CHARGE
    return _make_result(
        "net_charge_too_low",
        passed,
        "error",
        f"Approximate net charge {charge:g} is below minimum {MIN_ALLOWED_NET_CHARGE}.",
    )


def filter_basic_fraction_too_high(metadata: dict[str, object]) -> FilterResult:
    fraction = float(metadata["basic_fraction"])
    passed = fraction <= MAX_BASIC_FRACTION
    return _make_result(
        "basic_fraction_too_high",
        passed,
        "error",
        f"Basic fraction {fraction:.2f} exceeds maximum {MAX_BASIC_FRACTION:.2f}.",
    )


def filter_hydrophobic_fraction_too_high(metadata: dict[str, object]) -> FilterResult:
    fraction = float(metadata["hydrophobic_fraction"])
    passed = fraction <= MAX_HYDROPHOBIC_FRACTION
    return _make_result(
        "hydrophobic_fraction_too_high",
        passed,
        "error",
        f"Hydrophobic fraction {fraction:.2f} exceeds maximum {MAX_HYDROPHOBIC_FRACTION:.2f}.",
    )


def filter_basic_run_too_long(metadata: dict[str, object]) -> FilterResult:
    longest_run = int(metadata["longest_basic_run"])
    passed = longest_run <= MAX_CONSECUTIVE_BASIC_RESIDUES
    return _make_result(
        "basic_run_too_long",
        passed,
        "error",
        f"Longest basic run {longest_run} exceeds maximum {MAX_CONSECUTIVE_BASIC_RESIDUES}.",
    )


def filter_hydrophobic_run_too_long(metadata: dict[str, object]) -> FilterResult:
    longest_run = int(metadata["longest_hydrophobic_run"])
    passed = longest_run <= MAX_CONSECUTIVE_HYDROPHOBIC_RESIDUES
    return _make_result(
        "hydrophobic_run_too_long",
        passed,
        "error",
        f"Longest hydrophobic run {longest_run} exceeds maximum {MAX_CONSECUTIVE_HYDROPHOBIC_RESIDUES}.",
    )


def filter_repeated_identical_residue_run(metadata: dict[str, object]) -> FilterResult:
    longest_run = int(metadata["longest_identical_run"])
    passed = longest_run <= MAX_REPEATED_IDENTICAL_RESIDUE_RUN
    return _make_result(
        "repeated_identical_residue_run_too_long",
        passed,
        "error",
        f"Longest identical residue run {longest_run} exceeds maximum {MAX_REPEATED_IDENTICAL_RESIDUE_RUN}.",
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
        filter_sequence_too_short(computed_metadata),
        filter_sequence_too_long(computed_metadata),
        filter_net_charge_too_high(computed_metadata),
        filter_net_charge_too_low(computed_metadata),
        filter_basic_fraction_too_high(computed_metadata),
        filter_hydrophobic_fraction_too_high(computed_metadata),
        filter_basic_run_too_long(computed_metadata),
        filter_hydrophobic_run_too_long(computed_metadata),
        filter_repeated_identical_residue_run(computed_metadata),
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
