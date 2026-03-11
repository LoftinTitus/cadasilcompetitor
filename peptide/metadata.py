"""Pure sequence metadata helpers."""

from __future__ import annotations

from collections import Counter

from .config import (
    ACIDIC_RESIDUES,
    ALLOWED_RESIDUES,
    BASIC_RESIDUES,
    HYDROPHOBIC_RESIDUES,
    POLAR_RESIDUES,
)

_CHARGE_TABLE = {"K": 1, "R": 1, "D": -1, "E": -1}


def normalize_sequence(sequence: str) -> str:
    normalized = "".join(sequence.split()).upper()
    invalid_chars = sorted(set(normalized) - set(ALLOWED_RESIDUES))
    if invalid_chars:
        raise ValueError(
            f"Sequence '{sequence}' contains unsupported residues: {', '.join(invalid_chars)}"
        )
    return normalized


def compute_residue_counts(sequence: str) -> dict[str, int]:
    normalized = normalize_sequence(sequence)
    counts = Counter(normalized)
    return {residue: counts[residue] for residue in sorted(counts)}


def compute_residue_fractions(sequence: str) -> dict[str, float]:
    normalized = normalize_sequence(sequence)
    length = len(normalized)
    counts = compute_residue_counts(normalized)
    return {residue: count / length for residue, count in counts.items()}


def compute_approx_net_charge(sequence: str) -> int:
    normalized = normalize_sequence(sequence)
    return sum(_CHARGE_TABLE.get(residue, 0) for residue in normalized)


def positions_of_residues(sequence: str, residues: set[str] | frozenset[str]) -> list[int]:
    normalized = normalize_sequence(sequence)
    return [index for index, residue in enumerate(normalized, start=1) if residue in residues]


def spacing_between_positions(positions: list[int]) -> list[int]:
    if len(positions) < 2:
        return []
    return [right - left for left, right in zip(positions, positions[1:])]


def longest_consecutive_run(sequence: str, residues: set[str] | frozenset[str]) -> int:
    normalized = normalize_sequence(sequence)
    longest = 0
    current = 0
    for residue in normalized:
        if residue in residues:
            current += 1
            longest = max(longest, current)
        else:
            current = 0
    return longest


def longest_identical_residue_run(sequence: str) -> int:
    normalized = normalize_sequence(sequence)
    longest = 0
    current = 0
    previous = ""
    for residue in normalized:
        if residue == previous:
            current += 1
        else:
            previous = residue
            current = 1
        longest = max(longest, current)
    return longest


def low_complexity_heuristic(sequence: str) -> bool:
    normalized = normalize_sequence(sequence)
    counts = compute_residue_counts(normalized)
    dominant_fraction = max(counts.values()) / len(normalized)
    unique_residues = len(counts)
    return dominant_fraction >= 0.5 or unique_residues <= 3 or longest_identical_residue_run(normalized) >= 4


def compute_sequence_metadata(sequence: str) -> dict[str, object]:
    normalized = normalize_sequence(sequence)
    length = len(normalized)
    residue_counts = compute_residue_counts(normalized)
    residue_fractions = compute_residue_fractions(normalized)
    basic_positions = positions_of_residues(normalized, BASIC_RESIDUES)

    hydrophobic_count = sum(residue_counts.get(residue, 0) for residue in HYDROPHOBIC_RESIDUES)
    polar_count = sum(residue_counts.get(residue, 0) for residue in POLAR_RESIDUES)
    acidic_count = sum(residue_counts.get(residue, 0) for residue in ACIDIC_RESIDUES)
    basic_count = sum(residue_counts.get(residue, 0) for residue in BASIC_RESIDUES)

    dominant_fraction = max(residue_counts.values()) / length

    return {
        "sequence": normalized,
        "length": length,
        "approx_net_charge": compute_approx_net_charge(normalized),
        "residue_counts": residue_counts,
        "residue_fractions": residue_fractions,
        "basic_fraction": basic_count / length,
        "hydrophobic_fraction": hydrophobic_count / length,
        "polar_fraction": polar_count / length,
        "acidic_fraction": acidic_count / length,
        "longest_basic_run": longest_consecutive_run(normalized, BASIC_RESIDUES),
        "longest_hydrophobic_run": longest_consecutive_run(normalized, HYDROPHOBIC_RESIDUES),
        "longest_identical_run": longest_identical_residue_run(normalized),
        "basic_positions": basic_positions,
        "charge_spacing": spacing_between_positions(basic_positions),
        "is_low_complexity": low_complexity_heuristic(normalized),
        "dominant_residue_fraction": dominant_fraction,
        "unique_residue_count": len(residue_counts),
    }
