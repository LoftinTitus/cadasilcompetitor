"""Feature extraction for surrogate models."""

from __future__ import annotations

from typing import Any

from peptide.metadata import compute_sequence_metadata

TRACKED_RESIDUES: tuple[str, ...] = ("K", "R", "D", "E", "G", "Q", "S", "A", "W", "Y", "H", "P")
FEATURE_NAMES: tuple[str, ...] = (
    "length",
    "approx_net_charge",
    "basic_fraction",
    "hydrophobic_fraction",
    "acidic_fraction",
    "polar_fraction",
    "longest_basic_run",
    "longest_hydrophobic_run",
    "longest_identical_run",
    "dominant_residue_fraction",
    "unique_residue_count",
    "basic_density",
    "charge_asymmetry",
    "kr_ratio",
    "hydrophobic_run_density",
) + tuple(f"frac_{residue}" for residue in TRACKED_RESIDUES)


def extract_feature_map(candidate: dict[str, Any]) -> dict[str, float]:
    sequence = str(candidate.get("sequence", "") or "")
    metadata = compute_sequence_metadata(sequence)
    residue_fractions = metadata["residue_fractions"]

    lys_fraction = float(residue_fractions.get("K", 0.0))
    arg_fraction = float(residue_fractions.get("R", 0.0))
    length = float(metadata["length"])
    basic_fraction = float(metadata["basic_fraction"])
    longest_basic_run = float(metadata["longest_basic_run"])
    longest_hydrophobic_run = float(metadata["longest_hydrophobic_run"])

    feature_map = {
        "length": length,
        "approx_net_charge": float(metadata["approx_net_charge"]),
        "basic_fraction": basic_fraction,
        "hydrophobic_fraction": float(metadata["hydrophobic_fraction"]),
        "acidic_fraction": float(metadata["acidic_fraction"]),
        "polar_fraction": float(metadata["polar_fraction"]),
        "longest_basic_run": longest_basic_run,
        "longest_hydrophobic_run": longest_hydrophobic_run,
        "longest_identical_run": float(metadata["longest_identical_run"]),
        "dominant_residue_fraction": float(metadata["dominant_residue_fraction"]),
        "unique_residue_count": float(metadata["unique_residue_count"]),
        "basic_density": basic_fraction * max(length, 1.0),
        "charge_asymmetry": abs(lys_fraction - arg_fraction),
        "kr_ratio": lys_fraction / max(arg_fraction, 1e-6),
        "hydrophobic_run_density": longest_hydrophobic_run / max(length, 1.0),
    }
    for residue in TRACKED_RESIDUES:
        feature_map[f"frac_{residue}"] = float(residue_fractions.get(residue, 0.0))
    return feature_map


def extract_feature_vector(candidate: dict[str, Any]) -> list[float]:
    feature_map = extract_feature_map(candidate)
    return [feature_map[name] for name in FEATURE_NAMES]

