"""Feature extraction for surrogate models."""

from __future__ import annotations

from typing import Any

from peptide.metadata import compute_sequence_metadata

TRACKED_RESIDUES: tuple[str, ...] = ("K", "R", "D", "E", "G", "Q", "S", "A", "W", "Y", "H", "P")
TRACKED_MOTIFS: tuple[str, ...] = (
    "KK",
    "KR",
    "RK",
    "RR",
    "KXK",
    "RXR",
    "BBXB",
    "BXXB",
)
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
    "n_terminal_basic_fraction_3",
    "c_terminal_basic_fraction_3",
    "n_terminal_hydrophobic_fraction_3",
    "c_terminal_hydrophobic_fraction_3",
    "terminal_charge_asymmetry_3",
    "max_basic_window_fraction_5",
    "max_hydrophobic_window_fraction_5",
    "max_positive_charge_window_5",
    "mean_charge_transition_rate",
    "basic_pair_density",
) + tuple(f"frac_{residue}" for residue in TRACKED_RESIDUES)
FEATURE_NAMES = FEATURE_NAMES + tuple(f"motif_{motif}" for motif in TRACKED_MOTIFS)


def extract_feature_map(candidate: dict[str, Any]) -> dict[str, float]:
    sequence = str(candidate.get("sequence", "") or "")
    metadata = compute_sequence_metadata(sequence)
    residue_fractions = metadata["residue_fractions"]

    lys_fraction = float(residue_fractions.get("K", 0.0))
    arg_fraction = float(residue_fractions.get("R", 0.0))
    length = float(metadata["length"])
    length_int = int(length)
    basic_fraction = float(metadata["basic_fraction"])
    longest_basic_run = float(metadata["longest_basic_run"])
    longest_hydrophobic_run = float(metadata["longest_hydrophobic_run"])
    n_terminal_3 = sequence[:3]
    c_terminal_3 = sequence[-3:]
    charge_classes = [_charge_class(residue) for residue in sequence]

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
        "n_terminal_basic_fraction_3": _residue_class_fraction(n_terminal_3, {"K", "R"}),
        "c_terminal_basic_fraction_3": _residue_class_fraction(c_terminal_3, {"K", "R"}),
        "n_terminal_hydrophobic_fraction_3": _residue_class_fraction(
            n_terminal_3,
            {"A", "F", "I", "L", "M", "V", "W", "Y"},
        ),
        "c_terminal_hydrophobic_fraction_3": _residue_class_fraction(
            c_terminal_3,
            {"A", "F", "I", "L", "M", "V", "W", "Y"},
        ),
        "terminal_charge_asymmetry_3": abs(
            sum(_charge_class(residue) for residue in n_terminal_3)
            - sum(_charge_class(residue) for residue in c_terminal_3)
        )
        / max(min(length_int, 3), 1),
        "max_basic_window_fraction_5": _max_window_fraction(sequence, {"K", "R"}, 5),
        "max_hydrophobic_window_fraction_5": _max_window_fraction(
            sequence,
            {"A", "F", "I", "L", "M", "V", "W", "Y"},
            5,
        ),
        "max_positive_charge_window_5": _max_positive_charge_window(sequence, 5),
        "mean_charge_transition_rate": _charge_transition_rate(charge_classes),
        "basic_pair_density": _count_basic_pairs(sequence) / max(length - 1.0, 1.0),
    }
    for residue in TRACKED_RESIDUES:
        feature_map[f"frac_{residue}"] = float(residue_fractions.get(residue, 0.0))
    for motif in TRACKED_MOTIFS:
        feature_map[f"motif_{motif}"] = _motif_density(sequence, motif)
    return feature_map


def extract_feature_vector(candidate: dict[str, Any]) -> list[float]:
    feature_map = extract_feature_map(candidate)
    return [feature_map[name] for name in FEATURE_NAMES]


def _residue_class_fraction(sequence: str, residues: set[str]) -> float:
    if not sequence:
        return 0.0
    return sum(1 for residue in sequence if residue in residues) / len(sequence)


def _max_window_fraction(sequence: str, residues: set[str], window_size: int) -> float:
    if not sequence:
        return 0.0
    effective_window = min(len(sequence), window_size)
    return max(
        _residue_class_fraction(sequence[index : index + effective_window], residues)
        for index in range(0, len(sequence) - effective_window + 1)
    )


def _max_positive_charge_window(sequence: str, window_size: int) -> float:
    if not sequence:
        return 0.0
    effective_window = min(len(sequence), window_size)
    return max(
        sum(max(_charge_class(residue), 0) for residue in sequence[index : index + effective_window])
        / effective_window
        for index in range(0, len(sequence) - effective_window + 1)
    )


def _charge_transition_rate(charge_classes: list[int]) -> float:
    if len(charge_classes) < 2:
        return 0.0
    transitions = sum(
        1
        for current, next_value in zip(charge_classes, charge_classes[1:])
        if current != next_value
    )
    return transitions / (len(charge_classes) - 1)


def _charge_class(residue: str) -> int:
    if residue in {"K", "R"}:
        return 1
    if residue in {"D", "E"}:
        return -1
    return 0


def _count_basic_pairs(sequence: str) -> int:
    return sum(
        1
        for left, right in zip(sequence, sequence[1:])
        if left in {"K", "R"} and right in {"K", "R"}
    )


def _motif_density(sequence: str, motif: str) -> float:
    if not sequence or len(sequence) < len(motif):
        return 0.0
    hits = 0
    for index in range(0, len(sequence) - len(motif) + 1):
        window = sequence[index : index + len(motif)]
        if _matches_motif(window, motif):
            hits += 1
    return hits / max(len(sequence) - len(motif) + 1, 1)


def _matches_motif(window: str, motif: str) -> bool:
    for residue, symbol in zip(window, motif):
        if symbol == "X":
            continue
        if symbol == "B":
            if residue not in {"K", "R"}:
                return False
            continue
        if residue != symbol:
            return False
    return True
