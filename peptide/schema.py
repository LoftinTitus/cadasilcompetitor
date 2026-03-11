"""Dataclasses used across the peptide package."""

from __future__ import annotations

from dataclasses import dataclass, field


def _serialize_list(values: list[object]) -> str:
    return ";".join(str(value) for value in values)


def _serialize_mapping(values: dict[str, object]) -> str:
    return ";".join(f"{key}:{values[key]}" for key in sorted(values))


@dataclass(slots=True)
class SeedMotif:
    motif_id: str
    raw_pattern: str
    normalized_pattern: str
    pattern_type: str
    source: str | None = None
    is_natural: bool | None = None
    notes: str | None = None

    def __post_init__(self) -> None:
        self.motif_id = self.motif_id.strip()
        self.raw_pattern = self.raw_pattern.strip().upper()
        self.normalized_pattern = self.normalized_pattern.strip().upper()
        self.pattern_type = self.pattern_type.strip().lower()
        self.source = self.source.strip() if self.source else None
        self.notes = self.notes.strip() if self.notes else None


@dataclass(slots=True)
class PeptideCandidate:
    candidate_id: str
    sequence: str
    length: int
    source_motif_ids: list[str] = field(default_factory=list)
    generation_method: str = ""
    parent_motif_pattern: str | None = None
    approx_net_charge: int | float = 0
    residue_counts: dict[str, int] = field(default_factory=dict)
    residue_fractions: dict[str, float] = field(default_factory=dict)
    basic_fraction: float = 0.0
    hydrophobic_fraction: float = 0.0
    polar_fraction: float = 0.0
    acidic_fraction: float = 0.0
    longest_basic_run: int = 0
    longest_hydrophobic_run: int = 0
    basic_positions: list[int] = field(default_factory=list)
    charge_spacing: list[int] = field(default_factory=list)
    motif_hits: list[str] = field(default_factory=list)
    passed_filters: bool = False
    filter_flags: list[str] = field(default_factory=list)
    warning_flags: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.sequence = self.sequence.strip().upper()
        self.length = len(self.sequence)

    def to_flat_dict(self) -> dict[str, object]:
        return {
            "candidate_id": self.candidate_id,
            "sequence": self.sequence,
            "length": self.length,
            "approx_net_charge": self.approx_net_charge,
            "basic_fraction": self.basic_fraction,
            "hydrophobic_fraction": self.hydrophobic_fraction,
            "acidic_fraction": self.acidic_fraction,
            "polar_fraction": self.polar_fraction,
            "longest_basic_run": self.longest_basic_run,
            "longest_hydrophobic_run": self.longest_hydrophobic_run,
            "source_motif_ids": _serialize_list(self.source_motif_ids),
            "generation_method": self.generation_method,
            "parent_motif_pattern": self.parent_motif_pattern or "",
            "passed_filters": self.passed_filters,
            "filter_flags": _serialize_list(self.filter_flags),
            "warning_flags": _serialize_list(self.warning_flags),
            "basic_positions": _serialize_list(self.basic_positions),
            "charge_spacing": _serialize_list(self.charge_spacing),
            "motif_hits": _serialize_list(self.motif_hits),
            "residue_counts": _serialize_mapping(self.residue_counts),
            "residue_fractions": _serialize_mapping(
                {key: f"{value:.4f}" for key, value in self.residue_fractions.items()}
            ),
        }
