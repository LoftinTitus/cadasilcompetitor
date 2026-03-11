"""Simple configuration constants for peptide generation."""

from __future__ import annotations

from typing import Final

PEPTIDE_MIN_LENGTH: Final[int] = 8
PEPTIDE_MAX_LENGTH: Final[int] = 20

ALLOWED_RESIDUES: Final[tuple[str, ...]] = tuple("ACDEFGHIKLMNPQRSTVWY")
BASIC_RESIDUES: Final[frozenset[str]] = frozenset({"K", "R"})
ACIDIC_RESIDUES: Final[frozenset[str]] = frozenset({"D", "E"})
HYDROPHOBIC_RESIDUES: Final[frozenset[str]] = frozenset({"A", "F", "I", "L", "M", "V", "W", "Y"})
POLAR_RESIDUES: Final[frozenset[str]] = frozenset({"C", "G", "H", "N", "P", "Q", "S", "T"})

ALLOWED_FLANKING_RESIDUES: Final[tuple[str, ...]] = ("G", "Q", "S", "A")
ALLOWED_SPACER_RESIDUES: Final[tuple[str, ...]] = ("G", "Q", "S", "A", "N")

MAX_VARIANTS_PER_MOTIF: Final[int] = 24
MAX_TOTAL_CANDIDATES: Final[int] = 500

MIN_ALLOWED_NET_CHARGE: Final[int] = 1
MAX_ALLOWED_NET_CHARGE: Final[int] = 12
MAX_BASIC_FRACTION: Final[float] = 0.70
MAX_HYDROPHOBIC_FRACTION: Final[float] = 0.55
MAX_CONSECUTIVE_BASIC_RESIDUES: Final[int] = 4
MAX_CONSECUTIVE_HYDROPHOBIC_RESIDUES: Final[int] = 4
MAX_REPEATED_IDENTICAL_RESIDUE_RUN: Final[int] = 3

RANDOM_SEED: Final[int] = 7
