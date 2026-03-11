"""Loading and normalization utilities for seed motifs."""

from __future__ import annotations

import csv
import re
from pathlib import Path

from .config import ALLOWED_RESIDUES
from .schema import SeedMotif

_CANONICAL_FIELDS: dict[str, set[str]] = {
    "motif_id": {"motif_id", "id", "identifier", "motif_name"},
    "pattern": {"motif", "pattern", "sequence", "raw_pattern", "seed_motif"},
    "source": {"source", "name", "gene", "protein", "origin"},
    "pattern_type": {"pattern_type", "motif_type", "type"},
    "notes": {"note", "notes", "comment", "comments"},
    "natural": {"natural", "is_natural"},
    "synthetic": {"synthetic", "is_synthetic"},
}

_SUPPORTED_PATTERN_CHARS = set(ALLOWED_RESIDUES) | {"B", "X"}
_ALIAS_TO_CANONICAL = {
    alias: canonical_name
    for canonical_name, aliases in _CANONICAL_FIELDS.items()
    for alias in aliases
}


def _normalize_column_name(name: str) -> str:
    cleaned = name.replace("\ufeff", "").strip().lower()
    cleaned = re.sub(r"[^a-z0-9]+", "_", cleaned)
    return cleaned.strip("_")


def _canonicalize_row(row: dict[str, str]) -> dict[str, str]:
    canonical: dict[str, str] = {}
    for key, value in row.items():
        normalized_key = _normalize_column_name(key)
        canonical_name = _ALIAS_TO_CANONICAL.get(normalized_key)
        if canonical_name:
            canonical[canonical_name] = value
    return canonical


def _is_blank_row(row: dict[str, str]) -> bool:
    return not any((value or "").strip() for value in row.values())


def _coerce_optional_text(value: str | None) -> str | None:
    if value is None:
        return None
    cleaned = value.strip()
    return cleaned or None


def _parse_bool(value: str | None) -> bool | None:
    if value is None:
        return None
    cleaned = value.strip().lower()
    if cleaned in {"true", "t", "yes", "y", "1", "natural"}:
        return True
    if cleaned in {"false", "f", "no", "n", "0", "synthetic"}:
        return False
    return None


def _build_motif_id(source: str | None, row_index: int) -> str:
    if source:
        stem = re.sub(r"[^A-Za-z0-9]+", "_", source.strip()).strip("_").lower()
        if stem:
            return f"{stem}_{row_index - 1:04d}"
    return f"motif_{row_index - 1:04d}"


def normalize_motif_pattern(pattern: str) -> str:
    normalized = "".join(pattern.split()).upper()
    if not normalized:
        raise ValueError("Encountered an empty motif pattern after normalization.")

    invalid_chars = sorted(set(normalized) - _SUPPORTED_PATTERN_CHARS)
    if invalid_chars:
        raise ValueError(
            f"Unsupported motif characters in pattern '{pattern}': {', '.join(invalid_chars)}"
        )
    return normalized


def infer_pattern_type(normalized_pattern: str, pattern_type: str | None = None) -> str:
    if "B" in normalized_pattern or "X" in normalized_pattern:
        return "pattern"
    if pattern_type:
        lowered = pattern_type.strip().lower()
        if lowered in {"exact", "sequence"}:
            return "exact"
        if lowered in {"pattern", "degenerate", "template"}:
            return "pattern"
    return "exact"


def validate_supported_pattern(normalized_pattern: str, pattern_type: str) -> None:
    if pattern_type == "exact" and ("B" in normalized_pattern or "X" in normalized_pattern):
        raise ValueError(
            f"Exact motif '{normalized_pattern}' contains pattern-only symbols B or X."
        )


def _infer_is_natural(row: dict[str, str]) -> bool | None:
    natural_value = _parse_bool(row.get("natural"))
    synthetic_value = _parse_bool(row.get("synthetic"))
    if natural_value is not None:
        return natural_value
    if synthetic_value is not None:
        return not synthetic_value
    source = (row.get("source") or "").strip().lower()
    if "synthetic" in source:
        return False
    return None


def load_seed_motifs(csv_path: str | Path) -> list[SeedMotif]:
    path = Path(csv_path)
    if not path.exists():
        raise FileNotFoundError(f"Seed motif CSV not found: {path}")

    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"Seed motif CSV has no header row: {path}")

        normalized_headers = {_normalize_column_name(name) for name in reader.fieldnames}
        pattern_headers = _CANONICAL_FIELDS["pattern"]
        if not normalized_headers.intersection(pattern_headers):
            raise ValueError(
                "Seed motif CSV must include a motif column such as "
                "'motif', 'pattern', or 'sequence'."
            )

        motifs: list[SeedMotif] = []
        seen: set[tuple[object, ...]] = set()

        for row_index, raw_row in enumerate(reader, start=2):
            if _is_blank_row(raw_row):
                continue

            row = _canonicalize_row(raw_row)
            raw_pattern = _coerce_optional_text(row.get("pattern"))
            if not raw_pattern:
                continue

            normalized_pattern = normalize_motif_pattern(raw_pattern)
            pattern_type = infer_pattern_type(normalized_pattern, row.get("pattern_type"))
            validate_supported_pattern(normalized_pattern, pattern_type)

            source = _coerce_optional_text(row.get("source"))
            motif_id = _coerce_optional_text(row.get("motif_id")) or _build_motif_id(
                source, row_index
            )
            notes = _coerce_optional_text(row.get("notes"))
            is_natural = _infer_is_natural(row)

            dedupe_key = (normalized_pattern, pattern_type, source, is_natural, notes)
            if dedupe_key in seen:
                continue
            seen.add(dedupe_key)

            motifs.append(
                SeedMotif(
                    motif_id=motif_id,
                    raw_pattern=raw_pattern,
                    normalized_pattern=normalized_pattern,
                    pattern_type=pattern_type,
                    source=source,
                    is_natural=is_natural,
                    notes=notes,
                )
            )

        if not motifs:
            raise ValueError(f"No usable motif rows were loaded from {path}.")

        return motifs
