"""Entry point for generating peptide candidates from a motif CSV."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

if __package__ in {None, ""}:
    project_root = Path(__file__).resolve().parents[1]
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    from peptide.generator import generate_candidates_from_motifs
    from peptide.motifs import load_seed_motifs
else:
    from .generator import generate_candidates_from_motifs
    from .motifs import load_seed_motifs

EXPORT_FIELDS = [
    "candidate_id",
    "sequence",
    "length",
    "approx_net_charge",
    "basic_fraction",
    "hydrophobic_fraction",
    "acidic_fraction",
    "polar_fraction",
    "longest_basic_run",
    "longest_hydrophobic_run",
    "source_motif_ids",
    "generation_method",
    "parent_motif_pattern",
    "passed_filters",
    "filter_flags",
    "warning_flags",
]


def _resolve_seed_motif_path() -> Path:
    package_dir = Path(__file__).resolve().parent
    preferred = package_dir / "seedMotifs.csv"
    if preferred.exists():
        return preferred

    legacy = package_dir.parent / "peptides" / "seedMotifs.csv"
    if legacy.exists():
        return legacy

    raise FileNotFoundError(
        "Could not find seedMotifs.csv in 'peptide/' or legacy 'peptides/' locations."
    )


def _write_candidates_csv(path: Path, candidates: list[object]) -> None:
    rows = [candidate.to_flat_dict() for candidate in candidates]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=EXPORT_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in EXPORT_FIELDS})


def main() -> int:
    input_path = _resolve_seed_motif_path()
    output_dir = Path(__file__).resolve().parent

    motifs = load_seed_motifs(input_path)
    candidates = generate_candidates_from_motifs(motifs)

    passed = [candidate for candidate in candidates if candidate.passed_filters]
    failed = [candidate for candidate in candidates if not candidate.passed_filters]

    _write_candidates_csv(output_dir / "generated_candidates_all.csv", candidates)
    _write_candidates_csv(output_dir / "generated_candidates_passed.csv", passed)
    _write_candidates_csv(output_dir / "generated_candidates_failed.csv", failed)

    print(f"Motifs loaded: {len(motifs)}")
    print(f"Candidates generated: {len(candidates)}")
    print(f"Candidates passed: {len(passed)}")
    print(f"Candidates failed: {len(failed)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
