"""CLI for fitting a lightweight surrogate and proposing next candidates."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Any

if __package__ in {None, ""}:
    project_root = Path(__file__).resolve().parents[1]
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    from models.optimization import DEFAULT_TARGET_FIELDS, propose_candidates
    from models.screening_calibrator import try_load_screening_calibrator
    from models.surrogate import fit_bootstrap_surrogate
    from scoring.screening import load_candidates_from_csv
else:
    from .optimization import DEFAULT_TARGET_FIELDS, propose_candidates
    from .screening_calibrator import try_load_screening_calibrator
    from .surrogate import fit_bootstrap_surrogate
    from scoring.screening import load_candidates_from_csv


PROPOSAL_FIELDS = [
    "proposal_rank",
    "candidate_id",
    "sequence",
    "proposal_source",
    "length",
    "approx_net_charge",
    "warning_flags",
    "filter_flags",
    "acquisition_score",
    "novelty_score",
    "uncertainty_score",
    "ml_screening_score",
    "ml_predicted_composite_score",
    "ml_screening_uncertainty",
    "predicted_composite_score",
    "pred_hs_affinity_reward",
    "pred_hs_selectivity_reward",
    "pred_transport_reward",
    "pred_residence_reward",
    "pred_anticoagulant_penalty",
    "pred_off_target_penalty",
    "pred_heparin_penalty",
    "pred_developability_penalty",
    "std_hs_selectivity_reward",
    "std_transport_reward",
    "std_anticoagulant_penalty",
]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Fit a lightweight surrogate on screened rows and propose next candidates.",
    )
    parser.add_argument(
        "--training",
        default="reports/ranked_candidates.csv",
        help="CSV of already screened candidates.",
    )
    parser.add_argument(
        "--candidate-pool",
        default="peptide/generated_candidates_passed.csv",
        help="CSV of candidate sequences to rank with the surrogate.",
    )
    parser.add_argument(
        "--output",
        default="reports/optimization_proposals.csv",
        help="Output CSV for ranked proposals.",
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=50,
        help="Number of proposals to write.",
    )
    parser.add_argument(
        "--screening-calibrator",
        default="reports/screening_calibrator.json",
        help="Optional learned screening calibrator JSON.",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    training_rows = _load_screened_rows(Path(args.training))
    candidate_pool_rows = load_candidates_from_csv(Path(args.candidate_pool))
    surrogate = fit_bootstrap_surrogate(
        training_rows,
        target_fields=DEFAULT_TARGET_FIELDS,
    )
    screening_calibrator = try_load_screening_calibrator(Path(args.screening_calibrator))
    proposals = propose_candidates(
        candidate_pool_rows,
        trained_surrogate=surrogate,
        reference_rows=training_rows,
        screening_calibrator=screening_calibrator,
        top_k=args.top_k,
    )
    if not proposals:
        proposals = propose_candidates(
            candidate_pool_rows,
            trained_surrogate=surrogate,
            reference_rows=training_rows,
            screening_calibrator=screening_calibrator,
            top_k=args.top_k,
            allow_seen=True,
        )
    _write_csv(Path(args.output), proposals, PROPOSAL_FIELDS)
    print(f"Training rows loaded: {len(training_rows)}")
    print(f"Candidate pool rows loaded: {len(candidate_pool_rows)}")
    print(f"Screening calibrator loaded: {screening_calibrator is not None}")
    print(f"Proposals written: {len(proposals)}")
    print(f"Output CSV written to: {args.output}")
    return 0


def _load_screened_rows(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = [dict(row) for row in reader]

    normalized_rows: list[dict[str, Any]] = []
    for row in rows:
        normalized: dict[str, Any] = {}
        for key, value in row.items():
            if value is None or value == "":
                normalized[key] = value
                continue
            if key in {"sequence", "candidate_id", "screening_status", "tier_a_best_variant_id", "transport_probe_variant_id"}:
                normalized[key] = value
                continue
            if key in {"warning_flags", "filter_flags", "risk_flags"}:
                normalized[key] = [item for item in value.split(";") if item]
                continue
            normalized[key] = _parse_scalar(value)
        normalized_rows.append(normalized)
    return normalized_rows


def _parse_scalar(value: str) -> Any:
    try:
        numeric_value = float(value)
    except ValueError:
        return value
    if numeric_value.is_integer():
        return int(numeric_value)
    return numeric_value


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


if __name__ == "__main__":
    raise SystemExit(main())
