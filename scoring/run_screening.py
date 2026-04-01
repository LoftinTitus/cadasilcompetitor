"""CLI entry point for screening peptide candidates and exporting ranked CSVs."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

if __package__ in {None, ""}:
    project_root = Path(__file__).resolve().parents[1]
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    from scoring.screening import screen_candidates_from_csv
else:
    from .screening import screen_candidates_from_csv


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run intro screening and export a ranked candidate CSV.",
    )
    parser.add_argument(
        "--input",
        default="peptide/generated_candidates_passed.csv",
        help="Path to the candidate CSV to score.",
    )
    parser.add_argument(
        "--output",
        default="reports/ranked_candidates.csv",
        help="Path for the ranked CSV output.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional limit on how many candidate rows to score.",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    ranked_rows = screen_candidates_from_csv(
        Path(args.input),
        output_csv_path=Path(args.output),
        limit=args.limit,
    )
    print(f"Candidates scored: {len(ranked_rows)}")
    print(f"Ranked CSV written to: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
