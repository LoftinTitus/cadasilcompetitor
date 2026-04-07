"""CLI for training a learned screening calibration model."""

from __future__ import annotations

import argparse

from .screening_calibrator import (
    DEFAULT_SCREENING_CALIBRATOR_PATH,
    fit_screening_calibrator,
    load_screening_training_rows,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Train a learned screening-status/composite-score calibrator.",
    )
    parser.add_argument(
        "--training",
        default="reports/ranked_candidates.csv",
        help="CSV of screened candidate rows.",
    )
    parser.add_argument(
        "--output",
        default=str(DEFAULT_SCREENING_CALIBRATOR_PATH),
        help="Output JSON path for the trained calibrator.",
    )
    parser.add_argument(
        "--ensemble-size",
        type=int,
        default=25,
        help="Number of bootstrap models per calibration target.",
    )
    parser.add_argument(
        "--ridge-alpha",
        type=float,
        default=0.25,
        help="Ridge regularization strength.",
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=7,
        help="Random seed for bootstrap resampling.",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    rows = load_screening_training_rows(args.training)
    calibrator = fit_screening_calibrator(
        rows,
        ensemble_size=args.ensemble_size,
        ridge_alpha=args.ridge_alpha,
        random_seed=args.random_seed,
        training_source=args.training,
    )
    calibrator.write_json(args.output)

    print(f"Training rows loaded: {len(rows)}")
    print(f"Screening calibrator written to: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
