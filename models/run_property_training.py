"""CLI for training a lightweight ML surrogate for candidate properties."""

from __future__ import annotations

import argparse

from .property_surrogate import (
    DEFAULT_PROPERTY_ESTIMATOR_PATH,
    fit_property_surrogate,
    load_property_training_rows,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Train a lightweight property surrogate from ranked candidate rows.",
    )
    parser.add_argument(
        "--training",
        default="reports/ranked_candidates.csv",
        help="CSV containing candidate features and estimated property targets.",
    )
    parser.add_argument(
        "--output",
        default=str(DEFAULT_PROPERTY_ESTIMATOR_PATH),
        help="JSON path for the trained property surrogate.",
    )
    parser.add_argument(
        "--ensemble-size",
        type=int,
        default=25,
        help="Number of bootstrap models to fit per property target.",
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

    rows = load_property_training_rows(args.training)
    model = fit_property_surrogate(
        rows,
        ensemble_size=args.ensemble_size,
        ridge_alpha=args.ridge_alpha,
        random_seed=args.random_seed,
        training_source=args.training,
    )
    model.write_json(args.output)

    print(f"Training rows loaded: {len(rows)}")
    print(f"Property surrogate written to: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
