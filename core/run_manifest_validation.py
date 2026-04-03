"""CLI to summarize manifest completeness and panel/manifest consistency."""

from __future__ import annotations

import argparse
from pathlib import Path

from core.manifest_loader import load_data_manifest
from core.manifest_registry import build_manifest_validation_summary
from core.config_loader import load_hs_variant_panel


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate manifest completeness and consistency with the configured HS panel.",
    )
    parser.add_argument(
        "--panel",
        default="config/hs_variant_panel.yaml",
        help="Path to the HS panel YAML.",
    )
    parser.add_argument(
        "--hs-manifest",
        default="data/manifests/hs_gag_panel.json",
        help="Path to the HS/GAG manifest JSON.",
    )
    parser.add_argument(
        "--off-target-manifest",
        default="data/manifests/off_target_proteins.json",
        help="Path to the off-target manifest JSON.",
    )
    parser.add_argument(
        "--aggregate-manifest",
        default="data/manifests/aggregate_structures.json",
        help="Path to the aggregate structure manifest JSON.",
    )
    parser.add_argument(
        "--output",
        default="reports/manifest_validation_summary.md",
        help="Markdown output path.",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    summary = build_manifest_validation_summary(
        hs_panel=load_hs_variant_panel(args.panel),
        hs_manifest=load_data_manifest(args.hs_manifest),
        off_target_manifest=load_data_manifest(args.off_target_manifest),
        aggregate_manifest=load_data_manifest(args.aggregate_manifest),
    )

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(summary, encoding="utf-8")

    print(f"Manifest validation summary written to: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
