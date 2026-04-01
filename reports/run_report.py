"""Generate markdown summaries, Pareto CSVs, and simple SVG figures."""

from __future__ import annotations

import argparse
import csv
import html
import sys
from pathlib import Path
from typing import Any

if __package__ in {None, ""}:
    project_root = Path(__file__).resolve().parents[1]
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    from core.manifest_loader import load_data_manifest, summarize_manifest
    from models.optimization import pareto_frontier
else:
    from core.manifest_loader import load_data_manifest, summarize_manifest
    from models.optimization import pareto_frontier


PARETO_FIELDS = [
    "rank",
    "candidate_id",
    "sequence",
    "screening_status",
    "composite_screen_score",
    "hs_affinity_reward",
    "hs_selectivity_reward",
    "transport_reward",
    "residence_reward",
    "anticoagulant_penalty",
    "off_target_penalty",
    "heparin_penalty",
    "developability_penalty",
]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate markdown summaries and Pareto artifacts from ranked candidates.",
    )
    parser.add_argument(
        "--ranked",
        default="reports/ranked_candidates.csv",
        help="Input ranked candidate CSV.",
    )
    parser.add_argument(
        "--proposals",
        default="reports/optimization_proposals.csv",
        help="Optional proposal CSV to summarize if present.",
    )
    parser.add_argument(
        "--output-dir",
        default="reports",
        help="Directory to write report artifacts into.",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    ranked_rows = _load_csv_rows(Path(args.ranked))
    proposal_path = Path(args.proposals)
    proposal_rows = _load_csv_rows(proposal_path) if proposal_path.exists() else []
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    score_floor = _quantile(
        [float(row.get("composite_screen_score", 0.0) or 0.0) for row in ranked_rows],
        0.60,
    )
    pareto_input_rows = [
        row for row in ranked_rows if float(row.get("composite_screen_score", 0.0) or 0.0) >= score_floor
    ]
    pareto_rows = pareto_frontier(
        pareto_input_rows,
        maximize_fields=("hs_selectivity_reward", "transport_reward", "residence_reward"),
        minimize_fields=("anticoagulant_penalty", "off_target_penalty", "heparin_penalty"),
        tolerance=0.003,
        max_points=15,
    )
    _write_csv(output_dir / "pareto_frontier.csv", pareto_rows, PARETO_FIELDS)

    _write_svg_scatter(
        output_dir / "pareto_affinity_selectivity.svg",
        ranked_rows,
        x_field="hs_affinity_reward",
        y_field="hs_selectivity_reward",
        title="Affinity vs Selectivity",
    )
    _write_svg_scatter(
        output_dir / "pareto_transport_safety.svg",
        ranked_rows,
        x_field="transport_reward",
        y_field="safety_margin",
        title="Transport vs Safety Margin",
        value_overrides={
            "safety_margin": lambda row: max(
                0.0,
                1.0
                - float(row.get("anticoagulant_penalty", 0.0))
                - (0.7 * float(row.get("off_target_penalty", 0.0)))
                - (0.3 * float(row.get("heparin_penalty", 0.0))),
            )
        },
    )

    manifest_summaries = [
        summarize_manifest(load_data_manifest("data/manifests/hs_gag_panel.json")),
        summarize_manifest(load_data_manifest("data/manifests/off_target_proteins.json")),
        summarize_manifest(load_data_manifest("data/manifests/aggregate_structures.json")),
    ]
    summary_markdown = _build_summary_markdown(
        ranked_rows=ranked_rows,
        pareto_rows=pareto_rows,
        proposal_rows=proposal_rows,
        manifest_summaries=manifest_summaries,
    )
    (output_dir / "screening_summary.md").write_text(summary_markdown, encoding="utf-8")

    print(f"Ranked rows loaded: {len(ranked_rows)}")
    print(f"Pareto rows written: {len(pareto_rows)}")
    print(f"Proposal rows summarized: {len(proposal_rows)}")
    print(f"Artifacts written under: {output_dir}")
    return 0


def _load_csv_rows(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return [_normalize_row(dict(row)) for row in reader]


def _normalize_row(row: dict[str, str]) -> dict[str, Any]:
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
    return normalized


def _parse_scalar(value: str) -> Any:
    try:
        numeric_value = float(value)
    except ValueError:
        return value
    if numeric_value.is_integer():
        return int(numeric_value)
    return numeric_value


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def _quantile(values: list[float], q: float) -> float:
    if not values:
        return 0.0
    sorted_values = sorted(values)
    index = int(round((len(sorted_values) - 1) * q))
    index = max(0, min(index, len(sorted_values) - 1))
    return sorted_values[index]


def _write_svg_scatter(
    path: Path,
    rows: list[dict[str, Any]],
    *,
    x_field: str,
    y_field: str,
    title: str,
    value_overrides: dict[str, Any] | None = None,
) -> None:
    width = 720
    height = 420
    padding = 50
    value_overrides = value_overrides or {}

    def value_for(row: dict[str, Any], field_name: str) -> float:
        override = value_overrides.get(field_name)
        if override is not None:
            return float(override(row))
        return float(row.get(field_name, 0.0) or 0.0)

    x_values = [value_for(row, x_field) for row in rows]
    y_values = [value_for(row, y_field) for row in rows]
    min_x, max_x = min(x_values, default=0.0), max(x_values, default=1.0)
    min_y, max_y = min(y_values, default=0.0), max(y_values, default=1.0)
    if max_x == min_x:
        max_x = min_x + 1.0
    if max_y == min_y:
        max_y = min_y + 1.0

    def scale_x(value: float) -> float:
        return padding + ((value - min_x) / (max_x - min_x)) * (width - (2 * padding))

    def scale_y(value: float) -> float:
        return height - padding - ((value - min_y) / (max_y - min_y)) * (height - (2 * padding))

    status_colors = {
        "pass": "#1d8f6a",
        "review": "#d98f00",
        "reject": "#c53929",
    }

    circles = []
    for row in rows:
        x_value = value_for(row, x_field)
        y_value = value_for(row, y_field)
        fill = status_colors.get(str(row.get("screening_status", "pass")), "#3b82f6")
        label = html.escape(str(row.get("candidate_id", row.get("sequence", ""))))
        circles.append(
            f'<circle cx="{scale_x(x_value):.2f}" cy="{scale_y(y_value):.2f}" r="5" fill="{fill}">'
            f"<title>{label}: {x_field}={x_value:.3f}, {y_field}={y_value:.3f}</title></circle>"
        )

    svg = "\n".join(
        [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
            '<rect width="100%" height="100%" fill="#fffdf8" />',
            f'<text x="{padding}" y="30" font-family="Helvetica, Arial, sans-serif" font-size="20" fill="#1f2937">{html.escape(title)}</text>',
            f'<line x1="{padding}" y1="{height - padding}" x2="{width - padding}" y2="{height - padding}" stroke="#6b7280" stroke-width="2" />',
            f'<line x1="{padding}" y1="{padding}" x2="{padding}" y2="{height - padding}" stroke="#6b7280" stroke-width="2" />',
            f'<text x="{width / 2:.2f}" y="{height - 12}" text-anchor="middle" font-family="Helvetica, Arial, sans-serif" font-size="13" fill="#374151">{html.escape(x_field)}</text>',
            f'<text x="18" y="{height / 2:.2f}" text-anchor="middle" transform="rotate(-90 18 {height / 2:.2f})" font-family="Helvetica, Arial, sans-serif" font-size="13" fill="#374151">{html.escape(y_field)}</text>',
            *circles,
            "</svg>",
        ]
    )
    path.write_text(svg, encoding="utf-8")


def _build_summary_markdown(
    *,
    ranked_rows: list[dict[str, Any]],
    pareto_rows: list[dict[str, Any]],
    proposal_rows: list[dict[str, Any]],
    manifest_summaries: list[dict[str, Any]],
) -> str:
    by_status = {
        status: sum(1 for row in ranked_rows if row.get("screening_status") == status)
        for status in ("pass", "review", "reject")
    }
    top_rows = ranked_rows[:5]
    top_proposals = proposal_rows[:5]

    lines = [
        "# Screening Summary",
        "",
        f"- Screened candidates: {len(ranked_rows)}",
        f"- Pareto-front candidates: {len(pareto_rows)}",
        f"- Passed: {by_status['pass']}",
        f"- Review: {by_status['review']}",
        f"- Reject: {by_status['reject']}",
        "",
        "## Top Ranked Candidates",
        "",
    ]

    for row in top_rows:
        lines.append(
            "- "
            f"{row.get('candidate_id', '')} | score={float(row.get('composite_screen_score', 0.0)):.2f} | "
            f"status={row.get('screening_status', '')} | hs_selectivity={float(row.get('hs_selectivity_reward', 0.0)):.3f} | "
            f"transport={float(row.get('transport_reward', 0.0)):.3f}"
        )

    lines.extend(["", "## Manifest Coverage", ""])
    for summary in manifest_summaries:
        lines.append(
            "- "
            f"{summary['manifest_name']} | records={summary['record_count']} | "
            f"source_dbs={summary['by_source_db']} | statuses={summary['by_status']}"
        )

    if top_proposals:
        lines.extend(["", "## Top Surrogate Proposals", ""])
        for row in top_proposals:
            lines.append(
                "- "
                f"{row.get('candidate_id', '')} | acquisition={float(row.get('acquisition_score', 0.0)):.3f} | "
                f"predicted_score={float(row.get('predicted_composite_score', 0.0)):.2f} | "
                f"novelty={float(row.get('novelty_score', 0.0)):.3f}"
            )

    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
