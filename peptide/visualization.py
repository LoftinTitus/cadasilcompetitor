"""Generate report-ready figures for the peptide pipeline."""

from __future__ import annotations

import argparse
import os
import re
import tempfile
from pathlib import Path
from typing import Iterable

_CACHE_ROOT = Path(tempfile.gettempdir()) / "cadasilcompetitor-plot-cache"
_CACHE_ROOT.mkdir(parents=True, exist_ok=True)
(_CACHE_ROOT / "matplotlib").mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_CACHE_ROOT / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(_CACHE_ROOT))

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns
except ModuleNotFoundError as exc:  # pragma: no cover - import guard
    raise SystemExit(
        f"Missing dependency '{exc.name}'. Activate the project virtualenv before running "
        "peptide/visualization.py."
    ) from exc


STATUS_PALETTE = {
    "pass": "#157f6b",
    "review": "#e09f3e",
    "reject": "#d1495b",
}


def build_parser() -> argparse.ArgumentParser:
    project_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description="Generate peptide pipeline figures from existing CSV artifacts.",
    )
    parser.add_argument(
        "--motifs",
        default=str(project_root / "peptide" / "seedMotifs.csv"),
        help="Path to the seed motif CSV.",
    )
    parser.add_argument(
        "--generated",
        default=str(project_root / "peptide" / "generated_candidates_all.csv"),
        help="Path to the generated candidate CSV.",
    )
    parser.add_argument(
        "--passed",
        default=str(project_root / "peptide" / "generated_candidates_passed.csv"),
        help="Path to the passed candidate CSV.",
    )
    parser.add_argument(
        "--ranked",
        default=str(project_root / "reports" / "ranked_candidates.csv"),
        help="Path to the ranked candidate CSV.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(project_root / "reports" / "figures"),
        help="Directory where figure files should be written.",
    )
    return parser


def _normalize_column_name(name: str) -> str:
    cleaned = name.replace("\ufeff", "").strip().lower()
    cleaned = re.sub(r"[^a-z0-9]+", "_", cleaned)
    return cleaned.strip("_")


def _pretty_label(name: str) -> str:
    return _normalize_column_name(name).replace("_", " ").title()


def _resolve_column(df: pd.DataFrame, *aliases: str) -> str:
    column_map = {_normalize_column_name(column): column for column in df.columns}
    for alias in aliases:
        resolved = column_map.get(_normalize_column_name(alias))
        if resolved is not None:
            return resolved
    raise KeyError(f"Could not find any of these columns: {', '.join(aliases)}")


def _optional_column(df: pd.DataFrame, *aliases: str) -> str | None:
    try:
        return _resolve_column(df, *aliases)
    except KeyError:
        return None


def _finalize_figure(fig: plt.Figure, save_path: str | Path | None) -> plt.Figure | Path:
    fig.tight_layout()
    if save_path is None:
        return fig

    output_path = Path(save_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output_path


def _set_plot_style() -> None:
    sns.set_theme(style="whitegrid", context="talk")
    plt.rcParams.update(
        {
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titleweight": "bold",
            "figure.facecolor": "#fcfcfc",
            "axes.facecolor": "#fcfcfc",
            "savefig.facecolor": "#fcfcfc",
        }
    )


def _load_csv(path: str | Path) -> pd.DataFrame:
    csv_path = Path(path)
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV not found: {csv_path}")
    return pd.read_csv(csv_path)


def _clean_string_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().replace({"": pd.NA, "nan": pd.NA, "None": pd.NA})


def _prepare_ranked_frame(metadata_df: pd.DataFrame) -> pd.DataFrame:
    ranked = metadata_df.copy()
    for column_name in (
        "hs_affinity_reward",
        "hs_selectivity_reward",
        "transport_reward",
        "anticoagulant_penalty",
        "off_target_penalty",
        "heparin_penalty",
        "composite_screen_score",
        "length",
    ):
        resolved = _optional_column(ranked, column_name)
        if resolved is not None:
            ranked[resolved] = pd.to_numeric(ranked[resolved], errors="coerce")

    anticoagulant = _optional_column(ranked, "anticoagulant_penalty")
    off_target = _optional_column(ranked, "off_target_penalty")
    heparin = _optional_column(ranked, "heparin_penalty")
    if anticoagulant and off_target and heparin:
        ranked["safety_margin"] = (
            1.0
            - ranked[anticoagulant].fillna(0.0)
            - (0.7 * ranked[off_target].fillna(0.0))
            - (0.3 * ranked[heparin].fillna(0.0))
        ).clip(lower=0.0)
    return ranked


def plot_pipeline_flowchart(save_path: str | Path | None = None) -> plt.Figure | Path:
    import matplotlib.patches as mpatches

    fig, ax = plt.subplots(figsize=(12, 3.6))
    steps = [
        ("Load Seed Motifs", "#d8f3dc"),
        ("Generate Candidates", "#b7e4c7"),
        ("Apply Filters", "#95d5b2"),
        ("Score And Rank", "#74c69d"),
        ("Write Reports", "#52b788"),
    ]

    for index, (label, color) in enumerate(steps):
        x_position = index * 2.3
        box = mpatches.FancyBboxPatch(
            (x_position, 0),
            1.9,
            1.0,
            boxstyle="round,pad=0.18,rounding_size=0.08",
            fc=color,
            ec="#1f2937",
            linewidth=1.2,
        )
        ax.add_patch(box)
        ax.text(
            x_position + 0.95,
            0.5,
            label,
            ha="center",
            va="center",
            fontsize=12,
            color="#102a43",
        )
        if index < len(steps) - 1:
            ax.annotate(
                "",
                xy=(x_position + 2.15, 0.5),
                xytext=(x_position + 1.9, 0.5),
                arrowprops={"arrowstyle": "->", "lw": 2, "color": "#334e68"},
            )

    ax.set_xlim(-0.3, (len(steps) * 2.3) - 0.1)
    ax.set_ylim(-0.3, 1.3)
    ax.axis("off")
    ax.set_title("Peptide Pipeline Overview", loc="left")
    return _finalize_figure(fig, save_path)


def plot_motif_distribution(
    motifs_df: pd.DataFrame,
    save_path: str | Path | None = None,
    *,
    top_n: int = 20,
) -> plt.Figure | Path:
    motif_col = _resolve_column(
        motifs_df,
        "motif",
        "pattern",
        "sequence",
        "raw_pattern",
        "seed_motif",
        "normalized_pattern",
    )
    counts = _clean_string_series(motifs_df[motif_col]).dropna().value_counts().head(top_n)

    fig_height = max(5.5, 0.35 * max(len(counts), 1) + 2.5)
    fig, ax = plt.subplots(figsize=(12, fig_height))
    bars = ax.barh(counts.index[::-1], counts.values[::-1], color="#2a9d8f")
    ax.bar_label(bars, padding=4, fontsize=10)
    ax.set_title("Most Frequent Seed Motifs", loc="left")
    ax.set_xlabel("Occurrences")
    ax.set_ylabel("Motif")
    return _finalize_figure(fig, save_path)


def plot_filtering_funnel(
    counts: Iterable[int],
    labels: Iterable[str],
    save_path: str | Path | None = None,
) -> plt.Figure | Path:
    stage_labels = list(labels)
    stage_counts = [int(count) for count in counts]

    fig, ax = plt.subplots(figsize=(10, 5))
    colors = ["#577590", "#43aa8b", "#277da1"][: len(stage_counts)]
    bars = ax.bar(stage_labels, stage_counts, color=colors)
    ax.bar_label(bars, padding=4, fontsize=11)

    ax.set_title("Candidate Funnel", loc="left")
    ax.set_ylabel("Sequences")
    ax.set_xlabel("Stage")
    max_count = max(stage_counts, default=0)
    ax.set_ylim(0, max_count * 1.18 if max_count else 1)
    for index in range(1, len(stage_counts)):
        previous = stage_counts[index - 1]
        current = stage_counts[index]
        retention = 0.0 if previous == 0 else (current / previous) * 100
        ax.text(
            index,
            current + max_count * 0.04,
            f"{retention:.1f}%",
            ha="center",
            va="bottom",
            fontsize=10,
            color="#334e68",
        )
    return _finalize_figure(fig, save_path)


def plot_screening_status_breakdown(
    metadata_df: pd.DataFrame,
    save_path: str | Path | None = None,
) -> plt.Figure | Path:
    status_col = _resolve_column(metadata_df, "screening_status", "status")
    status_counts = (
        _clean_string_series(metadata_df[status_col]).dropna().str.lower().value_counts()
    )
    ordered_statuses = [status for status in ("pass", "review", "reject") if status in status_counts]

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(
        ordered_statuses,
        [int(status_counts[status]) for status in ordered_statuses],
        color=[STATUS_PALETTE[status] for status in ordered_statuses],
    )
    ax.bar_label(bars, padding=4, fontsize=11)
    ax.set_title("Screening Status Breakdown", loc="left")
    ax.set_xlabel("Status")
    ax.set_ylabel("Candidates")
    return _finalize_figure(fig, save_path)


def plot_property_distribution(
    metadata_df: pd.DataFrame,
    property_name: str,
    by: str | None = None,
    save_path: str | Path | None = None,
) -> plt.Figure | Path:
    value_col = _resolve_column(metadata_df, property_name)
    plot_df = metadata_df.copy()
    plot_df[value_col] = pd.to_numeric(plot_df[value_col], errors="coerce")
    plot_df = plot_df.dropna(subset=[value_col])

    fig, ax = plt.subplots(figsize=(10, 5.5))
    if by:
        category_col = _resolve_column(plot_df, by)
        plot_df[category_col] = _clean_string_series(plot_df[category_col]).fillna("unknown")
        sns.violinplot(
            data=plot_df,
            x=category_col,
            y=value_col,
            hue=category_col,
            palette=STATUS_PALETTE,
            inner="quartile",
            cut=0,
            ax=ax,
        )
        legend = ax.get_legend()
        if legend is not None:
            legend.remove()
        ax.set_xlabel(_pretty_label(category_col))
    else:
        sns.histplot(plot_df[value_col], kde=True, bins=30, color="#277da1", ax=ax)
        ax.set_xlabel(_pretty_label(value_col))

    ax.set_title(f"{_pretty_label(value_col)} Distribution", loc="left")
    ax.set_ylabel(_pretty_label(value_col) if by else "Count")
    return _finalize_figure(fig, save_path)


def plot_scatter_properties(
    metadata_df: pd.DataFrame,
    x: str,
    y: str,
    hue: str | None = None,
    save_path: str | Path | None = None,
) -> plt.Figure | Path:
    x_col = _resolve_column(metadata_df, x)
    y_col = _resolve_column(metadata_df, y)
    plot_df = metadata_df.copy()
    plot_df[x_col] = pd.to_numeric(plot_df[x_col], errors="coerce")
    plot_df[y_col] = pd.to_numeric(plot_df[y_col], errors="coerce")
    plot_df = plot_df.dropna(subset=[x_col, y_col])

    hue_col = _resolve_column(plot_df, hue) if hue else None
    if hue_col:
        plot_df[hue_col] = _clean_string_series(plot_df[hue_col]).fillna("unknown").str.lower()

    fig, ax = plt.subplots(figsize=(9.5, 6.5))
    sns.scatterplot(
        data=plot_df,
        x=x_col,
        y=y_col,
        hue=hue_col,
        palette=STATUS_PALETTE if hue_col else None,
        s=85,
        linewidth=0.5,
        edgecolor="white",
        ax=ax,
    )

    candidate_col = _optional_column(plot_df, "candidate_id", "candidate")
    score_col = _optional_column(plot_df, "composite_screen_score", "score")
    if candidate_col and score_col:
        plot_df[score_col] = pd.to_numeric(plot_df[score_col], errors="coerce")
        annotations = plot_df.dropna(subset=[score_col]).nlargest(8, score_col)
        for _, row in annotations.iterrows():
            ax.text(
                row[x_col],
                row[y_col],
                str(row[candidate_col]),
                fontsize=8,
                color="#1f2937",
                alpha=0.85,
            )

    ax.set_title(f"{_pretty_label(y_col)} vs {_pretty_label(x_col)}", loc="left")
    ax.set_xlabel(_pretty_label(x_col))
    ax.set_ylabel(_pretty_label(y_col))
    return _finalize_figure(fig, save_path)


def plot_motif_heatmap(
    candidates_df: pd.DataFrame,
    motif_col: str = "parent_motif_pattern",
    save_path: str | Path | None = None,
    *,
    group_col: str = "generation_method",
    top_n: int = 18,
) -> plt.Figure | Path:
    resolved_motif_col = _resolve_column(
        candidates_df,
        motif_col,
        "parent_motif_pattern",
        "motif",
        "pattern",
        "motif_hits",
    )
    resolved_group_col = _resolve_column(
        candidates_df,
        group_col,
        "generation_method",
        "screening_status",
        "passed_filters",
    )

    plot_df = candidates_df[[resolved_motif_col, resolved_group_col]].copy()
    plot_df[resolved_motif_col] = _clean_string_series(plot_df[resolved_motif_col])
    plot_df[resolved_group_col] = _clean_string_series(plot_df[resolved_group_col])
    plot_df = plot_df.dropna()

    top_motifs = plot_df[resolved_motif_col].value_counts().head(top_n).index
    heatmap_data = pd.crosstab(
        plot_df.loc[plot_df[resolved_motif_col].isin(top_motifs), resolved_motif_col],
        plot_df.loc[plot_df[resolved_motif_col].isin(top_motifs), resolved_group_col],
    )
    heatmap_data = heatmap_data.loc[heatmap_data.sum(axis=1).sort_values(ascending=False).index]

    fig_height = max(6.5, 0.38 * max(len(heatmap_data), 1) + 2.5)
    fig, ax = plt.subplots(figsize=(11.5, fig_height))
    sns.heatmap(
        heatmap_data,
        cmap="YlGnBu",
        linewidths=0.5,
        annot=True,
        fmt="d",
        cbar_kws={"label": "Candidates"},
        ax=ax,
    )
    ax.set_title("Top Motifs By Generation Method", loc="left")
    ax.set_xlabel(_pretty_label(resolved_group_col))
    ax.set_ylabel(_pretty_label(resolved_motif_col))
    return _finalize_figure(fig, save_path)


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    _set_plot_style()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    motifs_df = _load_csv(args.motifs)
    generated_df = _load_csv(args.generated)
    passed_df = _load_csv(args.passed)
    ranked_df = _prepare_ranked_frame(_load_csv(args.ranked))

    status_col = _resolve_column(ranked_df, "screening_status", "status")
    passing_status_count = int(
        _clean_string_series(ranked_df[status_col]).fillna("").str.lower().eq("pass").sum()
    )

    figure_paths = [
        plot_pipeline_flowchart(output_dir / "pipeline_overview.png"),
        plot_motif_distribution(motifs_df, output_dir / "seed_motif_frequency.png"),
        plot_filtering_funnel(
            [len(generated_df), len(passed_df), passing_status_count],
            ["Generated", "Passed peptide filters", "Screening pass"],
            output_dir / "candidate_funnel.png",
        ),
        plot_screening_status_breakdown(ranked_df, output_dir / "screening_status_breakdown.png"),
        plot_property_distribution(
            ranked_df,
            "length",
            by="screening_status",
            save_path=output_dir / "candidate_length_by_status.png",
        ),
        plot_scatter_properties(
            ranked_df,
            "hs_affinity_reward",
            "hs_selectivity_reward",
            hue="screening_status",
            save_path=output_dir / "affinity_selectivity_tradeoff.png",
        ),
        plot_scatter_properties(
            ranked_df,
            "transport_reward",
            "safety_margin",
            hue="screening_status",
            save_path=output_dir / "transport_safety_tradeoff.png",
        ),
        plot_motif_heatmap(generated_df, save_path=output_dir / "motif_generation_heatmap.png"),
    ]

    print(f"Figures written under: {output_dir}")
    for figure_path in figure_paths:
        print(f"- {figure_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
