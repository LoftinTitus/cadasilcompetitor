"""Helpers for linking configured panels to JSON-backed data manifests."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from core.config_loader import HSVariantPanelConfig, load_hs_variant_panel
from core.manifest_loader import DataManifest, ManifestRecord, load_data_manifest, summarize_manifest

DEFAULT_HS_PANEL_PATH = Path("config/hs_variant_panel.yaml")
DEFAULT_HS_MANIFEST_PATH = Path("data/manifests/hs_gag_panel.json")
DEFAULT_OFF_TARGET_MANIFEST_PATH = Path("data/manifests/off_target_proteins.json")
DEFAULT_AGGREGATE_MANIFEST_PATH = Path("data/manifests/aggregate_structures.json")


@dataclass(frozen=True, slots=True)
class HSVariantManifestAnnotation:
    variant_id: str
    display_name: str
    manifest_status: str
    manifest_source_accession: str | None
    manifest_notes: str | None
    atomistic_preparation_status: str
    glytoucan_accession: str | None
    provenance_source_status: str


def load_default_manifests() -> dict[str, DataManifest]:
    return {
        "hs_gag_panel": load_data_manifest(DEFAULT_HS_MANIFEST_PATH),
        "off_target_proteins": load_data_manifest(DEFAULT_OFF_TARGET_MANIFEST_PATH),
        "aggregate_structures": load_data_manifest(DEFAULT_AGGREGATE_MANIFEST_PATH),
    }


def index_manifest_records(manifest: DataManifest) -> dict[str, ManifestRecord]:
    return {record.record_id: record for record in manifest.records}


def annotate_hs_panel_variants(
    panel: HSVariantPanelConfig,
    manifest: DataManifest,
) -> dict[str, HSVariantManifestAnnotation]:
    records_by_id = index_manifest_records(manifest)
    annotations: dict[str, HSVariantManifestAnnotation] = {}
    for variant in panel.variants:
        record = records_by_id.get(variant.variant_id)
        structural_annotations = variant.structural_annotations
        provenance = variant.provenance
        annotations[variant.variant_id] = HSVariantManifestAnnotation(
            variant_id=variant.variant_id,
            display_name=variant.display_name,
            manifest_status=record.status if record is not None else "missing_manifest_record",
            manifest_source_accession=(
                record.source_accession
                if record is not None
                else _normalize_optional_text(structural_annotations.get("glytoucan_accession"))
            ),
            manifest_notes=record.notes if record is not None else None,
            atomistic_preparation_status=str(
                structural_annotations.get("atomistic_preparation_status", "unknown")
            ),
            glytoucan_accession=_normalize_optional_text(
                structural_annotations.get("glytoucan_accession")
            )
            or (record.source_accession if record is not None else None),
            provenance_source_status=str(
                provenance.get("source_status", panel.source_status or "unknown")
            ),
        )
    return annotations


def summarize_hs_panel_readiness(
    panel: HSVariantPanelConfig,
    manifest: DataManifest,
) -> dict[str, Any]:
    annotations = annotate_hs_panel_variants(panel, manifest)
    manifest_records = index_manifest_records(manifest)
    panel_variant_ids = {variant.variant_id for variant in panel.variants}
    manifest_variant_ids = set(manifest_records)

    resolved_variant_ids = [
        variant_id
        for variant_id, annotation in annotations.items()
        if annotation.glytoucan_accession
        and annotation.atomistic_preparation_status == "registered_structure_resolved"
    ]
    unresolved_variant_ids = [
        variant_id for variant_id in panel_variant_ids if variant_id not in resolved_variant_ids
    ]
    accession_mismatches = [
        variant_id
        for variant_id, annotation in annotations.items()
        if (
            manifest_records.get(variant_id) is not None
            and _normalize_optional_text(
                panel_variant_accession(panel, variant_id)
            )
            and manifest_records[variant_id].source_accession
            and panel_variant_accession(panel, variant_id)
            != manifest_records[variant_id].source_accession
        )
    ]

    return {
        "panel_variant_count": len(panel.variants),
        "manifest_record_count": len(manifest.records),
        "matched_variant_count": len(panel_variant_ids & manifest_variant_ids),
        "missing_in_manifest": sorted(panel_variant_ids - manifest_variant_ids),
        "missing_in_panel": sorted(manifest_variant_ids - panel_variant_ids),
        "resolved_variant_ids": sorted(resolved_variant_ids),
        "unresolved_variant_ids": sorted(unresolved_variant_ids),
        "accession_mismatches": sorted(accession_mismatches),
        "readiness_fraction": _safe_fraction(len(resolved_variant_ids), len(panel.variants)),
    }


def panel_variant_accession(panel: HSVariantPanelConfig, variant_id: str) -> str | None:
    for variant in panel.variants:
        if variant.variant_id == variant_id:
            return _normalize_optional_text(variant.structural_annotations.get("glytoucan_accession"))
    return None


def build_manifest_validation_summary(
    *,
    hs_panel: HSVariantPanelConfig,
    hs_manifest: DataManifest,
    off_target_manifest: DataManifest,
    aggregate_manifest: DataManifest,
) -> str:
    hs_readiness = summarize_hs_panel_readiness(hs_panel, hs_manifest)
    hs_manifest_summary = summarize_manifest(hs_manifest)
    off_target_summary = summarize_manifest(off_target_manifest)
    aggregate_summary = summarize_manifest(aggregate_manifest)

    lines = [
        "# Manifest Validation Summary",
        "",
        "## HS Panel Alignment",
        "",
        f"- Panel variants: {hs_readiness['panel_variant_count']}",
        f"- Manifest records: {hs_readiness['manifest_record_count']}",
        f"- Matched IDs: {hs_readiness['matched_variant_count']}",
        f"- Resolved atomistic-ready variants: {len(hs_readiness['resolved_variant_ids'])}",
        f"- Unresolved variants: {len(hs_readiness['unresolved_variant_ids'])}",
        f"- Readiness fraction: {hs_readiness['readiness_fraction']:.3f}",
        f"- Missing in manifest: {', '.join(hs_readiness['missing_in_manifest']) or 'none'}",
        f"- Missing in panel: {', '.join(hs_readiness['missing_in_panel']) or 'none'}",
        f"- Accession mismatches: {', '.join(hs_readiness['accession_mismatches']) or 'none'}",
        "",
        "## Manifest Status Snapshots",
        "",
        f"- {hs_manifest_summary['manifest_name']} | statuses={hs_manifest_summary['by_status']}",
        f"- {off_target_summary['manifest_name']} | statuses={off_target_summary['by_status']}",
        f"- {aggregate_summary['manifest_name']} | statuses={aggregate_summary['by_status']}",
        "",
        "## Unresolved HS Variants",
        "",
    ]

    annotations = annotate_hs_panel_variants(hs_panel, hs_manifest)
    for variant_id in hs_readiness["unresolved_variant_ids"]:
        annotation = annotations[variant_id]
        accession = annotation.glytoucan_accession or "unassigned"
        lines.append(
            "- "
            f"{variant_id} | manifest_status={annotation.manifest_status} | "
            f"atomistic_status={annotation.atomistic_preparation_status} | accession={accession}"
        )

    return "\n".join(lines) + "\n"


def load_default_hs_resources() -> tuple[HSVariantPanelConfig, DataManifest]:
    return (
        load_hs_variant_panel(str(DEFAULT_HS_PANEL_PATH)),
        load_data_manifest(DEFAULT_HS_MANIFEST_PATH),
    )


def _normalize_optional_text(value: Any) -> str | None:
    if value is None:
        return None
    cleaned = str(value).strip()
    return cleaned or None


def _safe_fraction(numerator: int, denominator: int) -> float:
    if denominator <= 0:
        return 0.0
    return numerator / denominator
