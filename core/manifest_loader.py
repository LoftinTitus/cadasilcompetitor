"""Helpers for loading simple JSON-backed data manifests."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True, slots=True)
class ManifestRecord:
    record_id: str
    display_name: str
    entity_type: str
    source_db: str | None
    source_accession: str | None
    status: str
    tags: tuple[str, ...]
    attributes: dict[str, Any]
    notes: str | None


@dataclass(frozen=True, slots=True)
class DataManifest:
    manifest_id: str
    manifest_name: str
    entity_type: str
    source_status: str
    generated_by: str
    records: tuple[ManifestRecord, ...]


class ManifestValidationError(ValueError):
    """Raised when a JSON manifest is malformed."""


def load_data_manifest(path: str | Path) -> DataManifest:
    manifest_path = Path(path)
    with manifest_path.open(encoding="utf-8") as handle:
        raw_manifest = json.load(handle)

    if not isinstance(raw_manifest, dict):
        raise ManifestValidationError(
            f"Manifest at '{manifest_path}' must contain a top-level object."
        )

    return DataManifest(
        manifest_id=_read_required_string(raw_manifest, "manifest_id", manifest_path),
        manifest_name=_read_required_string(raw_manifest, "manifest_name", manifest_path),
        entity_type=_read_required_string(raw_manifest, "entity_type", manifest_path),
        source_status=_read_required_string(raw_manifest, "source_status", manifest_path),
        generated_by=_read_required_string(raw_manifest, "generated_by", manifest_path),
        records=tuple(
            _build_record(index, raw_record, manifest_path)
            for index, raw_record in enumerate(
                _read_required_list(raw_manifest, "records", manifest_path),
                start=1,
            )
        ),
    )


def summarize_manifest(manifest: DataManifest) -> dict[str, Any]:
    by_status: dict[str, int] = {}
    by_source_db: dict[str, int] = {}
    for record in manifest.records:
        by_status[record.status] = by_status.get(record.status, 0) + 1
        source_key = record.source_db or "unassigned"
        by_source_db[source_key] = by_source_db.get(source_key, 0) + 1

    return {
        "manifest_id": manifest.manifest_id,
        "manifest_name": manifest.manifest_name,
        "entity_type": manifest.entity_type,
        "record_count": len(manifest.records),
        "by_status": by_status,
        "by_source_db": by_source_db,
    }


def _build_record(index: int, raw_record: Any, manifest_path: Path) -> ManifestRecord:
    if not isinstance(raw_record, dict):
        raise ManifestValidationError(
            f"Record {index} in '{manifest_path}' must be an object."
        )

    tags = raw_record.get("tags", [])
    if not isinstance(tags, list) or any(not isinstance(tag, str) for tag in tags):
        raise ManifestValidationError(
            f"Record {index} in '{manifest_path}' has an invalid 'tags' field."
        )

    attributes = raw_record.get("attributes", {})
    if not isinstance(attributes, dict):
        raise ManifestValidationError(
            f"Record {index} in '{manifest_path}' has an invalid 'attributes' field."
        )

    notes = raw_record.get("notes")
    if notes is not None and not isinstance(notes, str):
        raise ManifestValidationError(
            f"Record {index} in '{manifest_path}' has an invalid 'notes' field."
        )

    return ManifestRecord(
        record_id=_read_required_string(raw_record, "record_id", manifest_path),
        display_name=_read_required_string(raw_record, "display_name", manifest_path),
        entity_type=_read_required_string(raw_record, "entity_type", manifest_path),
        source_db=_read_optional_string(raw_record.get("source_db"), "source_db", manifest_path),
        source_accession=_read_optional_string(
            raw_record.get("source_accession"),
            "source_accession",
            manifest_path,
        ),
        status=_read_required_string(raw_record, "status", manifest_path),
        tags=tuple(tags),
        attributes=attributes,
        notes=notes,
    )


def _read_required_string(raw_mapping: dict[str, Any], field_name: str, manifest_path: Path) -> str:
    value = raw_mapping.get(field_name)
    if not isinstance(value, str) or not value.strip():
        raise ManifestValidationError(
            f"Field '{field_name}' in '{manifest_path}' must be a non-empty string."
        )
    return value.strip()


def _read_optional_string(value: Any, field_name: str, manifest_path: Path) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise ManifestValidationError(
            f"Field '{field_name}' in '{manifest_path}' must be a string or null."
        )
    stripped = value.strip()
    return stripped or None


def _read_required_list(raw_mapping: dict[str, Any], field_name: str, manifest_path: Path) -> list[Any]:
    value = raw_mapping.get(field_name)
    if not isinstance(value, list):
        raise ManifestValidationError(
            f"Field '{field_name}' in '{manifest_path}' must be a list."
        )
    return value
