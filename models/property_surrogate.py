"""Trainable surrogate for candidate transport and kinetic properties."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from core.config_loader import CandidateProperties

from .surrogate import BootstrapSurrogate, fit_bootstrap_surrogate

PROPERTY_MODEL_VERSION = 1
DEFAULT_PROPERTY_ESTIMATOR_PATH = Path("reports/property_estimator.json")
PROPERTY_TARGET_FIELDS: tuple[str, ...] = (
    "estimated_diffusion_coefficient_um2_s",
    "estimated_clearance_rate_per_s",
    "estimated_association_rate_M_inv_s",
    "estimated_dissociation_rate_per_s",
    "estimated_barrier_permeability_cm_s",
    "estimated_half_life_s",
    "estimated_enzymatic_degradation_rate_per_s",
    "estimated_spontaneous_degradation_rate_per_s",
)
PROPERTY_FIELD_BOUNDS: dict[str, tuple[float, float]] = {
    "estimated_diffusion_coefficient_um2_s": (20.0, 140.0),
    "estimated_clearance_rate_per_s": (1.0e-4, 0.01),
    "estimated_association_rate_M_inv_s": (2.0e4, 8.0e5),
    "estimated_dissociation_rate_per_s": (0.03, 1.5),
    "estimated_barrier_permeability_cm_s": (1.0e-8, 5.0e-6),
    "estimated_half_life_s": (900.0, 43200.0),
    "estimated_enzymatic_degradation_rate_per_s": (1.0e-5, 0.003),
    "estimated_spontaneous_degradation_rate_per_s": (1.0e-6, 5.0e-4),
}
PROPERTY_FIELD_TO_CANDIDATE_FIELD: dict[str, str] = {
    "estimated_diffusion_coefficient_um2_s": "diffusion_coefficient_um2_s",
    "estimated_clearance_rate_per_s": "clearance_rate_per_s",
    "estimated_association_rate_M_inv_s": "association_rate_M_inv_s",
    "estimated_dissociation_rate_per_s": "dissociation_rate_per_s",
    "estimated_barrier_permeability_cm_s": "barrier_permeability_cm_s",
    "estimated_half_life_s": "half_life_s",
    "estimated_enzymatic_degradation_rate_per_s": "enzymatic_degradation_rate_per_s",
    "estimated_spontaneous_degradation_rate_per_s": "spontaneous_degradation_rate_per_s",
}


@dataclass(frozen=True, slots=True)
class PropertySurrogateMetadata:
    training_row_count: int
    target_fields: tuple[str, ...]
    ensemble_size: int
    ridge_alpha: float
    random_seed: int
    training_source: str | None = None


@dataclass(frozen=True, slots=True)
class PropertySurrogate:
    surrogate: BootstrapSurrogate
    metadata: PropertySurrogateMetadata

    def predict(self, candidate: dict[str, Any]) -> dict[str, dict[str, float]]:
        predictions = self.surrogate.predict(candidate)
        bounded: dict[str, dict[str, float]] = {}
        for field_name, summary in predictions.items():
            minimum, maximum = PROPERTY_FIELD_BOUNDS[field_name]
            bounded[field_name] = {
                key: _clamp(float(value), minimum, maximum)
                for key, value in summary.items()
            }
        return bounded

    def predict_candidate_properties(self, candidate: dict[str, Any]) -> CandidateProperties:
        predictions = self.predict(candidate)
        return CandidateProperties(
            **{
                PROPERTY_FIELD_TO_CANDIDATE_FIELD[field_name]: float(predictions[field_name]["mean"])
                for field_name in PROPERTY_TARGET_FIELDS
            }
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "model_type": "property_surrogate",
            "model_version": PROPERTY_MODEL_VERSION,
            "metadata": asdict(self.metadata),
            "surrogate": self.surrogate.to_dict(),
        }

    def write_json(self, path: str | Path) -> None:
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> PropertySurrogate:
        metadata_payload = payload["metadata"]
        return cls(
            surrogate=BootstrapSurrogate.from_dict(payload["surrogate"]),
            metadata=PropertySurrogateMetadata(
                training_row_count=int(metadata_payload["training_row_count"]),
                target_fields=tuple(str(value) for value in metadata_payload["target_fields"]),
                ensemble_size=int(metadata_payload["ensemble_size"]),
                ridge_alpha=float(metadata_payload["ridge_alpha"]),
                random_seed=int(metadata_payload["random_seed"]),
                training_source=metadata_payload.get("training_source"),
            ),
        )

    @classmethod
    def read_json(cls, path: str | Path) -> PropertySurrogate:
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        return cls.from_dict(payload)


def fit_property_surrogate(
    rows: list[dict[str, Any]],
    *,
    ensemble_size: int = 25,
    ridge_alpha: float = 0.25,
    random_seed: int = 7,
    training_source: str | None = None,
) -> PropertySurrogate:
    _validate_training_rows(rows)
    surrogate = fit_bootstrap_surrogate(
        rows,
        target_fields=PROPERTY_TARGET_FIELDS,
        ensemble_size=ensemble_size,
        ridge_alpha=ridge_alpha,
        random_seed=random_seed,
    )
    return PropertySurrogate(
        surrogate=surrogate,
        metadata=PropertySurrogateMetadata(
            training_row_count=len(rows),
            target_fields=PROPERTY_TARGET_FIELDS,
            ensemble_size=ensemble_size,
            ridge_alpha=ridge_alpha,
            random_seed=random_seed,
            training_source=training_source,
        ),
    )


def load_property_training_rows(path: str | Path) -> list[dict[str, Any]]:
    csv_path = Path(path)
    with csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = [_normalize_row(dict(row)) for row in reader]
    _validate_training_rows(rows)
    return rows


def try_load_property_surrogate(path: str | Path = DEFAULT_PROPERTY_ESTIMATOR_PATH) -> PropertySurrogate | None:
    model_path = Path(path)
    if not model_path.exists():
        return None
    return PropertySurrogate.read_json(model_path)


def _normalize_row(row: dict[str, str]) -> dict[str, Any]:
    normalized: dict[str, Any] = {}
    for key, value in row.items():
        if value is None or value == "":
            normalized[key] = value
            continue
        if key in {"sequence", "candidate_id", "screening_status"}:
            normalized[key] = value
            continue
        if key in {"warning_flags", "filter_flags", "risk_flags"}:
            if value.strip() in {"", "[]"}:
                normalized[key] = []
            else:
                normalized[key] = [item for item in value.split(";") if item]
            continue
        try:
            numeric_value = float(value)
        except ValueError:
            normalized[key] = value
            continue
        normalized[key] = int(numeric_value) if numeric_value.is_integer() else numeric_value
    return normalized


def _validate_training_rows(rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("At least one training row is required to fit the property surrogate.")
    missing_fields = [
        field_name
        for field_name in PROPERTY_TARGET_FIELDS
        if any(field_name not in row or row[field_name] in {"", None} for row in rows)
    ]
    if missing_fields:
        missing = ", ".join(missing_fields)
        raise ValueError(f"Training rows are missing required property targets: {missing}")


def _clamp(value: float, minimum: float, maximum: float) -> float:
    return max(minimum, min(value, maximum))
