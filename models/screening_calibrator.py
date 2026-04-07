"""Learned calibration model for screening status and composite score."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from .surrogate import BootstrapSurrogate, fit_bootstrap_surrogate

SCREENING_MODEL_VERSION = 1
DEFAULT_SCREENING_CALIBRATOR_PATH = Path("reports/screening_calibrator.json")
SCREENING_STATUS_SCORE_FIELD = "screening_status_score"
COMPOSITE_SCORE_SCALED_FIELD = "composite_screen_score_scaled"
SCREENING_TARGET_FIELDS: tuple[str, ...] = (
    SCREENING_STATUS_SCORE_FIELD,
    COMPOSITE_SCORE_SCALED_FIELD,
)


@dataclass(frozen=True, slots=True)
class ScreeningCalibratorMetadata:
    training_row_count: int
    target_fields: tuple[str, ...]
    ensemble_size: int
    ridge_alpha: float
    random_seed: int
    training_source: str | None = None


@dataclass(frozen=True, slots=True)
class ScreeningCalibrator:
    surrogate: BootstrapSurrogate
    metadata: ScreeningCalibratorMetadata

    def predict(self, candidate: dict[str, Any]) -> dict[str, dict[str, float]]:
        raw = self.surrogate.predict(candidate)
        return {
            target_field: {
                key: _clamp01(float(value))
                for key, value in summary.items()
            }
            for target_field, summary in raw.items()
        }

    def predict_summary(self, candidate: dict[str, Any]) -> dict[str, float]:
        predictions = self.predict(candidate)
        status = predictions[SCREENING_STATUS_SCORE_FIELD]
        composite = predictions[COMPOSITE_SCORE_SCALED_FIELD]
        return {
            "ml_screening_score": float(status["mean"]),
            "ml_predicted_composite_score": 100.0 * float(composite["mean"]),
            "ml_screening_uncertainty": float(status["std"]),
        }

    def to_dict(self) -> dict[str, Any]:
        return {
            "model_type": "screening_calibrator",
            "model_version": SCREENING_MODEL_VERSION,
            "metadata": asdict(self.metadata),
            "surrogate": self.surrogate.to_dict(),
        }

    def write_json(self, path: str | Path) -> None:
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> ScreeningCalibrator:
        metadata_payload = payload["metadata"]
        return cls(
            surrogate=BootstrapSurrogate.from_dict(payload["surrogate"]),
            metadata=ScreeningCalibratorMetadata(
                training_row_count=int(metadata_payload["training_row_count"]),
                target_fields=tuple(str(value) for value in metadata_payload["target_fields"]),
                ensemble_size=int(metadata_payload["ensemble_size"]),
                ridge_alpha=float(metadata_payload["ridge_alpha"]),
                random_seed=int(metadata_payload["random_seed"]),
                training_source=metadata_payload.get("training_source"),
            ),
        )

    @classmethod
    def read_json(cls, path: str | Path) -> ScreeningCalibrator:
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        return cls.from_dict(payload)


def fit_screening_calibrator(
    rows: list[dict[str, Any]],
    *,
    ensemble_size: int = 25,
    ridge_alpha: float = 0.25,
    random_seed: int = 7,
    training_source: str | None = None,
) -> ScreeningCalibrator:
    training_rows = [_ensure_calibration_targets(row) for row in rows]
    surrogate = fit_bootstrap_surrogate(
        training_rows,
        target_fields=SCREENING_TARGET_FIELDS,
        ensemble_size=ensemble_size,
        ridge_alpha=ridge_alpha,
        random_seed=random_seed,
    )
    return ScreeningCalibrator(
        surrogate=surrogate,
        metadata=ScreeningCalibratorMetadata(
            training_row_count=len(training_rows),
            target_fields=SCREENING_TARGET_FIELDS,
            ensemble_size=ensemble_size,
            ridge_alpha=ridge_alpha,
            random_seed=random_seed,
            training_source=training_source,
        ),
    )


def load_screening_training_rows(path: str | Path) -> list[dict[str, Any]]:
    csv_path = Path(path)
    with csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = [_normalize_row(dict(row)) for row in reader]
    if not rows:
        raise ValueError("At least one training row is required to fit the screening calibrator.")
    return [_ensure_calibration_targets(row) for row in rows]


def try_load_screening_calibrator(
    path: str | Path = DEFAULT_SCREENING_CALIBRATOR_PATH,
) -> ScreeningCalibrator | None:
    model_path = Path(path)
    if not model_path.exists():
        return None
    return ScreeningCalibrator.read_json(model_path)


def _ensure_calibration_targets(row: dict[str, Any]) -> dict[str, Any]:
    if "sequence" not in row or not str(row.get("sequence", "")).strip():
        raise ValueError("Screening calibrator rows require a non-empty 'sequence'.")
    calibrated = dict(row)
    calibrated[SCREENING_STATUS_SCORE_FIELD] = _status_to_score(
        str(row.get("screening_status", ""))
    )
    composite_score = float(row.get("composite_screen_score", 0.0) or 0.0)
    calibrated[COMPOSITE_SCORE_SCALED_FIELD] = _clamp01(composite_score / 100.0)
    return calibrated


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
            normalized[key] = [item for item in value.split(";") if item]
            continue
        try:
            numeric_value = float(value)
        except ValueError:
            normalized[key] = value
            continue
        normalized[key] = int(numeric_value) if numeric_value.is_integer() else numeric_value
    return normalized


def _status_to_score(status: str) -> float:
    if status == "pass":
        return 1.0
    if status == "review":
        return 0.5
    if status == "reject":
        return 0.0
    return 0.25


def _clamp01(value: float) -> float:
    return max(0.0, min(value, 1.0))
