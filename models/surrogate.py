"""Bootstrap ridge-regression surrogates with simple uncertainty estimates."""

from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import Any

from .features import FEATURE_NAMES, extract_feature_vector


@dataclass(frozen=True, slots=True)
class FeatureScaler:
    means: tuple[float, ...]
    scales: tuple[float, ...]

    def transform(self, values: list[float]) -> list[float]:
        return [
            (value - mean) / scale
            for value, mean, scale in zip(values, self.means, self.scales)
        ]


@dataclass(frozen=True, slots=True)
class RidgeModel:
    scaler: FeatureScaler
    intercept: float
    coefficients: tuple[float, ...]

    def predict(self, candidate: dict[str, Any]) -> float:
        transformed = self.scaler.transform(extract_feature_vector(candidate))
        return self.intercept + sum(
            coefficient * value
            for coefficient, value in zip(self.coefficients, transformed)
        )


@dataclass(frozen=True, slots=True)
class BootstrapSurrogate:
    feature_names: tuple[str, ...]
    target_models: dict[str, tuple[RidgeModel, ...]]

    def predict(self, candidate: dict[str, Any]) -> dict[str, dict[str, float]]:
        predictions: dict[str, dict[str, float]] = {}
        for target_field, models in self.target_models.items():
            values = [model.predict(candidate) for model in models]
            mean_value = sum(values) / len(values)
            variance = sum((value - mean_value) ** 2 for value in values) / len(values)
            predictions[target_field] = {
                "mean": mean_value,
                "std": math.sqrt(max(variance, 0.0)),
                "min": min(values),
                "max": max(values),
            }
        return predictions


def fit_bootstrap_surrogate(
    rows: list[dict[str, Any]],
    *,
    target_fields: tuple[str, ...],
    ensemble_size: int = 25,
    ridge_alpha: float = 0.25,
    random_seed: int = 7,
) -> BootstrapSurrogate:
    if not rows:
        raise ValueError("At least one training row is required to fit a surrogate.")

    rng = random.Random(random_seed)
    target_models: dict[str, tuple[RidgeModel, ...]] = {}
    for target_field in target_fields:
        models: list[RidgeModel] = []
        for _ in range(max(ensemble_size, 1)):
            bootstrap_rows = [rng.choice(rows) for _ in range(max(len(rows), 1))]
            models.append(
                _fit_single_target_model(
                    bootstrap_rows,
                    target_field=target_field,
                    ridge_alpha=ridge_alpha,
                )
            )
        target_models[target_field] = tuple(models)

    return BootstrapSurrogate(
        feature_names=FEATURE_NAMES,
        target_models=target_models,
    )


def _fit_single_target_model(
    rows: list[dict[str, Any]],
    *,
    target_field: str,
    ridge_alpha: float,
) -> RidgeModel:
    features = [extract_feature_vector(row) for row in rows]
    targets = [float(row[target_field]) for row in rows]
    scaler = _fit_scaler(features)
    transformed = [scaler.transform(values) for values in features]

    dimension = len(FEATURE_NAMES) + 1
    normal_matrix = [[0.0 for _ in range(dimension)] for _ in range(dimension)]
    normal_rhs = [0.0 for _ in range(dimension)]

    for values, target in zip(transformed, targets):
        augmented = [1.0] + values
        for row_index in range(dimension):
            normal_rhs[row_index] += augmented[row_index] * target
            for column_index in range(dimension):
                normal_matrix[row_index][column_index] += (
                    augmented[row_index] * augmented[column_index]
                )

    for index in range(1, dimension):
        normal_matrix[index][index] += ridge_alpha

    solution = _solve_linear_system(normal_matrix, normal_rhs)
    return RidgeModel(
        scaler=scaler,
        intercept=solution[0],
        coefficients=tuple(solution[1:]),
    )


def _fit_scaler(features: list[list[float]]) -> FeatureScaler:
    means: list[float] = []
    scales: list[float] = []
    feature_count = len(FEATURE_NAMES)
    for feature_index in range(feature_count):
        values = [row[feature_index] for row in features]
        mean_value = sum(values) / len(values)
        variance = sum((value - mean_value) ** 2 for value in values) / len(values)
        scale = math.sqrt(variance) or 1.0
        means.append(mean_value)
        scales.append(scale)
    return FeatureScaler(means=tuple(means), scales=tuple(scales))


def _solve_linear_system(matrix: list[list[float]], rhs: list[float]) -> list[float]:
    dimension = len(rhs)
    augmented = [row[:] + [rhs_value] for row, rhs_value in zip(matrix, rhs)]

    for pivot_index in range(dimension):
        pivot_row = max(
            range(pivot_index, dimension),
            key=lambda row_index: abs(augmented[row_index][pivot_index]),
        )
        if abs(augmented[pivot_row][pivot_index]) < 1e-12:
            augmented[pivot_row][pivot_index] = 1e-12
        if pivot_row != pivot_index:
            augmented[pivot_index], augmented[pivot_row] = (
                augmented[pivot_row],
                augmented[pivot_index],
            )

        pivot_value = augmented[pivot_index][pivot_index]
        augmented[pivot_index] = [value / pivot_value for value in augmented[pivot_index]]

        for row_index in range(dimension):
            if row_index == pivot_index:
                continue
            factor = augmented[row_index][pivot_index]
            if factor == 0.0:
                continue
            augmented[row_index] = [
                current - factor * pivot
                for current, pivot in zip(augmented[row_index], augmented[pivot_index])
            ]

    return [augmented[index][-1] for index in range(dimension)]
