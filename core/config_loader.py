"""Load shared physiology configuration for simulation modules.

Expected YAML structure:

physiology:
  pH: 7.4
  salt_concentration_mM: 150.0
  temperature_K: 310.15
  exposure_duration_s: 3600.0

Required fields and units:
- pH: dimensionless
- salt_concentration_mM: millimolar
- temperature_K: Kelvin
- exposure_duration_s: seconds

To extend the configuration, add new keys under the top-level ``physiology``
mapping and then add matching fields to ``PhysiologyConfig`` when simulation
code is ready to consume them. Unknown keys are ignored by the current loader,
so adding future parameters does not break existing modules.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True, slots=True)
class PhysiologyConfig:
    pH: float
    salt_concentration_mM: float
    temperature_K: float
    exposure_duration_s: float


class ConfigValidationError(ValueError):
    """Raised when a physiology config file is missing required fields or values."""


def load_physiology_config(path: str) -> PhysiologyConfig:
    config_path = Path(path)
    raw_config = _load_yaml_mapping(config_path)

    if not isinstance(raw_config, dict):
        raise ConfigValidationError(
            f"Configuration at '{config_path}' must contain a top-level mapping."
        )

    physiology_section = raw_config.get("physiology")
    if not isinstance(physiology_section, dict):
        raise ConfigValidationError(
            f"Configuration at '{config_path}' must contain a 'physiology' mapping."
        )

    missing_fields = [
        field_name
        for field_name in (
            "pH",
            "salt_concentration_mM",
            "temperature_K",
            "exposure_duration_s",
        )
        if field_name not in physiology_section
    ]
    if missing_fields:
        missing = ", ".join(missing_fields)
        raise ConfigValidationError(
            f"Missing required physiology field(s) in '{config_path}': {missing}."
        )

    physiology = PhysiologyConfig(
        pH=_coerce_float(physiology_section["pH"], "physiology.pH"),
        salt_concentration_mM=_coerce_float(
            physiology_section["salt_concentration_mM"],
            "physiology.salt_concentration_mM",
        ),
        temperature_K=_coerce_float(
            physiology_section["temperature_K"],
            "physiology.temperature_K",
        ),
        exposure_duration_s=_coerce_float(
            physiology_section["exposure_duration_s"],
            "physiology.exposure_duration_s",
        ),
    )
    _validate_physiology_config(physiology)
    return physiology


def _coerce_float(value: Any, field_name: str) -> float:
    if isinstance(value, bool):
        raise ConfigValidationError(
            f"Field '{field_name}' must be numeric, but received a boolean."
        )

    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ConfigValidationError(
            f"Field '{field_name}' must be numeric, but received {value!r}."
        ) from exc


def _validate_physiology_config(physiology: PhysiologyConfig) -> None:
    if not 0.0 <= physiology.pH <= 14.0:
        raise ConfigValidationError("Field 'physiology.pH' must be between 0 and 14.")
    if physiology.salt_concentration_mM <= 0.0:
        raise ConfigValidationError(
            "Field 'physiology.salt_concentration_mM' must be positive."
        )
    if physiology.temperature_K <= 0.0:
        raise ConfigValidationError("Field 'physiology.temperature_K' must be positive.")
    if physiology.exposure_duration_s < 0.0:
        raise ConfigValidationError(
            "Field 'physiology.exposure_duration_s' must be greater than or equal to 0."
        )


def _load_yaml_mapping(path: Path) -> dict[str, Any]:
    try:
        import yaml  # type: ignore[import-not-found]
    except ModuleNotFoundError:
        return _load_simple_yaml_mapping(path)

    with path.open("r", encoding="utf-8") as handle:
        parsed = yaml.safe_load(handle)

    if parsed is None:
        raise ConfigValidationError(f"Configuration file '{path}' is empty.")
    if not isinstance(parsed, dict):
        raise ConfigValidationError(f"Configuration file '{path}' must contain a mapping.")
    return parsed


def _load_simple_yaml_mapping(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"Configuration file not found: '{path}'.")

    root: dict[str, Any] = {}
    stack: list[tuple[int, dict[str, Any]]] = [(-1, root)]

    for line_number, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        line = _strip_yaml_comment(raw_line).rstrip()
        if not line.strip():
            continue

        indent = len(line) - len(line.lstrip(" "))
        if indent % 2 != 0:
            raise ConfigValidationError(
                f"Invalid indentation in '{path}' on line {line_number}: use multiples of two spaces."
            )

        content = line.lstrip(" ")
        if ":" not in content:
            raise ConfigValidationError(
                f"Invalid YAML entry in '{path}' on line {line_number}: expected 'key: value'."
            )

        key, raw_value = content.split(":", maxsplit=1)
        key = key.strip()
        raw_value = raw_value.strip()
        if not key:
            raise ConfigValidationError(
                f"Invalid YAML entry in '{path}' on line {line_number}: missing key name."
            )

        while len(stack) > 1 and indent <= stack[-1][0]:
            stack.pop()

        parent = stack[-1][1]
        if raw_value == "":
            child: dict[str, Any] = {}
            parent[key] = child
            stack.append((indent, child))
            continue

        parent[key] = _parse_scalar(raw_value)

    if not root:
        raise ConfigValidationError(f"Configuration file '{path}' is empty.")
    return root


def _strip_yaml_comment(raw_line: str) -> str:
    for index, character in enumerate(raw_line):
        if character == "#":
            return raw_line[:index]
    return raw_line


def _parse_scalar(raw_value: str) -> Any:
    if len(raw_value) >= 2 and raw_value[0] == raw_value[-1] and raw_value[0] in {'"', "'"}:
        return raw_value[1:-1]

    normalized = raw_value.lower()
    if normalized == "true":
        return True
    if normalized == "false":
        return False
    if normalized in {"null", "none"}:
        return None

    try:
        if any(marker in raw_value for marker in (".", "e", "E")):
            return float(raw_value)
        return int(raw_value)
    except ValueError:
        return raw_value
