"""Load shared simulation configuration and validate candidate-specific inputs.

The YAML file stores only environment and assay settings that are shared across
candidate peptides. Properties that vary per candidate should be computed during
testing or prediction and passed separately as ``CandidateProperties``.

Expected YAML structure:

physiology:
  pH: 7.4
  salt_concentration_mM: 150.0
  temperature_K: 310.15
  exposure_duration_s: 3600.0

target:
  compartment: glycocalyx
  exposure_side: luminal
  biological_rationale_label: endothelial_surface_hspg_engagement

transport:
  flow_velocity_um_s: 100.0
  shear_stress_Pa: 1.5
  convective_washout_enabled: true

binding:
  binding_site_density_per_um2: 1000.0
  maximum_occupancy_fraction: 1.0

barrier:
  thickness_um: 2.0
  porosity_fraction: 0.2
  reflection_coefficient: 0.5

degradation:
  degradation_enabled: true

neurovascular_environment:
  extracellular_matrix_density_fraction: 0.15
  hspg_site_density_per_um3: 500.0
  interstitial_pressure_mmHg: 2.0
  oxygen_fraction: 0.21

Shared YAML fields and units:
- physiology.pH: dimensionless
- physiology.salt_concentration_mM: millimolar
- physiology.temperature_K: Kelvin
- physiology.exposure_duration_s: seconds
- target.compartment: one of glycocalyx, basement_membrane, perivascular_ecm
- target.exposure_side: one of luminal, abluminal
- target.biological_rationale_label: string label for screening intent
- transport.flow_velocity_um_s: micrometers per second
- transport.shear_stress_Pa: Pascal
- transport.convective_washout_enabled: boolean flag
- binding.binding_site_density_per_um2: sites per square micrometer
- binding.maximum_occupancy_fraction: unitless fraction
- barrier.thickness_um: micrometers
- barrier.porosity_fraction: unitless fraction
- barrier.reflection_coefficient: unitless fraction
- degradation.degradation_enabled: boolean flag
- neurovascular_environment.extracellular_matrix_density_fraction: unitless fraction
- neurovascular_environment.hspg_site_density_per_um3: sites per cubic micrometer
- neurovascular_environment.interstitial_pressure_mmHg: millimeters of mercury
- neurovascular_environment.oxygen_fraction: unitless fraction

Candidate-specific computed properties and units:
- diffusion_coefficient_um2_s: square micrometers per second
- clearance_rate_per_s: inverse seconds
- association_rate_M_inv_s: inverse molar inverse seconds
- dissociation_rate_per_s: inverse seconds
- barrier_permeability_cm_s: centimeters per second
- half_life_s: seconds
- enzymatic_degradation_rate_per_s: inverse seconds
- spontaneous_degradation_rate_per_s: inverse seconds
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

SUPPORTED_TARGET_COMPARTMENTS: tuple[str, ...] = (
    "glycocalyx",
    "basement_membrane",
    "perivascular_ecm",
)
SUPPORTED_EXPOSURE_SIDES: tuple[str, ...] = ("luminal", "abluminal")


@dataclass(frozen=True, slots=True)
class PhysiologyConfig:
    pH: float
    salt_concentration_mM: float
    temperature_K: float
    exposure_duration_s: float


@dataclass(frozen=True, slots=True)
class TargetConfig:
    compartment: str
    exposure_side: str
    biological_rationale_label: str


@dataclass(frozen=True, slots=True)
class TransportConfig:
    flow_velocity_um_s: float
    shear_stress_Pa: float
    convective_washout_enabled: bool


@dataclass(frozen=True, slots=True)
class BindingConfig:
    binding_site_density_per_um2: float
    maximum_occupancy_fraction: float


@dataclass(frozen=True, slots=True)
class BarrierConfig:
    thickness_um: float
    porosity_fraction: float
    reflection_coefficient: float


@dataclass(frozen=True, slots=True)
class DegradationConfig:
    degradation_enabled: bool


@dataclass(frozen=True, slots=True)
class NeurovascularEnvironmentConfig:
    extracellular_matrix_density_fraction: float
    hspg_site_density_per_um3: float
    interstitial_pressure_mmHg: float
    oxygen_fraction: float


@dataclass(frozen=True, slots=True)
class SimulationConfig:
    physiology: PhysiologyConfig
    target: TargetConfig
    transport: TransportConfig
    binding: BindingConfig
    barrier: BarrierConfig
    degradation: DegradationConfig
    neurovascular_environment: NeurovascularEnvironmentConfig


@dataclass(frozen=True, slots=True)
class CandidateProperties:
    diffusion_coefficient_um2_s: float
    clearance_rate_per_s: float
    association_rate_M_inv_s: float
    dissociation_rate_per_s: float
    barrier_permeability_cm_s: float
    half_life_s: float
    enzymatic_degradation_rate_per_s: float
    spontaneous_degradation_rate_per_s: float


class ConfigValidationError(ValueError):
    """Raised when configuration values or candidate properties are invalid."""


def load_simulation_config(path: str) -> SimulationConfig:
    config_path = Path(path)
    raw_config = _load_root_mapping(config_path)

    return SimulationConfig(
        physiology=_build_physiology_config(
            _require_mapping(raw_config, "physiology", config_path),
            config_path,
        ),
        target=_build_target_config(
            _require_mapping(raw_config, "target", config_path),
            config_path,
        ),
        transport=_build_transport_config(
            _require_mapping(raw_config, "transport", config_path),
            config_path,
        ),
        binding=_build_binding_config(
            _require_mapping(raw_config, "binding", config_path),
            config_path,
        ),
        barrier=_build_barrier_config(
            _require_mapping(raw_config, "barrier", config_path),
            config_path,
        ),
        degradation=_build_degradation_config(
            _require_mapping(raw_config, "degradation", config_path),
            config_path,
        ),
        neurovascular_environment=_build_neurovascular_environment_config(
            _require_mapping(raw_config, "neurovascular_environment", config_path),
            config_path,
        ),
    )


def load_physiology_config(path: str) -> PhysiologyConfig:
    config_path = Path(path)
    raw_config = _load_root_mapping(config_path)
    return _build_physiology_config(
        _require_mapping(raw_config, "physiology", config_path),
        config_path,
    )


def validate_candidate_properties(
    candidate_properties: CandidateProperties,
) -> CandidateProperties:
    if candidate_properties.diffusion_coefficient_um2_s <= 0.0:
        raise ConfigValidationError(
            "Field 'candidate_properties.diffusion_coefficient_um2_s' must be positive."
        )
    if candidate_properties.clearance_rate_per_s < 0.0:
        raise ConfigValidationError(
            "Field 'candidate_properties.clearance_rate_per_s' must be greater than or equal to 0."
        )
    if candidate_properties.association_rate_M_inv_s <= 0.0:
        raise ConfigValidationError(
            "Field 'candidate_properties.association_rate_M_inv_s' must be positive."
        )
    if candidate_properties.dissociation_rate_per_s < 0.0:
        raise ConfigValidationError(
            "Field 'candidate_properties.dissociation_rate_per_s' must be greater than or equal to 0."
        )
    if candidate_properties.barrier_permeability_cm_s < 0.0:
        raise ConfigValidationError(
            "Field 'candidate_properties.barrier_permeability_cm_s' must be greater than or equal to 0."
        )
    if candidate_properties.half_life_s <= 0.0:
        raise ConfigValidationError(
            "Field 'candidate_properties.half_life_s' must be positive."
        )
    if candidate_properties.enzymatic_degradation_rate_per_s < 0.0:
        raise ConfigValidationError(
            "Field 'candidate_properties.enzymatic_degradation_rate_per_s' must be greater than or equal to 0."
        )
    if candidate_properties.spontaneous_degradation_rate_per_s < 0.0:
        raise ConfigValidationError(
            "Field 'candidate_properties.spontaneous_degradation_rate_per_s' must be greater than or equal to 0."
        )
    return candidate_properties


def _load_root_mapping(path: Path) -> dict[str, Any]:
    raw_config = _load_yaml_mapping(path)
    if not isinstance(raw_config, dict):
        raise ConfigValidationError(
            f"Configuration at '{path}' must contain a top-level mapping."
        )
    return raw_config


def _require_mapping(
    raw_config: dict[str, Any],
    section_name: str,
    config_path: Path,
) -> dict[str, Any]:
    section = raw_config.get(section_name)
    if not isinstance(section, dict):
        raise ConfigValidationError(
            f"Configuration at '{config_path}' must contain a '{section_name}' mapping."
        )
    return section


def _require_fields(
    section: dict[str, Any],
    section_name: str,
    required_fields: tuple[str, ...],
    config_path: Path,
) -> None:
    missing_fields = [field_name for field_name in required_fields if field_name not in section]
    if missing_fields:
        missing = ", ".join(missing_fields)
        raise ConfigValidationError(
            f"Missing required field(s) in '{section_name}' section of '{config_path}': {missing}."
        )


def _build_physiology_config(
    section: dict[str, Any],
    config_path: Path,
) -> PhysiologyConfig:
    _require_fields(
        section,
        "physiology",
        (
            "pH",
            "salt_concentration_mM",
            "temperature_K",
            "exposure_duration_s",
        ),
        config_path,
    )

    physiology = PhysiologyConfig(
        pH=_coerce_float(section["pH"], "physiology.pH"),
        salt_concentration_mM=_coerce_float(
            section["salt_concentration_mM"],
            "physiology.salt_concentration_mM",
        ),
        temperature_K=_coerce_float(
            section["temperature_K"],
            "physiology.temperature_K",
        ),
        exposure_duration_s=_coerce_float(
            section["exposure_duration_s"],
            "physiology.exposure_duration_s",
        ),
    )
    _validate_physiology_config(physiology)
    return physiology


def _build_target_config(
    section: dict[str, Any],
    config_path: Path,
) -> TargetConfig:
    _require_fields(
        section,
        "target",
        (
            "compartment",
            "exposure_side",
            "biological_rationale_label",
        ),
        config_path,
    )

    target = TargetConfig(
        compartment=_coerce_string(section["compartment"], "target.compartment"),
        exposure_side=_coerce_string(section["exposure_side"], "target.exposure_side"),
        biological_rationale_label=_coerce_string(
            section["biological_rationale_label"],
            "target.biological_rationale_label",
        ),
    )
    _validate_target_config(target)
    return target


def _build_transport_config(
    section: dict[str, Any],
    config_path: Path,
) -> TransportConfig:
    _require_fields(
        section,
        "transport",
        (
            "flow_velocity_um_s",
            "shear_stress_Pa",
            "convective_washout_enabled",
        ),
        config_path,
    )

    transport = TransportConfig(
        flow_velocity_um_s=_coerce_float(
            section["flow_velocity_um_s"],
            "transport.flow_velocity_um_s",
        ),
        shear_stress_Pa=_coerce_float(
            section["shear_stress_Pa"],
            "transport.shear_stress_Pa",
        ),
        convective_washout_enabled=_coerce_bool(
            section["convective_washout_enabled"],
            "transport.convective_washout_enabled",
        ),
    )
    _validate_transport_config(transport)
    return transport


def _build_binding_config(
    section: dict[str, Any],
    config_path: Path,
) -> BindingConfig:
    _require_fields(
        section,
        "binding",
        (
            "binding_site_density_per_um2",
            "maximum_occupancy_fraction",
        ),
        config_path,
    )

    binding = BindingConfig(
        binding_site_density_per_um2=_coerce_float(
            section["binding_site_density_per_um2"],
            "binding.binding_site_density_per_um2",
        ),
        maximum_occupancy_fraction=_coerce_float(
            section["maximum_occupancy_fraction"],
            "binding.maximum_occupancy_fraction",
        ),
    )
    _validate_binding_config(binding)
    return binding


def _build_barrier_config(
    section: dict[str, Any],
    config_path: Path,
) -> BarrierConfig:
    _require_fields(
        section,
        "barrier",
        (
            "thickness_um",
            "porosity_fraction",
            "reflection_coefficient",
        ),
        config_path,
    )

    barrier = BarrierConfig(
        thickness_um=_coerce_float(section["thickness_um"], "barrier.thickness_um"),
        porosity_fraction=_coerce_float(
            section["porosity_fraction"],
            "barrier.porosity_fraction",
        ),
        reflection_coefficient=_coerce_float(
            section["reflection_coefficient"],
            "barrier.reflection_coefficient",
        ),
    )
    _validate_barrier_config(barrier)
    return barrier


def _build_degradation_config(
    section: dict[str, Any],
    config_path: Path,
) -> DegradationConfig:
    _require_fields(
        section,
        "degradation",
        ("degradation_enabled",),
        config_path,
    )

    return DegradationConfig(
        degradation_enabled=_coerce_bool(
            section["degradation_enabled"],
            "degradation.degradation_enabled",
        ),
    )


def _build_neurovascular_environment_config(
    section: dict[str, Any],
    config_path: Path,
) -> NeurovascularEnvironmentConfig:
    _require_fields(
        section,
        "neurovascular_environment",
        (
            "extracellular_matrix_density_fraction",
            "hspg_site_density_per_um3",
            "interstitial_pressure_mmHg",
            "oxygen_fraction",
        ),
        config_path,
    )

    neurovascular_environment = NeurovascularEnvironmentConfig(
        extracellular_matrix_density_fraction=_coerce_float(
            section["extracellular_matrix_density_fraction"],
            "neurovascular_environment.extracellular_matrix_density_fraction",
        ),
        hspg_site_density_per_um3=_coerce_float(
            section["hspg_site_density_per_um3"],
            "neurovascular_environment.hspg_site_density_per_um3",
        ),
        interstitial_pressure_mmHg=_coerce_float(
            section["interstitial_pressure_mmHg"],
            "neurovascular_environment.interstitial_pressure_mmHg",
        ),
        oxygen_fraction=_coerce_float(
            section["oxygen_fraction"],
            "neurovascular_environment.oxygen_fraction",
        ),
    )
    _validate_neurovascular_environment_config(neurovascular_environment)
    return neurovascular_environment


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


def _coerce_bool(value: Any, field_name: str) -> bool:
    if isinstance(value, bool):
        return value

    raise ConfigValidationError(
        f"Field '{field_name}' must be a boolean, but received {value!r}."
    )


def _coerce_string(value: Any, field_name: str) -> str:
    if not isinstance(value, str):
        raise ConfigValidationError(
            f"Field '{field_name}' must be a string, but received {value!r}."
        )

    stripped_value = value.strip()
    if not stripped_value:
        raise ConfigValidationError(f"Field '{field_name}' must not be empty.")
    return stripped_value


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


def _validate_target_config(target: TargetConfig) -> None:
    if target.compartment not in SUPPORTED_TARGET_COMPARTMENTS:
        supported = ", ".join(SUPPORTED_TARGET_COMPARTMENTS)
        raise ConfigValidationError(
            f"Field 'target.compartment' must be one of: {supported}."
        )
    if target.exposure_side not in SUPPORTED_EXPOSURE_SIDES:
        supported = ", ".join(SUPPORTED_EXPOSURE_SIDES)
        raise ConfigValidationError(
            f"Field 'target.exposure_side' must be one of: {supported}."
        )


def _validate_transport_config(transport: TransportConfig) -> None:
    if transport.flow_velocity_um_s < 0.0:
        raise ConfigValidationError(
            "Field 'transport.flow_velocity_um_s' must be greater than or equal to 0."
        )
    if transport.shear_stress_Pa < 0.0:
        raise ConfigValidationError(
            "Field 'transport.shear_stress_Pa' must be greater than or equal to 0."
        )


def _validate_binding_config(binding: BindingConfig) -> None:
    if binding.binding_site_density_per_um2 <= 0.0:
        raise ConfigValidationError(
            "Field 'binding.binding_site_density_per_um2' must be positive."
        )
    if not 0.0 <= binding.maximum_occupancy_fraction <= 1.0:
        raise ConfigValidationError(
            "Field 'binding.maximum_occupancy_fraction' must be between 0 and 1."
        )


def _validate_barrier_config(barrier: BarrierConfig) -> None:
    if barrier.thickness_um <= 0.0:
        raise ConfigValidationError("Field 'barrier.thickness_um' must be positive.")
    if not 0.0 <= barrier.porosity_fraction <= 1.0:
        raise ConfigValidationError(
            "Field 'barrier.porosity_fraction' must be between 0 and 1."
        )
    if not 0.0 <= barrier.reflection_coefficient <= 1.0:
        raise ConfigValidationError(
            "Field 'barrier.reflection_coefficient' must be between 0 and 1."
        )


def _validate_neurovascular_environment_config(
    neurovascular_environment: NeurovascularEnvironmentConfig,
) -> None:
    if not 0.0 <= neurovascular_environment.extracellular_matrix_density_fraction <= 1.0:
        raise ConfigValidationError(
            "Field 'neurovascular_environment.extracellular_matrix_density_fraction' must be between 0 and 1."
        )
    if neurovascular_environment.hspg_site_density_per_um3 < 0.0:
        raise ConfigValidationError(
            "Field 'neurovascular_environment.hspg_site_density_per_um3' must be greater than or equal to 0."
        )
    if neurovascular_environment.interstitial_pressure_mmHg < 0.0:
        raise ConfigValidationError(
            "Field 'neurovascular_environment.interstitial_pressure_mmHg' must be greater than or equal to 0."
        )
    if not 0.0 <= neurovascular_environment.oxygen_fraction <= 1.0:
        raise ConfigValidationError(
            "Field 'neurovascular_environment.oxygen_fraction' must be between 0 and 1."
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
