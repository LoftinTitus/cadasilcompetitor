"""Tier A occupancy model for HS-targeting peptide simulations."""

from __future__ import annotations

import math
from typing import Any

from core.config_loader import (
    CandidateProperties,
    DEFAULT_CONCENTRATION_GRID_UM,
    DEFAULT_OCCUPANCY_THRESHOLD_FRACTION,
    DEFAULT_TIME_STEP_S,
    HSVariantConfig,
    PhysiologyConfig,
    SimulationConfig,
)
from simulation.environment import (
    build_candidate_context,
    build_candidate_property_state,
    build_environment_state,
    build_hs_variant_context,
    build_target_context,
)

MICROMOLAR_TO_MOLAR = 1e-6
PATH_LENGTH_MULTIPLIER_BY_COMPARTMENT = {
    "glycocalyx": 1.0,
    "basement_membrane": 4.0,
    "perivascular_ecm": 8.0,
}
EXPOSURE_ACCESSIBILITY_BY_ROUTE = {
    ("glycocalyx", "luminal"): 1.0,
    ("glycocalyx", "abluminal"): 0.35,
    ("basement_membrane", "luminal"): 0.55,
    ("basement_membrane", "abluminal"): 0.75,
    ("perivascular_ecm", "luminal"): 0.2,
    ("perivascular_ecm", "abluminal"): 0.9,
}
WASHOUT_FACTOR_BY_COMPARTMENT = {
    "glycocalyx": 1.0,
    "basement_membrane": 0.55,
    "perivascular_ecm": 0.3,
}
EXPECTED_AFFINITY_SCALE = {
    "low": 0.65,
    "low_to_medium": 0.85,
    "medium": 1.0,
    "medium_to_high": 1.2,
    "high": 1.45,
}


def prepare_occupancy_model_inputs(
    candidate: Any,
    simulation_config: SimulationConfig | PhysiologyConfig,
    candidate_properties: CandidateProperties | None = None,
    hs_variant: HSVariantConfig | None = None,
) -> dict[str, Any]:
    physiology = (
        simulation_config.physiology
        if isinstance(simulation_config, SimulationConfig)
        else simulation_config
    )
    target_context = (
        build_target_context(simulation_config.target)
        if isinstance(simulation_config, SimulationConfig)
        else None
    )
    controls = _build_simulation_controls(simulation_config)
    candidate_property_state = build_candidate_property_state(candidate_properties)

    model_parameters = None
    if candidate_properties is not None:
        model_parameters = _build_effective_model_parameters(
            simulation_config,
            candidate_properties,
            hs_variant=hs_variant,
        )

    return {
        "candidate": build_candidate_context(candidate),
        "environment": build_environment_state(simulation_config),
        "target": target_context,
        "candidate_properties": candidate_property_state,
        "hs_variant": build_hs_variant_context(hs_variant),
        "simulation_controls": controls,
        "model_parameters": model_parameters,
        "exposure_duration_s": physiology.exposure_duration_s,
        "screening_focus": (
            target_context["screening_focus"] if target_context is not None else None
        ),
    }


def simulate_occupancy_sweep(model_inputs: dict[str, Any]) -> dict[str, Any]:
    model_parameters = model_inputs.get("model_parameters")
    if model_parameters is None:
        raise ValueError(
            "Candidate properties are required to simulate occupancy. "
            "Call prepare_occupancy_model_inputs() with CandidateProperties."
        )

    controls = model_inputs["simulation_controls"]
    threshold = float(controls["occupancy_threshold_fraction"])
    time_step_s = float(controls["time_step_s"])
    exposure_duration_s = float(model_inputs["exposure_duration_s"])

    concentration_results = [
        simulate_single_concentration(
            model_parameters=model_parameters,
            input_concentration_uM=float(concentration_uM),
            exposure_duration_s=exposure_duration_s,
            occupancy_threshold_fraction=threshold,
            time_step_s=time_step_s,
        )
        for concentration_uM in controls["concentration_grid_uM"]
    ]
    operating_point = _select_operating_point(
        concentration_results,
        occupancy_threshold_fraction=threshold,
    )

    return {
        "concentration_sweep": concentration_results,
        "operating_point": operating_point,
        "summary": {
            "screening_focus": model_inputs.get("screening_focus"),
            "occupancy_threshold_fraction": threshold,
            "recommended_concentration_uM": (
                operating_point["input_concentration_uM"] if operating_point is not None else None
            ),
            "lowest_concentration_meeting_threshold_uM": _lowest_concentration_meeting_threshold(
                concentration_results,
                threshold,
            ),
            "best_case_max_occupancy_fraction": max(
                result["summary"]["max_occupancy_fraction"] for result in concentration_results
            ),
            "best_time_above_threshold_s": max(
                result["summary"]["time_above_threshold_s"] for result in concentration_results
            ),
        },
    }


def simulate_single_concentration(
    *,
    model_parameters: dict[str, float],
    input_concentration_uM: float,
    exposure_duration_s: float,
    occupancy_threshold_fraction: float,
    time_step_s: float,
) -> dict[str, Any]:
    time_points_s = _build_time_grid(exposure_duration_s, time_step_s)
    source_concentration_M = input_concentration_uM * MICROMOLAR_TO_MOLAR

    delivery_rate_per_s = model_parameters["delivery_rate_per_s"]
    total_loss_rate_per_s = model_parameters["total_loss_rate_per_s"]
    effective_association_rate_M_inv_s = model_parameters[
        "effective_association_rate_M_inv_s"
    ]
    effective_dissociation_rate_per_s = model_parameters[
        "effective_dissociation_rate_per_s"
    ]
    maximum_occupancy_fraction = model_parameters["maximum_occupancy_fraction"]

    free_concentration_M = [
        _free_concentration_at_time(
            time_s,
            source_concentration_M=source_concentration_M,
            delivery_rate_per_s=delivery_rate_per_s,
            total_loss_rate_per_s=total_loss_rate_per_s,
        )
        for time_s in time_points_s
    ]

    occupancy_fraction = [0.0]
    for previous_time_s, next_time_s in zip(time_points_s, time_points_s[1:]):
        dt_s = next_time_s - previous_time_s
        midpoint_time_s = previous_time_s + (0.5 * dt_s)
        midpoint_concentration_M = _free_concentration_at_time(
            midpoint_time_s,
            source_concentration_M=source_concentration_M,
            delivery_rate_per_s=delivery_rate_per_s,
            total_loss_rate_per_s=total_loss_rate_per_s,
        )
        occupancy_fraction.append(
            _propagate_occupancy(
                occupancy_fraction[-1],
                dt_s=dt_s,
                concentration_M=midpoint_concentration_M,
                effective_association_rate_M_inv_s=effective_association_rate_M_inv_s,
                effective_dissociation_rate_per_s=effective_dissociation_rate_per_s,
                maximum_occupancy_fraction=maximum_occupancy_fraction,
            )
        )

    free_concentration_uM = [value / MICROMOLAR_TO_MOLAR for value in free_concentration_M]
    summary = _summarize_trajectory(
        time_points_s,
        free_concentration_uM,
        occupancy_fraction,
        occupancy_threshold_fraction=occupancy_threshold_fraction,
        estimated_kd_uM=model_parameters["estimated_kd_uM"],
        residence_time_s=model_parameters["residence_time_s"],
    )

    return {
        "input_concentration_uM": input_concentration_uM,
        "summary": summary,
        "trajectory": {
            "time_s": time_points_s,
            "free_concentration_uM": free_concentration_uM,
            "occupancy_fraction": occupancy_fraction,
        },
    }


def _build_simulation_controls(
    simulation_config: SimulationConfig | PhysiologyConfig,
) -> dict[str, Any]:
    if isinstance(simulation_config, SimulationConfig):
        controls = simulation_config.simulation_controls
        concentration_grid_uM = tuple(float(value) for value in controls.concentration_grid_uM)
        occupancy_threshold_fraction = float(controls.occupancy_threshold_fraction)
        time_step_s = float(controls.time_step_s)
    else:
        concentration_grid_uM = DEFAULT_CONCENTRATION_GRID_UM
        occupancy_threshold_fraction = DEFAULT_OCCUPANCY_THRESHOLD_FRACTION
        time_step_s = DEFAULT_TIME_STEP_S

    return {
        "concentration_grid_uM": concentration_grid_uM,
        "occupancy_threshold_fraction": occupancy_threshold_fraction,
        "time_step_s": time_step_s,
    }


def _build_effective_model_parameters(
    simulation_config: SimulationConfig | PhysiologyConfig,
    candidate_properties: CandidateProperties,
    *,
    hs_variant: HSVariantConfig | None,
) -> dict[str, float]:
    physiology = (
        simulation_config.physiology
        if isinstance(simulation_config, SimulationConfig)
        else simulation_config
    )

    if isinstance(simulation_config, SimulationConfig):
        target = simulation_config.target
        barrier = simulation_config.barrier
        transport = simulation_config.transport
        degradation = simulation_config.degradation
        neurovascular_environment = simulation_config.neurovascular_environment
        maximum_occupancy_fraction = simulation_config.binding.maximum_occupancy_fraction
    else:
        target = None
        barrier = None
        transport = None
        degradation = None
        neurovascular_environment = None
        maximum_occupancy_fraction = 1.0

    compartment = target.compartment if target is not None else "glycocalyx"
    exposure_side = target.exposure_side if target is not None else "luminal"
    route_factor = EXPOSURE_ACCESSIBILITY_BY_ROUTE[(compartment, exposure_side)]
    path_length_multiplier = PATH_LENGTH_MULTIPLIER_BY_COMPARTMENT[compartment]

    base_path_length_um = barrier.thickness_um if barrier is not None else 2.0
    effective_path_length_um = max(base_path_length_um * path_length_multiplier, 1e-6)
    porosity_fraction = barrier.porosity_fraction if barrier is not None else 0.2

    permeability_um_s = candidate_properties.barrier_permeability_cm_s * 1e4
    permeability_factor = (
        1.0
        if compartment == "glycocalyx" and exposure_side == "luminal"
        else permeability_um_s / (permeability_um_s + 0.1)
    )

    delivery_rate_per_s = (
        candidate_properties.diffusion_coefficient_um2_s
        * porosity_fraction
        * route_factor
        * max(permeability_factor, 1e-6)
        / (effective_path_length_um**2)
    )

    intrinsic_degradation_rate_per_s = (
        math.log(2.0) / candidate_properties.half_life_s
        if candidate_properties.half_life_s > 0.0
        else 0.0
    )
    degradation_rate_per_s = 0.0
    if degradation is None or degradation.degradation_enabled:
        degradation_rate_per_s = (
            intrinsic_degradation_rate_per_s
            + candidate_properties.enzymatic_degradation_rate_per_s
            + candidate_properties.spontaneous_degradation_rate_per_s
        )

    convective_washout_rate_per_s = 0.0
    if transport is not None and transport.convective_washout_enabled:
        convective_washout_rate_per_s = (
            (transport.flow_velocity_um_s * 2e-5) + (transport.shear_stress_Pa * 2e-3)
        ) * WASHOUT_FACTOR_BY_COMPARTMENT[compartment]

    pressure_loss_rate_per_s = 0.0
    if neurovascular_environment is not None:
        pressure_loss_rate_per_s = (
            neurovascular_environment.interstitial_pressure_mmHg * 1e-4
        ) * (0.5 if compartment == "glycocalyx" else 1.0)

    total_loss_rate_per_s = (
        candidate_properties.clearance_rate_per_s
        + degradation_rate_per_s
        + convective_washout_rate_per_s
        + pressure_loss_rate_per_s
    )

    variant_affinity_multiplier = _compute_variant_affinity_multiplier(hs_variant)
    affinity_rate_scale = math.sqrt(variant_affinity_multiplier)
    salt_sensitivity_factor = _compute_salt_sensitivity_factor(
        physiology.salt_concentration_mM
    )

    effective_association_rate_M_inv_s = (
        candidate_properties.association_rate_M_inv_s
        * affinity_rate_scale
        * salt_sensitivity_factor
    )
    effective_dissociation_rate_per_s = (
        candidate_properties.dissociation_rate_per_s / affinity_rate_scale
    )

    estimated_kd_uM = (
        (effective_dissociation_rate_per_s / effective_association_rate_M_inv_s) * 1e6
        if effective_association_rate_M_inv_s > 0.0
        else math.inf
    )
    residence_time_s = (
        math.inf if effective_dissociation_rate_per_s == 0.0 else 1.0 / effective_dissociation_rate_per_s
    )

    return {
        "delivery_rate_per_s": delivery_rate_per_s,
        "route_accessibility_factor": route_factor,
        "effective_path_length_um": effective_path_length_um,
        "permeability_factor": permeability_factor,
        "degradation_rate_per_s": degradation_rate_per_s,
        "convective_washout_rate_per_s": convective_washout_rate_per_s,
        "pressure_loss_rate_per_s": pressure_loss_rate_per_s,
        "total_loss_rate_per_s": total_loss_rate_per_s,
        "variant_affinity_multiplier": variant_affinity_multiplier,
        "salt_sensitivity_factor": salt_sensitivity_factor,
        "effective_association_rate_M_inv_s": effective_association_rate_M_inv_s,
        "effective_dissociation_rate_per_s": effective_dissociation_rate_per_s,
        "estimated_kd_uM": estimated_kd_uM,
        "residence_time_s": residence_time_s,
        "maximum_occupancy_fraction": maximum_occupancy_fraction,
    }


def _compute_variant_affinity_multiplier(hs_variant: HSVariantConfig | None) -> float:
    if hs_variant is None:
        return 1.0

    affinity_label = str(
        hs_variant.simulation_annotations.get("expected_relative_affinity", "medium")
    ).strip()
    label_scale = EXPECTED_AFFINITY_SCALE.get(affinity_label, 1.0)

    total_sulfates = float(hs_variant.sulfation_pattern.get("total_sulfates", 0) or 0)
    sulfation_scale = 1.0 + (0.06 * total_sulfates)

    degree_of_polymerization = float(
        hs_variant.chain_length.get("degree_of_polymerization", 4) or 4
    )
    length_scale = 1.0 + (0.04 * max(degree_of_polymerization - 4.0, 0.0))

    anticoagulant_sentinel_scale = (
        1.35 if hs_variant.sulfation_pattern.get("contains_3_o_sulfation") else 1.0
    )
    return _clamp(
        label_scale * sulfation_scale * length_scale * anticoagulant_sentinel_scale,
        0.25,
        6.0,
    )


def _compute_salt_sensitivity_factor(salt_concentration_mM: float) -> float:
    if salt_concentration_mM <= 0.0:
        return 1.0
    return _clamp(math.sqrt(150.0 / salt_concentration_mM), 0.5, 1.5)


def _free_concentration_at_time(
    time_s: float,
    *,
    source_concentration_M: float,
    delivery_rate_per_s: float,
    total_loss_rate_per_s: float,
) -> float:
    combined_rate_per_s = delivery_rate_per_s + total_loss_rate_per_s
    if combined_rate_per_s <= 0.0:
        return source_concentration_M

    steady_state_fraction = delivery_rate_per_s / combined_rate_per_s
    return source_concentration_M * steady_state_fraction * (
        1.0 - math.exp(-combined_rate_per_s * time_s)
    )


def _propagate_occupancy(
    previous_occupancy_fraction: float,
    *,
    dt_s: float,
    concentration_M: float,
    effective_association_rate_M_inv_s: float,
    effective_dissociation_rate_per_s: float,
    maximum_occupancy_fraction: float,
) -> float:
    association_flux_per_s = effective_association_rate_M_inv_s * concentration_M
    total_rate_per_s = association_flux_per_s + effective_dissociation_rate_per_s
    if total_rate_per_s <= 0.0:
        return previous_occupancy_fraction

    steady_state_fraction = (
        maximum_occupancy_fraction * association_flux_per_s / total_rate_per_s
    )
    next_occupancy_fraction = steady_state_fraction + (
        previous_occupancy_fraction - steady_state_fraction
    ) * math.exp(-total_rate_per_s * dt_s)
    return _clamp(next_occupancy_fraction, 0.0, maximum_occupancy_fraction)


def _summarize_trajectory(
    time_points_s: list[float],
    free_concentration_uM: list[float],
    occupancy_fraction: list[float],
    *,
    occupancy_threshold_fraction: float,
    estimated_kd_uM: float,
    residence_time_s: float,
) -> dict[str, Any]:
    max_occupancy_fraction = max(occupancy_fraction)
    occupancy_auc_fraction_s = _trapezoid_area(time_points_s, occupancy_fraction)
    exposure_duration_s = time_points_s[-1]
    time_to_threshold_s = _find_crossing_time(
        time_points_s,
        occupancy_fraction,
        occupancy_threshold_fraction,
    )
    time_to_half_max_s = _find_crossing_time(
        time_points_s,
        occupancy_fraction,
        max_occupancy_fraction * 0.5,
    )

    return {
        "max_occupancy_fraction": max_occupancy_fraction,
        "final_occupancy_fraction": occupancy_fraction[-1],
        "mean_occupancy_fraction": (
            occupancy_auc_fraction_s / exposure_duration_s if exposure_duration_s > 0.0 else 0.0
        ),
        "occupancy_auc_fraction_s": occupancy_auc_fraction_s,
        "time_above_threshold_s": _time_above_threshold(
            time_points_s,
            occupancy_fraction,
            occupancy_threshold_fraction,
        ),
        "time_to_threshold_s": time_to_threshold_s,
        "time_to_half_max_s": time_to_half_max_s,
        "peak_free_concentration_uM": max(free_concentration_uM),
        "steady_state_free_concentration_uM": free_concentration_uM[-1],
        "mean_free_concentration_uM": (
            _trapezoid_area(time_points_s, free_concentration_uM) / exposure_duration_s
            if exposure_duration_s > 0.0
            else 0.0
        ),
        "estimated_kd_uM": estimated_kd_uM,
        "residence_time_s": residence_time_s,
    }


def _build_time_grid(exposure_duration_s: float, time_step_s: float) -> list[float]:
    if exposure_duration_s <= 0.0:
        return [0.0]

    time_points_s = [0.0]
    current_time_s = 0.0
    while current_time_s + time_step_s < exposure_duration_s:
        current_time_s += time_step_s
        time_points_s.append(current_time_s)
    if time_points_s[-1] != exposure_duration_s:
        time_points_s.append(exposure_duration_s)
    return time_points_s


def _find_crossing_time(
    time_points_s: list[float],
    values: list[float],
    threshold: float,
) -> float | None:
    if threshold <= 0.0:
        return 0.0

    for left_time_s, right_time_s, left_value, right_value in zip(
        time_points_s,
        time_points_s[1:],
        values,
        values[1:],
    ):
        if left_value >= threshold:
            return left_time_s
        if left_value < threshold <= right_value:
            if right_value == left_value:
                return right_time_s
            fraction = (threshold - left_value) / (right_value - left_value)
            return left_time_s + (fraction * (right_time_s - left_time_s))
    return 0.0 if values and values[0] >= threshold else None


def _time_above_threshold(
    time_points_s: list[float],
    values: list[float],
    threshold: float,
) -> float:
    total_time_s = 0.0
    for left_time_s, right_time_s, left_value, right_value in zip(
        time_points_s,
        time_points_s[1:],
        values,
        values[1:],
    ):
        dt_s = right_time_s - left_time_s
        if left_value >= threshold and right_value >= threshold:
            total_time_s += dt_s
            continue
        if left_value < threshold and right_value < threshold:
            continue
        if right_value == left_value:
            continue

        crossing_fraction = (threshold - left_value) / (right_value - left_value)
        crossing_fraction = _clamp(crossing_fraction, 0.0, 1.0)
        crossing_time_s = left_time_s + (crossing_fraction * dt_s)
        if left_value < threshold <= right_value:
            total_time_s += right_time_s - crossing_time_s
        elif left_value >= threshold > right_value:
            total_time_s += crossing_time_s - left_time_s
    return total_time_s


def _trapezoid_area(time_points_s: list[float], values: list[float]) -> float:
    return sum(
        0.5 * (left_value + right_value) * (right_time_s - left_time_s)
        for left_time_s, right_time_s, left_value, right_value in zip(
            time_points_s,
            time_points_s[1:],
            values,
            values[1:],
        )
    )


def _select_operating_point(
    concentration_results: list[dict[str, Any]],
    *,
    occupancy_threshold_fraction: float,
) -> dict[str, Any] | None:
    threshold_hits = [
        result
        for result in concentration_results
        if result["summary"]["max_occupancy_fraction"] >= occupancy_threshold_fraction
    ]
    if threshold_hits:
        threshold_hits.sort(
            key=lambda result: (
                result["input_concentration_uM"],
                -result["summary"]["time_above_threshold_s"],
            )
        )
        return threshold_hits[0]

    if not concentration_results:
        return None

    return max(
        concentration_results,
        key=lambda result: (
            result["summary"]["max_occupancy_fraction"],
            result["summary"]["time_above_threshold_s"],
            -result["input_concentration_uM"],
        ),
    )


def _lowest_concentration_meeting_threshold(
    concentration_results: list[dict[str, Any]],
    occupancy_threshold_fraction: float,
) -> float | None:
    matches = [
        result["input_concentration_uM"]
        for result in concentration_results
        if result["summary"]["max_occupancy_fraction"] >= occupancy_threshold_fraction
    ]
    return min(matches) if matches else None


def _clamp(value: float, minimum: float, maximum: float) -> float:
    return max(minimum, min(value, maximum))
