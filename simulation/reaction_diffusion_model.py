"""Tier B 1D reaction-diffusion transport model."""

from __future__ import annotations

import math
from typing import Any

from simulation.occupancy_model import (
    MICROMOLAR_TO_MOLAR,
    _build_time_grid,
    _clamp,
    _find_crossing_time,
    _time_above_threshold,
    _trapezoid_area,
)


def simulate_reaction_diffusion_sweep(model_inputs: dict[str, Any]) -> dict[str, Any]:
    model_parameters = model_inputs.get("model_parameters")
    if model_parameters is None:
        raise ValueError(
            "Candidate properties are required to simulate Tier B reaction-diffusion transport."
        )

    controls = model_inputs["simulation_controls"]
    threshold = float(controls["occupancy_threshold_fraction"])
    time_step_s = float(controls["time_step_s"])
    exposure_duration_s = float(model_inputs["exposure_duration_s"])
    spatial_grid_points = int(controls["spatial_grid_points"])

    concentration_results = [
        simulate_reaction_diffusion_single_concentration(
            model_parameters=model_parameters,
            input_concentration_uM=float(concentration_uM),
            exposure_duration_s=exposure_duration_s,
            occupancy_threshold_fraction=threshold,
            time_step_s=time_step_s,
            spatial_grid_points=spatial_grid_points,
        )
        for concentration_uM in controls["concentration_grid_uM"]
    ]
    operating_point = _select_reaction_diffusion_operating_point(
        concentration_results,
        occupancy_threshold_fraction=threshold,
    )

    return {
        "concentration_sweep": concentration_results,
        "operating_point": operating_point,
        "summary": {
            "screening_focus": model_inputs.get("screening_focus"),
            "model_tier": "B",
            "occupancy_threshold_fraction": threshold,
            "recommended_concentration_uM": (
                operating_point["input_concentration_uM"] if operating_point is not None else None
            ),
            "lowest_concentration_meeting_threshold_uM": _lowest_concentration_meeting_threshold(
                concentration_results,
                threshold,
            ),
            "best_case_max_mean_occupancy_fraction": max(
                result["summary"]["max_mean_occupancy_fraction"]
                for result in concentration_results
            ),
            "best_case_max_distal_occupancy_fraction": max(
                result["summary"]["max_distal_occupancy_fraction"]
                for result in concentration_results
            ),
            "best_case_max_penetration_depth_um": max(
                result["summary"]["max_penetration_depth_um_at_threshold"]
                for result in concentration_results
            ),
            "selection_peak_occupancy_fraction": max(
                result["summary"]["max_mean_occupancy_fraction"]
                for result in concentration_results
            ),
        },
    }


def simulate_reaction_diffusion_single_concentration(
    *,
    model_parameters: dict[str, float],
    input_concentration_uM: float,
    exposure_duration_s: float,
    occupancy_threshold_fraction: float,
    time_step_s: float,
    spatial_grid_points: int,
) -> dict[str, Any]:
    domain_length_um = max(model_parameters["effective_path_length_um"], 1e-6)
    grid_points = max(spatial_grid_points, 3)
    dx_um = domain_length_um / (grid_points - 1)
    depth_um = [index * dx_um for index in range(grid_points)]
    time_points_s = _build_time_grid(exposure_duration_s, time_step_s)

    source_concentration_M = input_concentration_uM * MICROMOLAR_TO_MOLAR
    effective_diffusion_coefficient_um2_s = model_parameters[
        "effective_diffusion_coefficient_um2_s"
    ]
    total_loss_rate_per_s = model_parameters["total_loss_rate_per_s"]
    effective_association_rate_M_inv_s = model_parameters[
        "effective_association_rate_M_inv_s"
    ]
    effective_dissociation_rate_per_s = model_parameters[
        "effective_dissociation_rate_per_s"
    ]
    binding_site_capacity_M = (
        model_parameters["binding_site_capacity_uM"] * MICROMOLAR_TO_MOLAR
    )

    free_concentration_M = [0.0] * grid_points
    free_concentration_M[0] = source_concentration_M
    bound_concentration_M = [0.0] * grid_points

    mean_free_concentration_uM = [sum(free_concentration_M) / grid_points / MICROMOLAR_TO_MOLAR]
    surface_free_concentration_uM = [free_concentration_M[0] / MICROMOLAR_TO_MOLAR]
    distal_free_concentration_uM = [free_concentration_M[-1] / MICROMOLAR_TO_MOLAR]
    mean_occupancy_fraction = [
        _mean_occupancy_fraction(bound_concentration_M, binding_site_capacity_M)
    ]
    surface_occupancy_fraction = [
        _occupancy_fraction(bound_concentration_M[0], binding_site_capacity_M)
    ]
    distal_occupancy_fraction = [
        _occupancy_fraction(bound_concentration_M[-1], binding_site_capacity_M)
    ]
    max_gradient_uM_per_um = [
        _max_gradient_uM_per_um(free_concentration_M, dx_um)
    ]
    penetration_depth_um_at_threshold = [
        _penetration_depth_um(
            depth_um,
            _occupancy_profile(bound_concentration_M, binding_site_capacity_M),
            occupancy_threshold_fraction,
        )
    ]

    for previous_time_s, next_time_s in zip(time_points_s, time_points_s[1:]):
        dt_s = next_time_s - previous_time_s
        free_concentration_M, bound_concentration_M = _apply_reaction_step(
            free_concentration_M,
            bound_concentration_M,
            dt_s=dt_s,
            total_loss_rate_per_s=total_loss_rate_per_s,
            effective_association_rate_M_inv_s=effective_association_rate_M_inv_s,
            effective_dissociation_rate_per_s=effective_dissociation_rate_per_s,
            binding_site_capacity_M=binding_site_capacity_M,
        )
        free_concentration_M = _apply_implicit_diffusion_step(
            free_concentration_M,
            dt_s=dt_s,
            dx_um=dx_um,
            effective_diffusion_coefficient_um2_s=effective_diffusion_coefficient_um2_s,
            source_concentration_M=source_concentration_M,
        )

        occupancy_profile = _occupancy_profile(bound_concentration_M, binding_site_capacity_M)
        mean_free_concentration_uM.append(
            sum(free_concentration_M) / grid_points / MICROMOLAR_TO_MOLAR
        )
        surface_free_concentration_uM.append(
            free_concentration_M[0] / MICROMOLAR_TO_MOLAR
        )
        distal_free_concentration_uM.append(
            free_concentration_M[-1] / MICROMOLAR_TO_MOLAR
        )
        mean_occupancy_fraction.append(sum(occupancy_profile) / grid_points)
        surface_occupancy_fraction.append(occupancy_profile[0])
        distal_occupancy_fraction.append(occupancy_profile[-1])
        max_gradient_uM_per_um.append(
            _max_gradient_uM_per_um(free_concentration_M, dx_um)
        )
        penetration_depth_um_at_threshold.append(
            _penetration_depth_um(
                depth_um,
                occupancy_profile,
                occupancy_threshold_fraction,
            )
        )

    summary = _summarize_reaction_diffusion_result(
        time_points_s=time_points_s,
        depth_um=depth_um,
        mean_free_concentration_uM=mean_free_concentration_uM,
        surface_free_concentration_uM=surface_free_concentration_uM,
        distal_free_concentration_uM=distal_free_concentration_uM,
        mean_occupancy_fraction=mean_occupancy_fraction,
        surface_occupancy_fraction=surface_occupancy_fraction,
        distal_occupancy_fraction=distal_occupancy_fraction,
        max_gradient_uM_per_um=max_gradient_uM_per_um,
        penetration_depth_um_at_threshold=penetration_depth_um_at_threshold,
        occupancy_threshold_fraction=occupancy_threshold_fraction,
        estimated_kd_uM=model_parameters["estimated_kd_uM"],
        residence_time_s=model_parameters["residence_time_s"],
    )

    return {
        "input_concentration_uM": input_concentration_uM,
        "summary": summary,
        "trajectory": {
            "time_s": time_points_s,
            "mean_free_concentration_uM": mean_free_concentration_uM,
            "surface_free_concentration_uM": surface_free_concentration_uM,
            "distal_free_concentration_uM": distal_free_concentration_uM,
            "mean_occupancy_fraction": mean_occupancy_fraction,
            "surface_occupancy_fraction": surface_occupancy_fraction,
            "distal_occupancy_fraction": distal_occupancy_fraction,
            "max_gradient_uM_per_um": max_gradient_uM_per_um,
            "penetration_depth_um_at_threshold": penetration_depth_um_at_threshold,
        },
        "spatial_final_state": {
            "depth_um": depth_um,
            "free_concentration_uM": [
                value / MICROMOLAR_TO_MOLAR for value in free_concentration_M
            ],
            "bound_concentration_uM": [
                value / MICROMOLAR_TO_MOLAR for value in bound_concentration_M
            ],
            "occupancy_fraction": _occupancy_profile(
                bound_concentration_M,
                binding_site_capacity_M,
            ),
        },
    }


def _apply_reaction_step(
    free_concentration_M: list[float],
    bound_concentration_M: list[float],
    *,
    dt_s: float,
    total_loss_rate_per_s: float,
    effective_association_rate_M_inv_s: float,
    effective_dissociation_rate_per_s: float,
    binding_site_capacity_M: float,
) -> tuple[list[float], list[float]]:
    if dt_s <= 0.0:
        return list(free_concentration_M), list(bound_concentration_M)

    loss_half_decay = math.exp(-0.5 * total_loss_rate_per_s * dt_s)
    next_free_concentration_M: list[float] = []
    next_bound_concentration_M: list[float] = []

    for free_value_M, bound_value_M in zip(free_concentration_M, bound_concentration_M):
        free_half_M = free_value_M * loss_half_decay
        association_flux_per_s = effective_association_rate_M_inv_s * free_half_M
        total_binding_rate_per_s = association_flux_per_s + effective_dissociation_rate_per_s

        if total_binding_rate_per_s <= 0.0 or binding_site_capacity_M <= 0.0:
            next_bound_M = bound_value_M
        else:
            equilibrium_bound_M = (
                binding_site_capacity_M * association_flux_per_s / total_binding_rate_per_s
            )
            next_bound_M = equilibrium_bound_M + (
                bound_value_M - equilibrium_bound_M
            ) * math.exp(-total_binding_rate_per_s * dt_s)
        next_bound_M = _clamp(next_bound_M, 0.0, binding_site_capacity_M)

        delta_bound_M = next_bound_M - bound_value_M
        free_after_binding_M = max(0.0, free_half_M - delta_bound_M)
        next_free_M = free_after_binding_M * loss_half_decay

        next_free_concentration_M.append(next_free_M)
        next_bound_concentration_M.append(next_bound_M)

    return next_free_concentration_M, next_bound_concentration_M


def _apply_implicit_diffusion_step(
    free_concentration_M: list[float],
    *,
    dt_s: float,
    dx_um: float,
    effective_diffusion_coefficient_um2_s: float,
    source_concentration_M: float,
) -> list[float]:
    next_free_concentration_M = list(free_concentration_M)
    next_free_concentration_M[0] = source_concentration_M

    if len(next_free_concentration_M) <= 2 or effective_diffusion_coefficient_um2_s <= 0.0:
        next_free_concentration_M[-1] = next_free_concentration_M[-2]
        return next_free_concentration_M

    diffusion_number = effective_diffusion_coefficient_um2_s * dt_s / (dx_um * dx_um)
    unknown_count = len(next_free_concentration_M) - 1

    lower = [0.0] * unknown_count
    diagonal = [0.0] * unknown_count
    upper = [0.0] * unknown_count
    rhs = [0.0] * unknown_count

    for index in range(unknown_count):
        original_index = index + 1
        rhs[index] = next_free_concentration_M[original_index]
        if original_index == 1:
            rhs[index] += diffusion_number * source_concentration_M

        if original_index == len(next_free_concentration_M) - 1:
            diagonal[index] = 1.0 + (2.0 * diffusion_number)
            if index > 0:
                lower[index] = -2.0 * diffusion_number
        else:
            diagonal[index] = 1.0 + (2.0 * diffusion_number)
            upper[index] = -diffusion_number
            if index > 0:
                lower[index] = -diffusion_number

    solved = _solve_tridiagonal(lower, diagonal, upper, rhs)
    updated = [source_concentration_M] + solved
    updated[-1] = updated[-2]
    return updated


def _solve_tridiagonal(
    lower: list[float],
    diagonal: list[float],
    upper: list[float],
    rhs: list[float],
) -> list[float]:
    size = len(diagonal)
    if size == 0:
        return []

    modified_upper = [0.0] * size
    modified_rhs = [0.0] * size

    modified_upper[0] = upper[0] / diagonal[0] if size > 1 else 0.0
    modified_rhs[0] = rhs[0] / diagonal[0]

    for index in range(1, size):
        denominator = diagonal[index] - (lower[index] * modified_upper[index - 1])
        modified_upper[index] = upper[index] / denominator if index < size - 1 else 0.0
        modified_rhs[index] = (
            rhs[index] - (lower[index] * modified_rhs[index - 1])
        ) / denominator

    solution = [0.0] * size
    solution[-1] = modified_rhs[-1]
    for index in range(size - 2, -1, -1):
        solution[index] = modified_rhs[index] - (
            modified_upper[index] * solution[index + 1]
        )
    return solution


def _summarize_reaction_diffusion_result(
    *,
    time_points_s: list[float],
    depth_um: list[float],
    mean_free_concentration_uM: list[float],
    surface_free_concentration_uM: list[float],
    distal_free_concentration_uM: list[float],
    mean_occupancy_fraction: list[float],
    surface_occupancy_fraction: list[float],
    distal_occupancy_fraction: list[float],
    max_gradient_uM_per_um: list[float],
    penetration_depth_um_at_threshold: list[float],
    occupancy_threshold_fraction: float,
    estimated_kd_uM: float,
    residence_time_s: float,
) -> dict[str, Any]:
    exposure_duration_s = time_points_s[-1]

    return {
        "max_mean_occupancy_fraction": max(mean_occupancy_fraction),
        "final_mean_occupancy_fraction": mean_occupancy_fraction[-1],
        "max_surface_occupancy_fraction": max(surface_occupancy_fraction),
        "final_surface_occupancy_fraction": surface_occupancy_fraction[-1],
        "max_distal_occupancy_fraction": max(distal_occupancy_fraction),
        "final_distal_occupancy_fraction": distal_occupancy_fraction[-1],
        "mean_occupancy_auc_fraction_s": _trapezoid_area(
            time_points_s,
            mean_occupancy_fraction,
        ),
        "time_above_threshold_s": _time_above_threshold(
            time_points_s,
            mean_occupancy_fraction,
            occupancy_threshold_fraction,
        ),
        "time_to_threshold_s": _find_crossing_time(
            time_points_s,
            mean_occupancy_fraction,
            occupancy_threshold_fraction,
        ),
        "max_penetration_depth_um_at_threshold": max(penetration_depth_um_at_threshold),
        "final_penetration_depth_um_at_threshold": penetration_depth_um_at_threshold[-1],
        "peak_mean_free_concentration_uM": max(mean_free_concentration_uM),
        "final_surface_free_concentration_uM": surface_free_concentration_uM[-1],
        "final_distal_free_concentration_uM": distal_free_concentration_uM[-1],
        "max_gradient_uM_per_um": max(max_gradient_uM_per_um),
        "domain_length_um": depth_um[-1],
        "estimated_kd_uM": estimated_kd_uM,
        "residence_time_s": residence_time_s,
        "selection_peak_occupancy_fraction": max(mean_occupancy_fraction),
        "mean_occupancy_fraction": (
            _trapezoid_area(time_points_s, mean_occupancy_fraction) / exposure_duration_s
            if exposure_duration_s > 0.0
            else 0.0
        ),
    }


def _occupancy_fraction(bound_value_M: float, binding_site_capacity_M: float) -> float:
    if binding_site_capacity_M <= 0.0:
        return 0.0
    return _clamp(bound_value_M / binding_site_capacity_M, 0.0, 1.0)


def _mean_occupancy_fraction(
    bound_concentration_M: list[float],
    binding_site_capacity_M: float,
) -> float:
    return sum(_occupancy_profile(bound_concentration_M, binding_site_capacity_M)) / len(
        bound_concentration_M
    )


def _occupancy_profile(
    bound_concentration_M: list[float],
    binding_site_capacity_M: float,
) -> list[float]:
    return [
        _occupancy_fraction(value_M, binding_site_capacity_M)
        for value_M in bound_concentration_M
    ]


def _penetration_depth_um(
    depth_um: list[float],
    occupancy_fraction: list[float],
    threshold: float,
) -> float:
    if threshold <= 0.0:
        return depth_um[-1]
    if not occupancy_fraction or occupancy_fraction[0] < threshold:
        return 0.0

    penetration_depth = 0.0
    for left_depth_um, right_depth_um, left_value, right_value in zip(
        depth_um,
        depth_um[1:],
        occupancy_fraction,
        occupancy_fraction[1:],
    ):
        if right_value >= threshold:
            penetration_depth = right_depth_um
            continue
        if left_value >= threshold > right_value and right_value != left_value:
            fraction = (threshold - left_value) / (right_value - left_value)
            fraction = _clamp(fraction, 0.0, 1.0)
            return left_depth_um + (fraction * (right_depth_um - left_depth_um))
    return penetration_depth


def _max_gradient_uM_per_um(
    free_concentration_M: list[float],
    dx_um: float,
) -> float:
    if dx_um <= 0.0 or len(free_concentration_M) < 2:
        return 0.0
    return max(
        abs(right_value - left_value) / MICROMOLAR_TO_MOLAR / dx_um
        for left_value, right_value in zip(
            free_concentration_M,
            free_concentration_M[1:],
        )
    )


def _select_reaction_diffusion_operating_point(
    concentration_results: list[dict[str, Any]],
    *,
    occupancy_threshold_fraction: float,
) -> dict[str, Any] | None:
    threshold_hits = [
        result
        for result in concentration_results
        if result["summary"]["max_mean_occupancy_fraction"] >= occupancy_threshold_fraction
    ]
    if threshold_hits:
        threshold_hits.sort(
            key=lambda result: (
                result["input_concentration_uM"],
                -result["summary"]["max_distal_occupancy_fraction"],
            )
        )
        return threshold_hits[0]

    if not concentration_results:
        return None

    return max(
        concentration_results,
        key=lambda result: (
            result["summary"]["max_mean_occupancy_fraction"],
            result["summary"]["max_distal_occupancy_fraction"],
            result["summary"]["max_penetration_depth_um_at_threshold"],
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
        if result["summary"]["max_mean_occupancy_fraction"] >= occupancy_threshold_fraction
    ]
    return min(matches) if matches else None
