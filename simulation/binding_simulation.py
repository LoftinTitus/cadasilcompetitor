"""Binding simulation entry points for Tier A and Tier B occupancy modeling."""

from __future__ import annotations

from typing import Any

from core.config_loader import (
    CandidateProperties,
    HSVariantConfig,
    HSVariantPanelConfig,
    PhysiologyConfig,
    SimulationConfig,
)
from simulation.candidate_property_estimation import estimate_candidate_properties
from simulation.occupancy_model import (
    prepare_occupancy_model_inputs,
    simulate_occupancy_sweep,
)
from simulation.reaction_diffusion_model import simulate_reaction_diffusion_sweep


def run_binding_simulation(
    candidate: Any,
    simulation_config: SimulationConfig | PhysiologyConfig,
    candidate_properties: CandidateProperties | None = None,
    hs_variant: HSVariantConfig | None = None,
    model_tier: str = "A",
) -> dict[str, Any]:
    resolved_candidate_properties = (
        candidate_properties
        if candidate_properties is not None
        else estimate_candidate_properties(candidate)
    )

    model_inputs = prepare_occupancy_model_inputs(
        candidate,
        simulation_config,
        resolved_candidate_properties,
        hs_variant=hs_variant,
    )
    normalized_model_tier = model_tier.strip().upper()
    if normalized_model_tier == "A":
        occupancy_results = simulate_occupancy_sweep(model_inputs)
    elif normalized_model_tier == "B":
        occupancy_results = simulate_reaction_diffusion_sweep(model_inputs)
    else:
        raise ValueError("model_tier must be 'A' or 'B'.")

    return {
        **model_inputs,
        **occupancy_results,
    }


def run_binding_simulation_panel(
    candidate: Any,
    simulation_config: SimulationConfig | PhysiologyConfig,
    candidate_properties: CandidateProperties | None = None,
    hs_variant_panel: HSVariantPanelConfig | list[HSVariantConfig] | tuple[HSVariantConfig, ...] | None = None,
    model_tier: str = "A",
) -> dict[str, Any]:
    if hs_variant_panel is None:
        raise ValueError("hs_variant_panel is required for panel simulation.")

    variants = (
        hs_variant_panel.variants
        if isinstance(hs_variant_panel, HSVariantPanelConfig)
        else tuple(hs_variant_panel)
    )

    per_variant_results = [
        run_binding_simulation(
            candidate,
            simulation_config,
            candidate_properties,
            hs_variant=variant,
            model_tier=model_tier,
        )
        for variant in variants
    ]

    return {
        "candidate": per_variant_results[0]["candidate"] if per_variant_results else None,
        "panel_id": (
            hs_variant_panel.panel_id if isinstance(hs_variant_panel, HSVariantPanelConfig) else None
        ),
        "panel_name": (
            hs_variant_panel.panel_name
            if isinstance(hs_variant_panel, HSVariantPanelConfig)
            else None
        ),
        "results": per_variant_results,
        "summary": _summarize_panel_results(per_variant_results),
    }


def run_binding_simulation_tier_b(
    candidate: Any,
    simulation_config: SimulationConfig | PhysiologyConfig,
    candidate_properties: CandidateProperties | None = None,
    hs_variant: HSVariantConfig | None = None,
) -> dict[str, Any]:
    return run_binding_simulation(
        candidate,
        simulation_config,
        candidate_properties,
        hs_variant=hs_variant,
        model_tier="B",
    )


def _summarize_panel_results(per_variant_results: list[dict[str, Any]]) -> dict[str, Any]:
    if not per_variant_results:
        return {
            "best_variant_id": None,
            "worst_variant_id": None,
            "highest_max_occupancy_fraction": None,
            "lowest_estimated_kd_uM": None,
        }

    ranked_results = sorted(
        per_variant_results,
        key=lambda result: (
            result["summary"].get(
                "selection_peak_occupancy_fraction",
                result["summary"].get("best_case_max_occupancy_fraction", 0.0),
            ),
            -result["model_parameters"]["estimated_kd_uM"],
        ),
        reverse=True,
    )
    best_result = ranked_results[0]
    worst_result = ranked_results[-1]

    return {
        "best_variant_id": _variant_id(best_result),
        "worst_variant_id": _variant_id(worst_result),
        "highest_max_occupancy_fraction": best_result["summary"][
            "selection_peak_occupancy_fraction"
        ],
        "lowest_estimated_kd_uM": min(
            result["model_parameters"]["estimated_kd_uM"] for result in per_variant_results
        ),
    }


def _variant_id(result: dict[str, Any]) -> str | None:
    hs_variant = result.get("hs_variant")
    if hs_variant is None:
        return None
    return hs_variant.get("variant_id")
