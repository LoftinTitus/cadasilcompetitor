"""Binding simulation entry points for Tier A occupancy modeling.

The implementation here follows the README's first transport milestone:
a compartment-style model that maps candidate kinetics and environment
configuration into time-dependent HS occupancy.
"""

from __future__ import annotations

from typing import Any

from core.config_loader import (
    CandidateProperties,
    HSVariantConfig,
    HSVariantPanelConfig,
    PhysiologyConfig,
    SimulationConfig,
)
from simulation.occupancy_model import (
    prepare_occupancy_model_inputs,
    simulate_occupancy_sweep,
)


def run_binding_simulation(
    candidate: Any,
    simulation_config: SimulationConfig | PhysiologyConfig,
    candidate_properties: CandidateProperties | None = None,
    hs_variant: HSVariantConfig | None = None,
) -> dict[str, Any]:
    if candidate_properties is None:
        raise ValueError(
            "candidate_properties are required for occupancy simulation because the "
            "README model depends on diffusion, clearance, kon, koff, and degradation inputs."
        )

    model_inputs = prepare_occupancy_model_inputs(
        candidate,
        simulation_config,
        candidate_properties,
        hs_variant=hs_variant,
    )
    occupancy_results = simulate_occupancy_sweep(model_inputs)

    return {
        **model_inputs,
        **occupancy_results,
    }


def run_binding_simulation_panel(
    candidate: Any,
    simulation_config: SimulationConfig | PhysiologyConfig,
    candidate_properties: CandidateProperties,
    hs_variant_panel: HSVariantPanelConfig | list[HSVariantConfig] | tuple[HSVariantConfig, ...],
) -> dict[str, Any]:
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
            result["summary"]["best_case_max_occupancy_fraction"],
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
            "best_case_max_occupancy_fraction"
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
