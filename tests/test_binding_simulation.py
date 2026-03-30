from __future__ import annotations

import unittest
from dataclasses import replace

from core.config_loader import CandidateProperties, TargetConfig, load_hs_variant_panel, load_simulation_config
from simulation.binding_simulation import (
    run_binding_simulation,
    run_binding_simulation_panel,
    run_binding_simulation_tier_b,
)
from simulation.candidate_property_estimation import (
    estimate_candidate_properties,
    estimate_candidate_properties_with_diagnostics,
)


class BindingSimulationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.simulation_config = load_simulation_config("config/physiology.yaml")
        self.hs_variant_panel = load_hs_variant_panel("config/hs_variant_panel.yaml")
        self.candidate = {
            "candidate_id": "pep_test",
            "sequence": "AKRKRQGK",
            "approx_net_charge": 5,
        }
        self.candidate_properties = CandidateProperties(
            diffusion_coefficient_um2_s=80.0,
            clearance_rate_per_s=0.001,
            association_rate_M_inv_s=2.0e5,
            dissociation_rate_per_s=0.2,
            barrier_permeability_cm_s=1.0e-6,
            half_life_s=7200.0,
            enzymatic_degradation_rate_per_s=5.0e-5,
            spontaneous_degradation_rate_per_s=2.0e-5,
        )

    def test_config_includes_simulation_controls(self) -> None:
        controls = self.simulation_config.simulation_controls
        self.assertEqual(tuple(controls.concentration_grid_uM), (0.01, 0.1, 1.0, 10.0))
        self.assertEqual(controls.occupancy_threshold_fraction, 0.5)
        self.assertEqual(controls.time_step_s, 5.0)
        self.assertEqual(controls.spatial_grid_points, 25)

    def test_single_variant_run_returns_occupancy_sweep(self) -> None:
        variant = self._variant("HS-dp4-NS-6S")
        result = run_binding_simulation(
            self.candidate,
            self.simulation_config,
            self.candidate_properties,
            hs_variant=variant,
        )

        self.assertEqual(len(result["concentration_sweep"]), 4)
        self.assertIsNotNone(result["operating_point"])
        self.assertGreater(
            result["summary"]["best_case_max_occupancy_fraction"],
            0.0,
        )
        self.assertIn("effective_association_rate_M_inv_s", result["model_parameters"])

    def test_candidate_property_estimator_returns_valid_properties(self) -> None:
        properties = estimate_candidate_properties(self.candidate)

        self.assertGreater(properties.diffusion_coefficient_um2_s, 0.0)
        self.assertGreater(properties.association_rate_M_inv_s, 0.0)
        self.assertGreater(properties.half_life_s, 0.0)
        self.assertGreaterEqual(properties.clearance_rate_per_s, 0.0)

    def test_candidate_property_estimator_uses_warning_flags(self) -> None:
        calm_candidate = {
            "candidate_id": "pep_calm",
            "sequence": "AKRKRQGK",
            "warning_flags": [],
        }
        risky_candidate = {
            "candidate_id": "pep_risky",
            "sequence": "RRRRKKKKWW",
            "warning_flags": ["cpp_like_uptake_risk", "aggregation_risk"],
        }

        calm_properties = estimate_candidate_properties(calm_candidate)
        risky_properties = estimate_candidate_properties(risky_candidate)

        self.assertGreater(
            risky_properties.clearance_rate_per_s,
            calm_properties.clearance_rate_per_s,
        )
        self.assertGreater(
            risky_properties.barrier_permeability_cm_s,
            calm_properties.barrier_permeability_cm_s,
        )

    def test_estimator_diagnostics_include_metadata_and_flat_dict(self) -> None:
        diagnostics = estimate_candidate_properties_with_diagnostics(self.candidate)

        self.assertIn("candidate_properties", diagnostics)
        self.assertIn("candidate_properties_dict", diagnostics)
        self.assertIn("sequence_metadata", diagnostics)
        self.assertEqual(diagnostics["sequence_metadata"]["sequence"], "AKRKRQGK")

    def test_glycocalyx_is_more_accessible_than_basement_membrane_for_same_candidate(self) -> None:
        glycocalyx_result = run_binding_simulation(
            self.candidate,
            self.simulation_config,
            self.candidate_properties,
        )
        basement_target = TargetConfig(
            compartment="basement_membrane",
            exposure_side="luminal",
            biological_rationale_label="subendothelial_retention",
        )
        basement_config = replace(self.simulation_config, target=basement_target)
        basement_result = run_binding_simulation(
            self.candidate,
            basement_config,
            self.candidate_properties,
        )

        self.assertGreater(
            glycocalyx_result["model_parameters"]["delivery_rate_per_s"],
            basement_result["model_parameters"]["delivery_rate_per_s"],
        )
        self.assertGreater(
            glycocalyx_result["summary"]["best_case_max_occupancy_fraction"],
            basement_result["summary"]["best_case_max_occupancy_fraction"],
        )

    def test_high_affinity_sentinel_variant_outperforms_unsulfated_control(self) -> None:
        unsulfated = self._variant("HS-dp4-NAc")
        sentinel = self._variant("HS-dp5-AT-motif")

        unsulfated_result = run_binding_simulation(
            self.candidate,
            self.simulation_config,
            self.candidate_properties,
            hs_variant=unsulfated,
        )
        sentinel_result = run_binding_simulation(
            self.candidate,
            self.simulation_config,
            self.candidate_properties,
            hs_variant=sentinel,
        )

        self.assertLess(
            sentinel_result["model_parameters"]["estimated_kd_uM"],
            unsulfated_result["model_parameters"]["estimated_kd_uM"],
        )
        self.assertGreater(
            sentinel_result["summary"]["best_case_max_occupancy_fraction"],
            unsulfated_result["summary"]["best_case_max_occupancy_fraction"],
        )

    def test_panel_run_summarizes_variant_results(self) -> None:
        result = run_binding_simulation_panel(
            self.candidate,
            self.simulation_config,
            self.candidate_properties,
            self.hs_variant_panel,
        )

        self.assertEqual(len(result["results"]), len(self.hs_variant_panel.variants))
        self.assertEqual(result["summary"]["best_variant_id"], "HS-dp5-AT-motif")

    def test_simulation_can_estimate_candidate_properties_automatically(self) -> None:
        result = run_binding_simulation(
            self.candidate,
            self.simulation_config,
            hs_variant=self._variant("HS-dp4-NS-6S"),
        )

        self.assertIsNotNone(result["candidate_properties"])
        self.assertGreater(
            result["candidate_properties"]["association_rate_M_inv_s"],
            0.0,
        )
        self.assertGreater(
            result["summary"]["best_case_max_occupancy_fraction"],
            0.0,
        )

    def test_panel_run_can_estimate_candidate_properties_automatically(self) -> None:
        result = run_binding_simulation_panel(
            self.candidate,
            self.simulation_config,
            hs_variant_panel=self.hs_variant_panel,
        )

        self.assertEqual(len(result["results"]), len(self.hs_variant_panel.variants))
        self.assertEqual(result["results"][0]["candidate"]["candidate_id"], "pep_test")

    def test_tier_b_run_returns_spatial_profile(self) -> None:
        variant = self._variant("HS-dp4-NS-6S")
        result = run_binding_simulation_tier_b(
            self.candidate,
            self.simulation_config,
            self.candidate_properties,
            hs_variant=variant,
        )

        self.assertEqual(result["summary"]["model_tier"], "B")
        self.assertEqual(len(result["concentration_sweep"]), 4)
        first_result = result["concentration_sweep"][0]
        self.assertIn("spatial_final_state", first_result)
        self.assertEqual(
            len(first_result["spatial_final_state"]["depth_um"]),
            self.simulation_config.simulation_controls.spatial_grid_points,
        )

    def test_tier_b_shows_spatial_gradient_and_surface_first_behavior(self) -> None:
        result = run_binding_simulation(
            self.candidate,
            self.simulation_config,
            self.candidate_properties,
            model_tier="B",
        )

        selected = result["operating_point"]
        self.assertIsNotNone(selected)
        self.assertGreater(
            selected["summary"]["max_gradient_uM_per_um"],
            0.0,
        )
        self.assertGreater(
            selected["summary"]["max_surface_occupancy_fraction"],
            selected["summary"]["max_distal_occupancy_fraction"],
        )
        self.assertGreaterEqual(
            selected["summary"]["max_penetration_depth_um_at_threshold"],
            0.0,
        )

    def test_tier_b_can_estimate_candidate_properties_automatically(self) -> None:
        result = run_binding_simulation_tier_b(
            self.candidate,
            self.simulation_config,
            hs_variant=self._variant("HS-dp4-NS-6S"),
        )

        self.assertEqual(result["summary"]["model_tier"], "B")
        self.assertIsNotNone(result["candidate_properties"])

    def _variant(self, variant_id: str):
        for variant in self.hs_variant_panel.variants:
            if variant.variant_id == variant_id:
                return variant
        raise AssertionError(f"Could not find variant {variant_id!r}")


if __name__ == "__main__":
    unittest.main()
