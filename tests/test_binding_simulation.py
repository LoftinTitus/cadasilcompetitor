from __future__ import annotations

import unittest
from dataclasses import replace

from core.config_loader import CandidateProperties, TargetConfig, load_hs_variant_panel, load_simulation_config
from simulation.binding_simulation import run_binding_simulation, run_binding_simulation_panel


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

    def _variant(self, variant_id: str):
        for variant in self.hs_variant_panel.variants:
            if variant.variant_id == variant_id:
                return variant
        raise AssertionError(f"Could not find variant {variant_id!r}")


if __name__ == "__main__":
    unittest.main()
