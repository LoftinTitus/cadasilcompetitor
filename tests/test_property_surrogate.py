from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from models.property_surrogate import (
    fit_property_surrogate,
    load_property_training_rows,
    try_load_property_surrogate,
)
from simulation.candidate_property_estimation import estimate_candidate_properties


def _training_row(
    *,
    candidate_id: str,
    sequence: str,
    estimated_diffusion_coefficient_um2_s: float,
    estimated_clearance_rate_per_s: float,
    estimated_association_rate_M_inv_s: float,
    estimated_dissociation_rate_per_s: float,
    estimated_barrier_permeability_cm_s: float,
    estimated_half_life_s: float,
    estimated_enzymatic_degradation_rate_per_s: float,
    estimated_spontaneous_degradation_rate_per_s: float,
) -> dict[str, object]:
    return {
        "candidate_id": candidate_id,
        "sequence": sequence,
        "warning_flags": [],
        "filter_flags": [],
        "estimated_diffusion_coefficient_um2_s": estimated_diffusion_coefficient_um2_s,
        "estimated_clearance_rate_per_s": estimated_clearance_rate_per_s,
        "estimated_association_rate_M_inv_s": estimated_association_rate_M_inv_s,
        "estimated_dissociation_rate_per_s": estimated_dissociation_rate_per_s,
        "estimated_barrier_permeability_cm_s": estimated_barrier_permeability_cm_s,
        "estimated_half_life_s": estimated_half_life_s,
        "estimated_enzymatic_degradation_rate_per_s": estimated_enzymatic_degradation_rate_per_s,
        "estimated_spontaneous_degradation_rate_per_s": estimated_spontaneous_degradation_rate_per_s,
    }


class PropertySurrogateTests(unittest.TestCase):
    def setUp(self) -> None:
        self.training_rows = [
            _training_row(
                candidate_id="row_1",
                sequence="AKRKRQGK",
                estimated_diffusion_coefficient_um2_s=111.0,
                estimated_clearance_rate_per_s=0.00072,
                estimated_association_rate_M_inv_s=2.1e5,
                estimated_dissociation_rate_per_s=0.24,
                estimated_barrier_permeability_cm_s=3.8e-7,
                estimated_half_life_s=5200.0,
                estimated_enzymatic_degradation_rate_per_s=7.2e-5,
                estimated_spontaneous_degradation_rate_per_s=1.1e-5,
            ),
            _training_row(
                candidate_id="row_2",
                sequence="AKRGRKRRKQGA",
                estimated_diffusion_coefficient_um2_s=104.0,
                estimated_clearance_rate_per_s=0.00081,
                estimated_association_rate_M_inv_s=2.3e5,
                estimated_dissociation_rate_per_s=0.21,
                estimated_barrier_permeability_cm_s=4.0e-7,
                estimated_half_life_s=6100.0,
                estimated_enzymatic_degradation_rate_per_s=8.4e-5,
                estimated_spontaneous_degradation_rate_per_s=1.2e-5,
            ),
            _training_row(
                candidate_id="row_3",
                sequence="GKKQGKQEGKKQGKQ",
                estimated_diffusion_coefficient_um2_s=96.0,
                estimated_clearance_rate_per_s=0.00089,
                estimated_association_rate_M_inv_s=1.9e5,
                estimated_dissociation_rate_per_s=0.28,
                estimated_barrier_permeability_cm_s=3.3e-7,
                estimated_half_life_s=7200.0,
                estimated_enzymatic_degradation_rate_per_s=9.0e-5,
                estimated_spontaneous_degradation_rate_per_s=1.4e-5,
            ),
        ]

    def test_property_surrogate_predicts_candidate_properties(self) -> None:
        model = fit_property_surrogate(
            self.training_rows,
            ensemble_size=5,
            random_seed=11,
        )

        prediction = model.predict_candidate_properties({"sequence": "AKRGRQKRRKA"})

        self.assertGreater(prediction.diffusion_coefficient_um2_s, 0.0)
        self.assertGreater(prediction.association_rate_M_inv_s, 0.0)
        self.assertGreater(prediction.half_life_s, 0.0)

    def test_property_surrogate_round_trips_through_json(self) -> None:
        model = fit_property_surrogate(
            self.training_rows,
            ensemble_size=4,
            random_seed=3,
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            output_path = Path(temp_dir) / "property_model.json"
            model.write_json(output_path)
            loaded = try_load_property_surrogate(output_path)

        self.assertIsNotNone(loaded)
        loaded_prediction = loaded.predict({"sequence": "AKRGRQKRRKA"})
        self.assertIn("estimated_half_life_s", loaded_prediction)
        self.assertIn("mean", loaded_prediction["estimated_half_life_s"])

    def test_load_property_training_rows_parses_csv_scalars(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            csv_path = Path(temp_dir) / "training.csv"
            csv_path.write_text(
                "\n".join(
                    [
                        ",".join(self.training_rows[0].keys()),
                        ",".join(str(value) for value in self.training_rows[0].values()),
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            rows = load_property_training_rows(csv_path)

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["candidate_id"], "row_1")
        self.assertEqual(rows[0]["warning_flags"], [])
        self.assertAlmostEqual(rows[0]["estimated_half_life_s"], 5200.0)

    def test_estimate_candidate_properties_can_use_learned_estimator(self) -> None:
        heuristic = estimate_candidate_properties(
            {"sequence": "AKRGRQKRRKA"},
            property_estimator_path=None,
        )

        class StubEstimator:
            def predict_candidate_properties(self, candidate):
                self.last_candidate = candidate
                return heuristic.__class__(
                    diffusion_coefficient_um2_s=80.0,
                    clearance_rate_per_s=0.0012,
                    association_rate_M_inv_s=300000.0,
                    dissociation_rate_per_s=0.18,
                    barrier_permeability_cm_s=5.0e-7,
                    half_life_s=8000.0,
                    enzymatic_degradation_rate_per_s=6.0e-5,
                    spontaneous_degradation_rate_per_s=1.6e-5,
                )

        learned = estimate_candidate_properties(
            {"sequence": "AKRGRQKRRKA"},
            learned_estimator=StubEstimator(),
            property_estimator_path=None,
        )

        self.assertNotEqual(
            learned.diffusion_coefficient_um2_s,
            heuristic.diffusion_coefficient_um2_s,
        )
        self.assertGreater(learned.half_life_s, heuristic.half_life_s)


if __name__ == "__main__":
    unittest.main()
