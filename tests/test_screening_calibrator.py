from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from models.features import extract_feature_map
from models.run_screening_calibration import main as run_screening_calibration_main
from models.screening_calibrator import (
    fit_screening_calibrator,
    load_screening_training_rows,
    try_load_screening_calibrator,
)


class ScreeningCalibratorTests(unittest.TestCase):
    def setUp(self) -> None:
        self.training_rows = [
            {
                "candidate_id": "pass_1",
                "sequence": "AKRKRQGK",
                "screening_status": "pass",
                "composite_screen_score": 55.0,
            },
            {
                "candidate_id": "review_1",
                "sequence": "GKKQGKQEGKKQGKQ",
                "screening_status": "review",
                "composite_screen_score": 42.0,
            },
            {
                "candidate_id": "reject_1",
                "sequence": "RRRRKKKKWW",
                "screening_status": "reject",
                "composite_screen_score": 20.0,
            },
        ]

    def test_feature_map_includes_window_and_motif_features(self) -> None:
        features = extract_feature_map({"sequence": "AKRKRQGK"})

        self.assertIn("max_basic_window_fraction_5", features)
        self.assertIn("motif_KR", features)
        self.assertIn("motif_BBXB", features)
        self.assertGreater(features["motif_KR"], 0.0)

    def test_fit_screening_calibrator_predicts_summary(self) -> None:
        calibrator = fit_screening_calibrator(
            self.training_rows,
            ensemble_size=5,
            random_seed=11,
        )

        summary = calibrator.predict_summary({"sequence": "AKRGRQKRRKA"})

        self.assertGreaterEqual(summary["ml_screening_score"], 0.0)
        self.assertLessEqual(summary["ml_screening_score"], 1.0)
        self.assertGreaterEqual(summary["ml_predicted_composite_score"], 0.0)
        self.assertLessEqual(summary["ml_predicted_composite_score"], 100.0)

    def test_screening_calibrator_round_trips_json(self) -> None:
        calibrator = fit_screening_calibrator(
            self.training_rows,
            ensemble_size=4,
            random_seed=19,
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            output_path = Path(temp_dir) / "screening_calibrator.json"
            calibrator.write_json(output_path)
            loaded = try_load_screening_calibrator(output_path)

        self.assertIsNotNone(loaded)
        self.assertIn("ml_screening_score", loaded.predict_summary({"sequence": "AKRGRQKRRKA"}))

    def test_load_screening_training_rows_derives_targets(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            training_path = Path(temp_dir) / "ranked.csv"
            training_path.write_text(
                "\n".join(
                    [
                        "candidate_id,sequence,screening_status,composite_screen_score",
                        "pass_1,AKRKRQGK,pass,55.0",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            rows = load_screening_training_rows(training_path)

        self.assertEqual(rows[0]["screening_status_score"], 1.0)
        self.assertEqual(rows[0]["composite_screen_score_scaled"], 0.55)

    def test_run_screening_calibration_writes_model(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            training_path = Path(temp_dir) / "ranked.csv"
            output_path = Path(temp_dir) / "screening_calibrator.json"
            training_path.write_text(
                "\n".join(
                    [
                        "candidate_id,sequence,screening_status,composite_screen_score",
                        "pass_1,AKRKRQGK,pass,55.0",
                        "review_1,GKKQGKQEGKKQGKQ,review,42.0",
                        "reject_1,RRRRKKKKWW,reject,20.0",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            with patch(
                "sys.argv",
                [
                    "run_screening_calibration.py",
                    "--training",
                    str(training_path),
                    "--output",
                    str(output_path),
                    "--ensemble-size",
                    "3",
                ],
            ):
                exit_code = run_screening_calibration_main()

            self.assertEqual(exit_code, 0)
            self.assertTrue(output_path.exists())


if __name__ == "__main__":
    unittest.main()
