from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from core.config_loader import load_hs_variant_panel, load_simulation_config
from scoring.screening import (
    load_candidates_from_csv,
    rank_candidates,
    screen_candidate,
    screen_candidates_from_csv,
    write_ranked_candidates_csv,
)


class ScreeningTests(unittest.TestCase):
    def setUp(self) -> None:
        self.simulation_config = load_simulation_config("config/physiology.yaml")
        self.hs_variant_panel = load_hs_variant_panel("config/hs_variant_panel.yaml")
        self.safe_candidate = {
            "candidate_id": "pep_safe",
            "sequence": "AKRKRQGK",
            "length": 8,
            "approx_net_charge": 5,
            "basic_fraction": 0.625,
            "hydrophobic_fraction": 0.125,
            "acidic_fraction": 0.0,
            "polar_fraction": 0.25,
            "longest_basic_run": 2,
            "longest_hydrophobic_run": 1,
            "warning_flags": [],
            "filter_flags": [],
        }
        self.risky_candidate = {
            "candidate_id": "pep_risky",
            "sequence": "RRRRKKKKWW",
            "length": 10,
            "approx_net_charge": 8,
            "basic_fraction": 0.8,
            "hydrophobic_fraction": 0.2,
            "acidic_fraction": 0.0,
            "polar_fraction": 0.0,
            "longest_basic_run": 8,
            "longest_hydrophobic_run": 2,
            "warning_flags": ["cpp_like_uptake_risk", "aggregation_risk"],
            "filter_flags": ["basic_fraction_out_of_range"],
        }

    def test_screen_candidate_returns_scoring_columns(self) -> None:
        row = screen_candidate(
            self.safe_candidate,
            simulation_config=self.simulation_config,
            hs_variant_panel=self.hs_variant_panel,
        )

        self.assertIn("composite_screen_score", row)
        self.assertIn("tier_a_best_variant_id", row)
        self.assertIn("tier_b_mixed_peak_mean_occupancy", row)
        self.assertIn("estimated_association_rate_M_inv_s", row)
        self.assertEqual(row["candidate_id"], "pep_safe")

    def test_rank_candidates_orders_candidates_by_score(self) -> None:
        ranked_rows = rank_candidates(
            [self.risky_candidate, self.safe_candidate],
            simulation_config=self.simulation_config,
            hs_variant_panel=self.hs_variant_panel,
        )

        self.assertEqual(ranked_rows[0]["rank"], 1)
        self.assertEqual(ranked_rows[1]["rank"], 2)
        self.assertGreaterEqual(
            float(ranked_rows[0]["composite_screen_score"]),
            float(ranked_rows[1]["composite_screen_score"]),
        )
        self.assertEqual(ranked_rows[0]["candidate_id"], "pep_safe")
        self.assertEqual(ranked_rows[1]["screening_status"], "reject")

    def test_write_ranked_candidates_csv_writes_expected_header(self) -> None:
        ranked_rows = rank_candidates(
            [self.safe_candidate],
            simulation_config=self.simulation_config,
            hs_variant_panel=self.hs_variant_panel,
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            output_path = Path(temp_dir) / "ranked.csv"
            write_ranked_candidates_csv(output_path, ranked_rows)

            with output_path.open(newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle)
                rows = list(reader)

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["candidate_id"], "pep_safe")
        self.assertIn("composite_screen_score", rows[0])

    def test_screen_candidates_from_csv_runs_end_to_end(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            input_path = Path(temp_dir) / "candidates.csv"
            output_path = Path(temp_dir) / "ranked.csv"
            with input_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "candidate_id",
                        "sequence",
                        "length",
                        "approx_net_charge",
                        "basic_fraction",
                        "hydrophobic_fraction",
                        "acidic_fraction",
                        "polar_fraction",
                        "longest_basic_run",
                        "longest_hydrophobic_run",
                        "warning_flags",
                        "filter_flags",
                    ],
                )
                writer.writeheader()
                writer.writerow(self.safe_candidate)
                writer.writerow(
                    {
                        **self.risky_candidate,
                        "warning_flags": ";".join(self.risky_candidate["warning_flags"]),
                        "filter_flags": ";".join(self.risky_candidate["filter_flags"]),
                    }
                )

            ranked_rows = screen_candidates_from_csv(
                input_path,
                output_csv_path=output_path,
                simulation_config=self.simulation_config,
                hs_variant_panel=self.hs_variant_panel,
            )
            loaded_rows = load_candidates_from_csv(input_path)
            output_exists = output_path.exists()

        self.assertEqual(len(loaded_rows), 2)
        self.assertEqual(len(ranked_rows), 2)
        self.assertTrue(output_exists)


if __name__ == "__main__":
    unittest.main()
