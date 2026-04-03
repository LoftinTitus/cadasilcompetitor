from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from reports.run_report import main as run_report_main


class ReportGenerationTests(unittest.TestCase):
    def test_run_report_writes_expected_artifacts(self) -> None:
        ranked_rows = [
            {
                "rank": 1,
                "candidate_id": "pep_00001",
                "sequence": "AKRKRQGK",
                "screening_status": "pass",
                "composite_screen_score": 58.2,
                "hs_affinity_reward": 0.81,
                "hs_selectivity_reward": 0.31,
                "transport_reward": 0.56,
                "residence_reward": 0.43,
                "anticoagulant_penalty": 0.08,
                "off_target_penalty": 0.09,
                "heparin_penalty": 0.04,
                "developability_penalty": 0.05,
                "tier_a_best_variant_id": "HS-dp4-NS-6S",
                "transport_probe_variant_id": "HS-dp6-mixed",
                "warning_flags": "",
                "filter_flags": "",
            },
            {
                "rank": 2,
                "candidate_id": "pep_00002",
                "sequence": "AKRGRKRRKQGA",
                "screening_status": "review",
                "composite_screen_score": 47.4,
                "hs_affinity_reward": 0.74,
                "hs_selectivity_reward": 0.22,
                "transport_reward": 0.51,
                "residence_reward": 0.39,
                "anticoagulant_penalty": 0.14,
                "off_target_penalty": 0.12,
                "heparin_penalty": 0.06,
                "developability_penalty": 0.07,
                "tier_a_best_variant_id": "HS-dp5-AT-motif",
                "transport_probe_variant_id": "HS-dp6-mixed",
                "warning_flags": "aggregation_risk",
                "filter_flags": "",
            },
            {
                "rank": 3,
                "candidate_id": "pep_00003",
                "sequence": "GKKQGKQEGKKQGKQ",
                "screening_status": "reject",
                "composite_screen_score": 26.8,
                "hs_affinity_reward": 0.44,
                "hs_selectivity_reward": 0.02,
                "transport_reward": 0.34,
                "residence_reward": 0.28,
                "anticoagulant_penalty": 0.34,
                "off_target_penalty": 0.22,
                "heparin_penalty": 0.18,
                "developability_penalty": 0.10,
                "tier_a_best_variant_id": "HS-dp4-NAc",
                "transport_probe_variant_id": "HS-dp6-mixed",
                "warning_flags": "",
                "filter_flags": "basic_fraction_out_of_range",
            },
        ]
        proposal_rows = [
            {
                "proposal_rank": 1,
                "candidate_id": "pep_next",
                "sequence": "GKRRKAKRGRR",
                "warning_flags": "",
                "filter_flags": "",
                "acquisition_score": 0.73,
                "predicted_composite_score": 55.1,
                "novelty_score": 0.42,
            }
        ]

        with tempfile.TemporaryDirectory() as temp_dir:
            ranked_path = Path(temp_dir) / "ranked.csv"
            proposals_path = Path(temp_dir) / "proposals.csv"
            output_dir = Path(temp_dir) / "reports"
            self._write_csv(ranked_path, ranked_rows)
            self._write_csv(proposals_path, proposal_rows)

            with patch(
                "sys.argv",
                [
                    "run_report.py",
                    "--ranked",
                    str(ranked_path),
                    "--proposals",
                    str(proposals_path),
                    "--output-dir",
                    str(output_dir),
                ],
            ):
                exit_code = run_report_main()

            self.assertEqual(exit_code, 0)
            self.assertTrue((output_dir / "screening_summary.md").exists())
            self.assertTrue((output_dir / "manifest_validation_summary.md").exists())
            self.assertTrue((output_dir / "pareto_frontier.csv").exists())
            self.assertTrue((output_dir / "pareto_affinity_selectivity.svg").exists())
            self.assertTrue((output_dir / "pareto_transport_safety.svg").exists())

            summary_text = (output_dir / "screening_summary.md").read_text(encoding="utf-8")
            self.assertIn("# Screening Summary", summary_text)
            self.assertIn("Manifest Readiness", summary_text)
            self.assertIn("Top Ranked Candidates", summary_text)
            self.assertIn("Manifest Coverage", summary_text)
            self.assertIn("Top Surrogate Proposals", summary_text)

            manifest_text = (output_dir / "manifest_validation_summary.md").read_text(
                encoding="utf-8"
            )
            self.assertIn("# Manifest Validation Summary", manifest_text)
            self.assertIn("HS Panel Alignment", manifest_text)

    def _write_csv(self, path: Path, rows: list[dict[str, object]]) -> None:
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    unittest.main()
