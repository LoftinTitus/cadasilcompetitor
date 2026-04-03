from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from core.config_loader import load_hs_variant_panel
from core.manifest_loader import load_data_manifest
from core.manifest_registry import annotate_hs_panel_variants, summarize_hs_panel_readiness
from core.run_manifest_validation import main as run_manifest_validation_main


class ManifestRegistryTests(unittest.TestCase):
    def setUp(self) -> None:
        self.panel = load_hs_variant_panel("config/hs_variant_panel.yaml")
        self.manifest = load_data_manifest("data/manifests/hs_gag_panel.json")

    def test_annotate_hs_panel_variants_merges_panel_and_manifest_fields(self) -> None:
        annotations = annotate_hs_panel_variants(self.panel, self.manifest)

        self.assertIn("HS-dp4-NS-2S", annotations)
        self.assertEqual(
            annotations["HS-dp4-NS-2S"].glytoucan_accession,
            "G13487SU",
        )
        self.assertEqual(
            annotations["HS-dp4-NS-2S"].atomistic_preparation_status,
            "registered_structure_resolved",
        )

    def test_summarize_hs_panel_readiness_reports_unresolved_variants(self) -> None:
        summary = summarize_hs_panel_readiness(self.panel, self.manifest)

        self.assertEqual(summary["panel_variant_count"], 11)
        self.assertEqual(summary["matched_variant_count"], 11)
        self.assertEqual(len(summary["resolved_variant_ids"]), 1)
        self.assertEqual(len(summary["unresolved_variant_ids"]), 10)
        self.assertAlmostEqual(summary["readiness_fraction"], 1.0 / 11.0)

    def test_run_manifest_validation_writes_summary(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            output_path = Path(temp_dir) / "manifest_summary.md"
            with patch(
                "sys.argv",
                [
                    "run_manifest_validation.py",
                    "--output",
                    str(output_path),
                ],
            ):
                exit_code = run_manifest_validation_main()

            self.assertEqual(exit_code, 0)
            self.assertTrue(output_path.exists())
            summary_text = output_path.read_text(encoding="utf-8")
        self.assertIn("# Manifest Validation Summary", summary_text)


if __name__ == "__main__":
    unittest.main()
