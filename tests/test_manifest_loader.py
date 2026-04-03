from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from core.manifest_loader import (
    ManifestValidationError,
    load_data_manifest,
    summarize_manifest,
)


class ManifestLoaderTests(unittest.TestCase):
    def test_load_data_manifest_reads_repo_manifest(self) -> None:
        manifest = load_data_manifest("data/manifests/hs_gag_panel.json")

        self.assertEqual(manifest.manifest_id, "hs_gag_panel_v1")
        self.assertEqual(len(manifest.records), 11)
        self.assertEqual(manifest.records[0].record_id, "HS-dp4-NAc")

    def test_summarize_manifest_counts_status_and_source(self) -> None:
        manifest = load_data_manifest("data/manifests/hs_gag_panel.json")

        summary = summarize_manifest(manifest)

        self.assertEqual(summary["record_count"], 11)
        self.assertEqual(summary["by_source_db"]["GlyTouCan"], 11)
        self.assertEqual(summary["by_status"]["needs_programmatic_resolution"], 9)
        self.assertEqual(summary["by_status"]["registered_representative"], 1)

    def test_load_data_manifest_rejects_missing_required_fields(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            manifest_path = Path(temp_dir) / "invalid_manifest.json"
            manifest_path.write_text(
                json.dumps(
                    {
                        "manifest_name": "Broken manifest",
                        "entity_type": "gag_oligo",
                        "source_status": "seed",
                        "generated_by": "test",
                        "records": [],
                    }
                ),
                encoding="utf-8",
            )

            with self.assertRaises(ManifestValidationError):
                load_data_manifest(manifest_path)


if __name__ == "__main__":
    unittest.main()
