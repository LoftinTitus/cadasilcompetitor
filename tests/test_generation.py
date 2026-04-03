from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from peptide.generator import generate_candidates_from_motifs
from peptide.motifs import load_seed_motifs
from peptide.schema import SeedMotif


class GenerationTests(unittest.TestCase):
    def test_load_seed_motifs_supports_alias_headers_and_dedupes(self) -> None:
        csv_text = "\n".join(
            [
                "Identifier,Sequence,Origin,Type,Comments,Natural",
                "motif_a, b x r ,Example Protein,pattern,first motif,yes",
                "motif_dup,BXR,Example Protein,pattern,first motif,yes",
                "motif_b,AKRK,Synthetic Panel,exact,exact motif,no",
                ",,, , ,",
            ]
        )

        with tempfile.TemporaryDirectory() as temp_dir:
            motif_path = Path(temp_dir) / "motifs.csv"
            motif_path.write_text(csv_text, encoding="utf-8")
            motifs = load_seed_motifs(motif_path)

        self.assertEqual(len(motifs), 2)
        self.assertEqual(motifs[0].normalized_pattern, "BXR")
        self.assertEqual(motifs[0].pattern_type, "pattern")
        self.assertTrue(motifs[0].is_natural)
        self.assertEqual(motifs[1].normalized_pattern, "AKRK")
        self.assertFalse(motifs[1].is_natural)

    def test_generate_candidates_deduplicates_sequences_and_merges_motif_ids(self) -> None:
        motifs = [
            SeedMotif(
                motif_id="motif_a",
                raw_pattern="AKRK",
                normalized_pattern="AKRK",
                pattern_type="exact",
                source="Seed A",
            ),
            SeedMotif(
                motif_id="motif_b",
                raw_pattern="AKRK",
                normalized_pattern="AKRK",
                pattern_type="exact",
                source="Seed B",
            ),
        ]

        candidates = generate_candidates_from_motifs(
            motifs,
            max_variants_per_motif=3,
            max_total_candidates=10,
            random_seed=3,
        )

        self.assertGreaterEqual(len(candidates), 1)
        self.assertEqual(candidates[0].candidate_id, "pep_00001")
        self.assertEqual(candidates[0].source_motif_ids, ["motif_a", "motif_b"])
        self.assertEqual(len({candidate.sequence for candidate in candidates}), len(candidates))

    def test_generate_candidates_expands_pattern_motifs(self) -> None:
        motifs = [
            SeedMotif(
                motif_id="motif_pattern",
                raw_pattern="BXR",
                normalized_pattern="BXR",
                pattern_type="pattern",
                source="Pattern seed",
            )
        ]

        candidates = generate_candidates_from_motifs(
            motifs,
            max_variants_per_motif=6,
            max_total_candidates=20,
            random_seed=7,
        )

        sequences = [candidate.sequence for candidate in candidates]
        self.assertTrue(any("KGR" in sequence for sequence in sequences))
        self.assertTrue(
            all(candidate.generation_method.startswith("pattern_expansion_") for candidate in candidates)
        )
        self.assertTrue(all(candidate.length >= 8 for candidate in candidates))


if __name__ == "__main__":
    unittest.main()
