from __future__ import annotations

import unittest

from models.optimization import DEFAULT_TARGET_FIELDS, pareto_frontier, propose_candidates
from models.surrogate import fit_bootstrap_surrogate


def _training_row(
    *,
    candidate_id: str,
    sequence: str,
    hs_affinity_reward: float,
    hs_selectivity_reward: float,
    transport_reward: float,
    residence_reward: float,
    anticoagulant_penalty: float,
    off_target_penalty: float,
    heparin_penalty: float,
    developability_penalty: float,
) -> dict[str, object]:
    composite = 100.0 * (
        (0.35 * hs_affinity_reward)
        + (0.20 * hs_selectivity_reward)
        + (0.20 * transport_reward)
        + (0.10 * residence_reward)
        - (0.10 * anticoagulant_penalty)
        - (0.07 * off_target_penalty)
        - (0.03 * heparin_penalty)
        - (0.05 * developability_penalty)
    )
    return {
        "candidate_id": candidate_id,
        "sequence": sequence,
        "warning_flags": [],
        "filter_flags": [],
        "hs_affinity_reward": hs_affinity_reward,
        "hs_selectivity_reward": hs_selectivity_reward,
        "transport_reward": transport_reward,
        "residence_reward": residence_reward,
        "anticoagulant_penalty": anticoagulant_penalty,
        "off_target_penalty": off_target_penalty,
        "heparin_penalty": heparin_penalty,
        "developability_penalty": developability_penalty,
        "composite_screen_score": composite,
    }


class OptimizationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.training_rows = [
            _training_row(
                candidate_id="train_1",
                sequence="AKRKRQGK",
                hs_affinity_reward=0.82,
                hs_selectivity_reward=0.30,
                transport_reward=0.58,
                residence_reward=0.44,
                anticoagulant_penalty=0.08,
                off_target_penalty=0.09,
                heparin_penalty=0.04,
                developability_penalty=0.05,
            ),
            _training_row(
                candidate_id="train_2",
                sequence="GKKQGKQEGKKQGKQ",
                hs_affinity_reward=0.56,
                hs_selectivity_reward=0.24,
                transport_reward=0.41,
                residence_reward=0.36,
                anticoagulant_penalty=0.12,
                off_target_penalty=0.11,
                heparin_penalty=0.05,
                developability_penalty=0.08,
            ),
            _training_row(
                candidate_id="train_3",
                sequence="AKRGRKRRKQGA",
                hs_affinity_reward=0.76,
                hs_selectivity_reward=0.27,
                transport_reward=0.53,
                residence_reward=0.40,
                anticoagulant_penalty=0.10,
                off_target_penalty=0.10,
                heparin_penalty=0.04,
                developability_penalty=0.06,
            ),
        ]

    def test_fit_bootstrap_surrogate_predicts_target_summaries(self) -> None:
        surrogate = fit_bootstrap_surrogate(
            self.training_rows,
            target_fields=DEFAULT_TARGET_FIELDS,
            ensemble_size=5,
            random_seed=11,
        )

        prediction = surrogate.predict({"sequence": "AKRGRQKRRKA"})

        self.assertEqual(set(prediction), set(DEFAULT_TARGET_FIELDS))
        self.assertIn("mean", prediction["hs_affinity_reward"])
        self.assertIn("std", prediction["transport_reward"])

    def test_propose_candidates_filters_seen_and_flagged_rows(self) -> None:
        surrogate = fit_bootstrap_surrogate(
            self.training_rows,
            target_fields=DEFAULT_TARGET_FIELDS,
            ensemble_size=5,
            random_seed=13,
        )
        candidate_pool = [
            {
                "candidate_id": "seen_duplicate",
                "sequence": "AKRKRQGK",
                "warning_flags": [],
                "filter_flags": [],
            },
            {
                "candidate_id": "filtered_out",
                "sequence": "RRRRKKKKWW",
                "warning_flags": ["aggregation_risk"],
                "filter_flags": ["basic_fraction_out_of_range"],
            },
            {
                "candidate_id": "proposal_a",
                "sequence": "AKRGRQKRRKA",
                "warning_flags": [],
                "filter_flags": [],
            },
            {
                "candidate_id": "proposal_b",
                "sequence": "GKRRKAKRGRR",
                "warning_flags": ["cpp_like_uptake_risk"],
                "filter_flags": [],
            },
        ]

        proposals = propose_candidates(
            candidate_pool,
            trained_surrogate=surrogate,
            reference_rows=self.training_rows,
            top_k=10,
        )

        self.assertEqual([row["proposal_rank"] for row in proposals], [1, 2])
        self.assertEqual([row["candidate_id"] for row in proposals], ["proposal_a", "proposal_b"])
        self.assertEqual(proposals[0]["filter_flags"], "")
        self.assertEqual(proposals[1]["warning_flags"], "cpp_like_uptake_risk")

    def test_pareto_frontier_returns_only_non_dominated_rows(self) -> None:
        frontier = pareto_frontier(
            [
                {
                    "candidate_id": "dominant",
                    "composite_screen_score": 60.0,
                    "hs_affinity_reward": 0.80,
                    "hs_selectivity_reward": 0.35,
                    "transport_reward": 0.62,
                    "residence_reward": 0.42,
                    "anticoagulant_penalty": 0.06,
                    "off_target_penalty": 0.08,
                    "heparin_penalty": 0.04,
                    "developability_penalty": 0.05,
                },
                {
                    "candidate_id": "dominated",
                    "composite_screen_score": 52.0,
                    "hs_affinity_reward": 0.72,
                    "hs_selectivity_reward": 0.30,
                    "transport_reward": 0.50,
                    "residence_reward": 0.38,
                    "anticoagulant_penalty": 0.10,
                    "off_target_penalty": 0.10,
                    "heparin_penalty": 0.05,
                    "developability_penalty": 0.06,
                },
                {
                    "candidate_id": "tradeoff",
                    "composite_screen_score": 58.0,
                    "hs_affinity_reward": 0.74,
                    "hs_selectivity_reward": 0.33,
                    "transport_reward": 0.66,
                    "residence_reward": 0.40,
                    "anticoagulant_penalty": 0.07,
                    "off_target_penalty": 0.11,
                    "heparin_penalty": 0.04,
                    "developability_penalty": 0.05,
                },
            ]
        )

        self.assertEqual(
            [row["candidate_id"] for row in frontier],
            ["dominant", "tradeoff"],
        )


if __name__ == "__main__":
    unittest.main()
