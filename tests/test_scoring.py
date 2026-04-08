"""
Unit tests for the scoring module.

Verifies score ranges, priority classification thresholds,
and ranking behaviour.
"""

import sys
import unittest
from pathlib import Path

# Ensure src is importable
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from scoring import (
    CandidateScore,
    StructuralScorer,
    EvolutionaryScorer,
    ExpressionScorer,
    LocalizationScorer,
    TopologyScorer,
    CandidateRanker,
)


class TestCandidateScore(unittest.TestCase):
    """Tests for CandidateScore dataclass."""

    def test_total_score_sums_all_components(self):
        cs = CandidateScore(
            protein_id="TEST",
            pocket_score=10.0,
            docking_score=8.0,
            structure_quality=5.0,
            phylo_pattern_score=6.0,
            expansion_score=4.0,
            coexpression_score=3.0,
            tissue_specificity=2.0,
            localization_score=7.0,
            signal_peptide_score=3.0,
            tm_domain_score=4.0,
            pore_architecture=2.0,
            topology_match=1.0,
        )
        self.assertAlmostEqual(cs.total_score, 55.0)

    def test_priority_high(self):
        cs = CandidateScore(protein_id="H", pocket_score=13, docking_score=10,
                            structure_quality=7, phylo_pattern_score=10,
                            expansion_score=10, coexpression_score=10)
        self.assertEqual(cs.priority_class, "HIGH")

    def test_priority_medium(self):
        cs = CandidateScore(protein_id="M", pocket_score=10, docking_score=8,
                            structure_quality=5, phylo_pattern_score=6,
                            expansion_score=5, coexpression_score=6)
        # total = 40
        self.assertEqual(cs.priority_class, "MEDIUM")

    def test_priority_low(self):
        cs = CandidateScore(protein_id="L", pocket_score=2, docking_score=1)
        self.assertEqual(cs.priority_class, "LOW")

    def test_zero_score(self):
        cs = CandidateScore(protein_id="ZERO")
        self.assertAlmostEqual(cs.total_score, 0.0)
        self.assertEqual(cs.priority_class, "LOW")


class TestStructuralScorer(unittest.TestCase):
    """Tests for StructuralScorer static methods."""

    def test_pocket_score_clamped(self):
        self.assertAlmostEqual(StructuralScorer.score_pocket(15.0), 13.0)

    def test_pocket_score_passthrough(self):
        self.assertAlmostEqual(StructuralScorer.score_pocket(8.0), 8.0)

    def test_docking_excellent(self):
        score = StructuralScorer.score_docking(-10.0)
        self.assertAlmostEqual(score, 10.0)

    def test_docking_poor(self):
        score = StructuralScorer.score_docking(-2.0)
        self.assertGreaterEqual(score, 0.0)
        self.assertLessEqual(score, 10.0)  # Must be within valid range

    def test_structure_quality_range(self):
        for plddt in [20, 50, 70, 85, 95]:
            score = StructuralScorer.score_structure_quality(float(plddt))
            self.assertGreaterEqual(score, 0.0)
            self.assertLessEqual(score, 7.0)


class TestEvolutionaryScorer(unittest.TestCase):
    def test_known_pattern(self):
        self.assertAlmostEqual(
            EvolutionaryScorer.score_phylogenetic_pattern("Chlamydomonas_loss"), 10.0
        )

    def test_unknown_pattern(self):
        self.assertAlmostEqual(
            EvolutionaryScorer.score_phylogenetic_pattern("Nonexistent"), 0.0
        )

    def test_gene_expansion_high(self):
        score = EvolutionaryScorer.score_gene_expansion(12, 1)
        self.assertAlmostEqual(score, 10.0)

    def test_gene_expansion_no_expansion(self):
        score = EvolutionaryScorer.score_gene_expansion(1, 1)
        self.assertAlmostEqual(score, 2.5)


class TestTopologyScorer(unittest.TestCase):
    def test_tm_in_range(self):
        self.assertAlmostEqual(TopologyScorer.score_tm_domains(6), 5.0)

    def test_tm_just_outside(self):
        self.assertAlmostEqual(TopologyScorer.score_tm_domains(3), 3.0)

    def test_pore_full(self):
        self.assertAlmostEqual(
            TopologyScorer.score_pore_architecture(True, True), 5.0
        )


class TestCandidateRanker(unittest.TestCase):
    def test_ranking_order(self):
        ranker = CandidateRanker()
        c1 = CandidateScore(protein_id="A", pocket_score=10)
        c2 = CandidateScore(protein_id="B", pocket_score=5)
        ranker.add_candidate(c1)
        ranker.add_candidate(c2)
        ranked = ranker.rank_candidates()
        self.assertEqual(ranked[0].protein_id, "A")

    def test_high_priority_filter(self):
        ranker = CandidateRanker()
        ranker.add_candidate(CandidateScore(protein_id="H", pocket_score=13,
                                            docking_score=10, structure_quality=7,
                                            phylo_pattern_score=10, expansion_score=10,
                                            coexpression_score=10))
        ranker.add_candidate(CandidateScore(protein_id="L", pocket_score=1))
        high = ranker.get_high_priority(threshold=60)
        self.assertEqual(len(high), 1)
        self.assertEqual(high[0].protein_id, "H")

    def test_to_dataframe(self):
        ranker = CandidateRanker()
        ranker.add_candidate(CandidateScore(protein_id="X"))
        df = ranker.to_dataframe()
        self.assertEqual(len(df), 1)
        self.assertIn("total_score", df.columns)


if __name__ == "__main__":
    unittest.main()
