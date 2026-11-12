"""Regression tests for bugs fixed during the web-interface work (2026-11-12)."""
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from sequence_analysis import TransmembraneDomainPredictor  # noqa: E402
from scoring import CandidateRanker, CandidateScore  # noqa: E402
from structure_prediction import AlphaFold2Interface  # noqa: E402


class TestTmDomainCoordinates(unittest.TestCase):
    """TM domain ends are exclusive and lie inside the sequence."""

    def test_end_inside_sequence_at_c_terminus(self):
        seq = "K" * 5 + "L" * 40
        tm = TransmembraneDomainPredictor.predict_tm_domains(seq)
        self.assertEqual(len(tm), 1)
        self.assertEqual(tm[0][1], len(seq))

    def test_internal_domain_end_is_last_window_end(self):
        seq = "K" * 20 + "L" * 30 + "K" * 20
        w = 19
        prof = [np.mean([TransmembraneDomainPredictor.HYDROPHOBICITY[a] for a in seq[i:i + w]])
                for i in range(len(seq) - w + 1)]
        above = [i for i, h in enumerate(prof) if h >= 1.6]
        tm = TransmembraneDomainPredictor.predict_tm_domains(seq)
        self.assertEqual(tm, [(above[0], above[-1] + w)])


class TestRankByField(unittest.TestCase):
    """rank_candidates accepts every score field, not only the *_score ones."""

    def test_rank_by_structure_quality(self):
        r = CandidateRanker()
        r.add_candidate(CandidateScore(protein_id="low", structure_quality=1))
        r.add_candidate(CandidateScore(protein_id="high", structure_quality=6))
        self.assertEqual(r.rank_candidates("structure_quality")[0].protein_id, "high")

    def test_rank_by_short_name_still_works(self):
        r = CandidateRanker()
        r.add_candidate(CandidateScore(protein_id="a", pocket_score=2))
        r.add_candidate(CandidateScore(protein_id="b", pocket_score=9))
        self.assertEqual(r.rank_candidates("pocket")[0].protein_id, "b")

    def test_unknown_field_raises(self):
        with self.assertRaises(ValueError):
            CandidateRanker().rank_candidates("nonsense")


class TestMockStructureIsLabelled(unittest.TestCase):
    """The placeholder predictor must not label its output as an AlphaFold2 prediction."""

    def test_mock_pdb_header(self):
        with tempfile.TemporaryDirectory() as d:
            s = AlphaFold2Interface(output_dir=d).predict_structure("MAAALLL" * 5, "X")
            text = Path(s.pdb_file).read_text()
        self.assertNotIn("AlphaFold2 PREDICTION", text)
        self.assertIn("PLACEHOLDER", text)

    def test_committed_structures_are_labelled(self):
        root = Path(__file__).resolve().parent.parent / "outputs" / "structures"
        for pdb in root.glob("*.pdb"):
            self.assertIn("PLACEHOLDER", pdb.read_text().splitlines()[1], pdb.name)


if __name__ == "__main__":
    unittest.main()


class TestLocalizationCase(unittest.TestCase):
    """score_localization lower-cases its input, so mixed-case keys must still match."""

    def test_er_scores_nine(self):
        from scoring import LocalizationScorer
        for name in ("ER", "er", "Er"):
            self.assertEqual(LocalizationScorer.score_localization(name), 9.0)

    def test_other_keys_unchanged(self):
        from scoring import LocalizationScorer
        self.assertEqual(LocalizationScorer.score_localization("Tonoplast"), 10.0)
        self.assertEqual(LocalizationScorer.score_localization("plasma_membrane"), 6.0)
        self.assertEqual(LocalizationScorer.score_localization("unknown"), 0.0)
