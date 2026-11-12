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


class TestDockingScoreDirection(unittest.TestCase):
    """score_docking maps -9..-7 kcal/mol to 10..7, -7..-5 to 7..3 and above -5 down to 0 (its docstring)."""

    def test_documented_anchor_points(self):
        from scoring import StructuralScorer as S
        for aff, pts in ((-10, 10), (-9, 10), (-8, 8.5), (-7, 7), (-6, 5), (-5, 3), (0, 0), (2, 0)):
            self.assertAlmostEqual(S.score_docking(aff), pts, msg=aff)

    def test_never_above_ten_and_monotonic(self):
        from scoring import StructuralScorer as S
        affs = [x / 10 for x in range(-120, 30)]
        vals = [S.score_docking(a) for a in affs]
        self.assertTrue(all(0 <= v <= 10 for v in vals))
        self.assertTrue(all(a >= b for a, b in zip(vals, vals[1:])), "a better (more negative) affinity must not score lower")


class TestEfHandProsite(unittest.TestCase):
    """EF_hand is PROSITE PS00018 (EF_HAND_1) and finds calmodulin's four Ca2+-binding loops."""

    LOOPS = ["DKDGDGTITTKE", "DADGNGTIDFPE", "DKDGNGYISAAE", "DIDGDGQVNYEE"]

    def _cam(self):
        from utils import read_fasta
        return read_fasta(str(Path(__file__).resolve().parent.parent / "data" / "example_sequences.fasta"))["AT2G23456"]

    def test_all_four_calmodulin_loops(self):
        from sequence_analysis import MotifScanner
        seq = self._cam()
        hits = MotifScanner.scan_sequence(seq, MotifScanner.MOTIFS["EF_hand"])
        self.assertEqual(hits, [seq.find(lp) for lp in self.LOOPS])

    def test_loop_without_first_aspartate_is_rejected(self):
        from sequence_analysis import MotifScanner
        seq = self._cam().replace("DKDGDGTITTKEL", "AKDGDGTITTKEL")
        self.assertEqual(len(MotifScanner.scan_sequence(seq, MotifScanner.MOTIFS["EF_hand"])), 3)


class TestKyteDoolittleRule(unittest.TestCase):
    """Standard Kyte & Doolittle (1982) rule: any 19-residue window with mean hydropathy >= 1.6 marks a membrane-spanning
    segment. Benchmarked on UniProt-annotated TM proteins and soluble decoys (tests/fixtures/uniprot_tm_benchmark.json)."""

    @staticmethod
    def _bench():
        import json
        p = Path(__file__).resolve().parent / "fixtures" / "uniprot_tm_benchmark.json"
        return {e["accession"]: e for e in json.loads(p.read_text())["entries"]}

    @staticmethod
    def _found(pred, annotated):
        return sum(any(a < e and s - 1 < b for a, b in pred) for s, e in annotated)

    def test_twenty_residue_helix_is_found(self):
        tm = TransmembraneDomainPredictor.predict_tm_domains("K" * 20 + "L" * 20 + "K" * 20)
        self.assertEqual(len(tm), 1)

    def test_known_tm_proteins(self):
        b = self._bench()
        expect = {"P02724": 1, "P02945": 5, "P29972": 5, "P0AEJ4": 2}  # annotated: 1, 7, 6, 2
        for acc, n in expect.items():
            pred = TransmembraneDomainPredictor.predict_tm_domains(b[acc]["sequence"])
            self.assertEqual(self._found(pred, b[acc]["transmembrane"]), n, acc)

    def test_soluble_decoys_have_none(self):
        b = self._bench()
        for acc in ("P0DP23", "P02144"):
            self.assertEqual(TransmembraneDomainPredictor.predict_tm_domains(b[acc]["sequence"]), [], acc)

    def test_old_rule_still_available(self):
        seq = "K" * 20 + "L" * 20 + "K" * 20
        self.assertEqual(TransmembraneDomainPredictor.predict_tm_domains(seq, min_windows=15), [])


class TestPhylogeneticPatterns(unittest.TestCase):
    """classify_distribution returns every documented pattern (docs/in_silico_calcium_projects.md, section 3.1)."""

    SPECIES = ["Arabidopsis_thaliana", "Oryza_sativa", "Zea_mays", "Marchantia_polymorpha", "Physcomitrella_patens",
               "Chlamydomonas_reinhardtii", "Selaginella_moellendorffii", "Amborella_trichopoda"]

    def _p(self, *present):
        from phylogenetics import classify_distribution
        return classify_distribution("G", {s: s in present for s in self.SPECIES}).pattern

    def test_conserved(self):
        self.assertEqual(self._p(*self.SPECIES), "Conserved")

    def test_chlamydomonas_loss(self):
        # present in the alga, absent from early land plants, present in angiosperms
        self.assertEqual(self._p("Chlamydomonas_reinhardtii", "Arabidopsis_thaliana", "Oryza_sativa"), "Chlamydomonas_loss")

    def test_land_plant_innovation(self):
        self.assertEqual(self._p("Marchantia_polymorpha", "Physcomitrella_patens", "Selaginella_moellendorffii",
                                 "Arabidopsis_thaliana", "Zea_mays"), "Land_plant_innovation")

    def test_angiosperm_specific(self):
        self.assertEqual(self._p("Arabidopsis_thaliana", "Oryza_sativa", "Amborella_trichopoda"), "Angiosperm_specific")

    def test_bryophyte_origin(self):
        self.assertEqual(self._p("Marchantia_polymorpha", "Physcomitrella_patens"), "Bryophyte_origin")

    def test_sporadic(self):
        self.assertEqual(self._p("Chlamydomonas_reinhardtii", "Marchantia_polymorpha"), "Sporadic")

    def test_unknown_when_absent_everywhere(self):
        self.assertEqual(self._p(), "Unknown")

    def test_every_pattern_is_reachable_and_scored(self):
        from itertools import product
        from phylogenetics import classify_distribution
        from scoring import EvolutionaryScorer
        seen = {classify_distribution("G", dict(zip(self.SPECIES, bits))).pattern
                for bits in product([False, True], repeat=len(self.SPECIES))}
        self.assertEqual(seen, set(EvolutionaryScorer.PATTERN_SCORES))


def _heavy_atoms(smiles):
    """Heavy-atom counts from a SMILES string (enough to compare with a molecular formula)."""
    import re
    counts = {}
    for tok in re.findall(r"\[[^\]]+\]|Cl|Br|[BCNOPSFI]|[cnops]", smiles):
        el = re.match(r"\[?\d*([A-Z][a-z]?|[cnops])", tok).group(1).capitalize()
        counts[el] = counts.get(el, 0) + 1
    return counts


class TestLigandLibrary(unittest.TestCase):
    """SMILES are PubChem's (see THIRD_PARTY.md); heavy atoms must match the PubChem molecular formula."""

    FORMULA = {  # PubChem CID: molecular formula without H
        "InsP3": {"C": 6, "O": 15, "P": 3},           # CID 439456 C6H15O15P3
        "InsP4": {"C": 6, "O": 18, "P": 4},           # CID 107758 C6H16O18P4
        "cADPR": {"C": 15, "N": 5, "O": 13, "P": 2},  # CID 123847 C15H21N5O13P2
        "NAADP": {"C": 21, "N": 6, "O": 18, "P": 3},  # CID 123952 C21H27N6O18P3
        "NAD": {"C": 21, "N": 7, "O": 14, "P": 2},    # CID 5892 C21H27N7O14P2
        "8-NH2-cADPR": {"C": 15, "N": 6, "O": 13, "P": 2},  # ChEMBL CHEMBL3350988, see THIRD_PARTY.md
    }

    def test_formulas(self):
        from docking import LigandLibrary
        for name, f in self.FORMULA.items():
            self.assertEqual(_heavy_atoms(LigandLibrary.LIGANDS[name]), f, name)

    def test_bases(self):
        from docking import LigandLibrary
        adenine_like = "N=CN"  # purine ring fragment present in every adenine-containing SMILES here
        for name in ("cADPR", "NAADP", "NAD", "8-NH2-cADPR"):
            smi = LigandLibrary.LIGANDS[name]
            self.assertIn(adenine_like, smi, name)
            self.assertNotIn("C=CC(=NC", smi, name)  # cytosine ring fragment of the old cADPR SMILES
        self.assertIn("C(=O)O", LigandLibrary.LIGANDS["NAADP"])  # nicotinic acid, not nitro
        self.assertNotIn("[N+](=O)[O-]", LigandLibrary.LIGANDS["NAADP"])
        self.assertTrue(LigandLibrary.LIGANDS["NAD"].endswith("C(=O)N"))  # nicotinamide

    def test_inositol_ring_is_carbocyclic(self):
        from docking import LigandLibrary
        for name in ("InsP3", "InsP4"):
            self.assertNotIn("O1", LigandLibrary.LIGANDS[name], name)
            self.assertNotIn("N", LigandLibrary.LIGANDS[name], name)
