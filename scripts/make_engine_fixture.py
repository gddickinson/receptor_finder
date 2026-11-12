"""Write tests/fixtures/engine_fixture.json: inputs and the repo's own Python outputs for the JS engine tests.

The node tests (web/static/*_engine.test.mjs) replay every case through the JS ports and require the same answers.
tests/test_engine_fixture.py checks that this file is still what the current Python code produces.

Run: python scripts/make_engine_fixture.py
"""
import json
import sys
from itertools import product
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from binding_pocket import BindingPocket, PocketCharacterizer, PocketScorer  # noqa: E402
from docking import DockingAnalyzer  # noqa: E402
from phylogenetics import PLANT_SPECIES, classify_distribution  # noqa: E402
from scoring import (EvolutionaryScorer, ExpressionScorer, LocalizationScorer,  # noqa: E402
                     StructuralScorer, TopologyScorer)
from sequence_analysis import MotifScanner, TransmembraneDomainPredictor  # noqa: E402
from utils import read_fasta  # noqa: E402

OUT = ROOT / "tests" / "fixtures" / "engine_fixture.json"
AA = "ACDEFGHIKLMNPQRSTVWY"
KD = TransmembraneDomainPredictor.HYDROPHOBICITY


def _profile(seq, window=19):
    s = seq.upper()
    return [float(np.mean([KD.get(a, 0) for a in s[i:i + window]])) for i in range(len(s) - window + 1)]


def sequence_cases():
    rng = np.random.default_rng(20261112)
    seqs = dict(read_fasta(str(ROOT / "data" / "example_sequences.fasta")))
    bench = json.loads((ROOT / "tests" / "fixtures" / "uniprot_tm_benchmark.json").read_text())
    seqs.update({e["accession"]: e["sequence"] for e in bench["entries"]})
    ef = json.loads((ROOT / "tests" / "fixtures" / "uniprot_efhand_benchmark.json").read_text())
    seqs.update({e["accession"]: e["sequence"] for e in ef["entries"]})
    for n in (15, 20, 22, 23, 33, 40):
        seqs[f"polyL_{n}"] = "K" * 20 + "L" * n + "K" * 20
    seqs["two_close"] = "K" * 10 + "L" * 35 + "KK" + "L" * 35 + "K" * 10
    seqs["c_term"] = "K" * 5 + "L" * 40
    seqs["lower_xb"] = "mkllvvaxbllllIIIIffffvvvaaaLLLGGGkkkrrrDDDEEE*"
    hyd, pol = "LIVFAMCG", "KRDENQSTHPYW"
    for k in range(12):  # random proteins with hydrophobic segments of varied length
        parts = []
        for _ in range(rng.integers(2, 6)):
            parts.append("".join(rng.choice(list(pol), rng.integers(5, 30))))
            parts.append("".join(rng.choice(list(hyd), rng.integers(10, 45))))
        seqs[f"random_{k}"] = "".join(parts)
    for k in range(4):
        seqs[f"uniform_{k}"] = "".join(rng.choice(list(AA), int(rng.integers(30, 400))))
    cases = []
    for sid, seq in seqs.items():
        cases.append({
            "id": sid, "sequence": seq, "profile": _profile(seq),
            "tm": [list(t) for t in TransmembraneDomainPredictor.predict_tm_domains(seq)],
            "tm_w15_t1": [list(t) for t in TransmembraneDomainPredictor.predict_tm_domains(seq, 15, 1.0)],
            "tm_minw15_gap0": [list(t) for t in TransmembraneDomainPredictor.predict_tm_domains(seq, 19, 1.6, 15, 0)],
            "motifs": {m: MotifScanner.scan_sequence(seq, p) for m, p in MotifScanner.MOTIFS.items()},
        })
    return cases


def score_cases():
    aff = [-12, -9.0, -8.99, -8.2, -7.0, -6.5, -5.0, -4.2, -1, 0, 3]
    pl = [(m, t) for m in (20, 50, 50.5, 63, 70, 70.1, 85, 90, 95) for t in (None, 40, 69.9, 70, 88)]
    cn = [(c, a) for c in (0, 1, 2, 3, 4, 5, 9, 10, 30) for a in (0, 1, 2, 3)]
    corr = [-0.5, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.81, 1]
    tissues = [(["guard_cell"], ["guard_cell", "root_hair"]), ([], []), (["a", "b"], ["a", "b", "b"]),
               (["x"], ["y"]), (["p", "q", "r"], ["p", "q", "r"])]
    locs = list(LocalizationScorer.LOCALIZATION_SCORES) + ["Tonoplast", "er", "unknown"]
    sig = [(h, t) for h in (False, True) for t in (None, *LocalizationScorer.SIGNAL_TYPE_SCORES, "other")]
    tm = [(n, r) for n in range(0, 13) for r in ((4, 8), (2, 10), (6, 6))]
    return {
        "docking": [[a, StructuralScorer.score_docking(a)] for a in aff],
        "quality": [[m, t, StructuralScorer.score_structure_quality(m, t)] for m, t in pl],
        "pocket": [[r, mx, StructuralScorer.score_pocket(r, mx)] for r in (0, 7.5, 13, 20) for mx in (12, 13)],
        "pattern": [[p, EvolutionaryScorer.score_phylogenetic_pattern(p)]
                    for p in list(EvolutionaryScorer.PATTERN_SCORES) + ["nonsense"]],
        "expansion": [[c, a, EvolutionaryScorer.score_gene_expansion(c, a)] for c, a in cn],
        "coexpression": [[r, ExpressionScorer.score_coexpression(r)] for r in corr],
        "tissue": [[t, e, ExpressionScorer.score_tissue_specificity(t, e)] for t, e in tissues],
        "localization": [[x, LocalizationScorer.score_localization(x)] for x in locs],
        "signal": [[h, t, LocalizationScorer.score_signal_peptide(h, t)] for h, t in sig],
        "tm": [[n, list(r), TopologyScorer.score_tm_domains(n, r)] for n, r in tm],
        "pore": [[a, b, TopologyScorer.score_pore_architecture(a, b)] for a in (False, True) for b in (False, True)],
    }


def pocket_cases():
    rng = np.random.default_rng(7)
    names = ["ARG", "LYS", "HIS", "ASP", "GLU", "PHE", "TRP", "TYR", "SER", "THR", "ASN", "GLN", "CYS", "LEU", "VAL"]
    out = []
    for k in range(40):
        p = BindingPocket(pocket_id=f"P{k}", protein_id="X",
                          volume=float(rng.choice([None, 150, 200, 250, 300, 301, 450]) or 0) or None,
                          depth=float(rng.choice([3, 5, 5.5, 9])), charge=float(rng.choice([-2, 0, 0.5, 2, 2.5, 4])),
                          hydrophobicity=float(rng.choice([0, 0.2, 0.3, 0.4, 0.5])),
                          residues=list(rng.choice(names, int(rng.integers(0, 12)))))
        i3, cd = PocketScorer.score_insp3_binding(p), PocketScorer.score_cadpr_binding(p)
        i4 = [None, i3 - 3, i3 - 2, i3][k % 4]
        nad = [None, cd - 3, cd - 2, cd, cd - 2.5][k % 5]
        out.append({"volume": p.volume, "depth": p.depth, "charge": p.charge, "hydrophobicity": p.hydrophobicity,
                    "residues": p.residues, "insp3": i3, "cadpr": cd, "insp4": i4, "nad": nad,
                    "counts": [p.basic_count, p.aromatic_count, p.acidic_count, p.polar_count],
                    "specificity": PocketCharacterizer.assess_specificity(p, i3, cd, i4, nad)})
    return out


def phylo_cases():
    species = list(PLANT_SPECIES)
    out = []
    for bits in product([False, True], repeat=len(species)):
        presence = dict(zip(species, bits))
        pat = classify_distribution("G", presence)
        out.append({"presence": presence, "pattern": pat.pattern, "confidence": pat.confidence})
    return out


def docking_cases():
    pairs = [(-9, -6), (-8, -7), (-7, -9.5), (-6.3, -4.2), (-10, -7.9), (-5, -5)]
    return [DockingAnalyzer.assess_specificity(a, b) for a, b in pairs]


def candidate_cases():
    """Whole candidates through POST /api/score's composition (web/app.py)."""
    sys.path.insert(0, str(ROOT))
    from web.app import ScoreIn, score
    rng = np.random.default_rng(11)
    pats, locs = list(EvolutionaryScorer.PATTERN_SCORES), list(LocalizationScorer.LOCALIZATION_SCORES)
    out = [{"input": {}, "output": score(ScoreIn())}]
    for _ in range(30):
        inp = {"affinity": float(rng.uniform(-12, 1)), "mean_plddt": float(rng.uniform(20, 98)),
               "tm_region_plddt": [None, float(rng.uniform(30, 95))][int(rng.integers(0, 2))],
               "pocket_raw": float(rng.uniform(0, 15)), "pattern": pats[int(rng.integers(0, len(pats)))],
               "copies": int(rng.integers(0, 15)), "ancestral": int(rng.integers(0, 4)),
               "correlation": float(rng.uniform(-0.3, 1)), "tissues": ["guard_cell", "pollen_tube"][:int(rng.integers(0, 3))],
               "expected_tissues": ["guard_cell", "root_hair", "pollen_tube"][:int(rng.integers(0, 4))],
               "localization": locs[int(rng.integers(0, len(locs)))], "has_signal": bool(rng.integers(0, 2)),
               "signal_type": [None, "ER", "vacuolar", "other"][int(rng.integers(0, 4))], "n_tm": int(rng.integers(0, 13)),
               "pore_loop": bool(rng.integers(0, 2)), "selectivity_filter": bool(rng.integers(0, 2)),
               "topology_match": float(rng.uniform(-1, 6))}
        out.append({"input": inp, "output": score(ScoreIn(**inp))})
    return out


def build():
    return {"sequences": sequence_cases(), "scores": score_cases(), "pockets": pocket_cases(),
            "phylo": phylo_cases(), "docking": docking_cases(), "candidates": candidate_cases()}


if __name__ == "__main__":
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(build(), separators=(",", ":")))
    print(f"wrote {OUT} ({OUT.stat().st_size / 1e3:.0f} kB)")
