"""Benchmarks of the shipped rules, computed by the repo's own Python code (evidence-audit: name what the instrument
cannot do). Each entry says what was tested, against what, and what the code returned. Nothing here changes src/.
"""
from itertools import product

from web import repo_data  # noqa: F401  (puts src/ on sys.path)
from phylogenetics import PLANT_SPECIES, classify_distribution
from sequence_analysis import MotifScanner, TransmembraneDomainPredictor
from utils import read_fasta

# The four 12-residue Ca2+-binding loops (positions 1-12 of each EF hand) inside example sequence AT2G23456, whose
# EF-hand region reproduces calmodulin's. The test suite checks each is a substring of the committed FASTA.
EF_LOOPS = ["DKDGDGTITTKE", "DADGNGTIDFPE", "DKDGNGYISAAE", "DIDGDGQVNYEE"]


def ef_hand():
    """Scan all of AT2G23456 (the PS00018 pattern spans 13 residues: the 12-residue loop plus one hydrophobic residue),
    then check each known loop start is among the hits. Decoy: each loop with its position-1 Asp mutated to Ala."""
    seq = read_fasta(str(repo_data.EXAMPLE_FASTA))["AT2G23456"]
    pat = MotifScanner.MOTIFS["EF_hand"]
    hits = MotifScanner.scan_sequence(seq, pat)
    loops = []
    for lp in EF_LOOPS:
        pos = seq.find(lp)
        mutant = seq[:pos] + "A" + seq[pos + 1:]
        loops.append({"loop": lp, "next": seq[pos + 12], "position": pos, "hit": pos in hits,
                      "d1a_mutant_hit": pos in MotifScanner.scan_sequence(mutant, pat)})
    return {"pattern": pat, "source": "PROSITE PS00018 (EF_HAND_1)", "loops": loops,
            "detected": sum(x["hit"] for x in loops), "decoys_detected": sum(x["d1a_mutant_hit"] for x in loops),
            "hits_in_AT2G23456": hits}


def tm_min_length(flank="K", core="L", flank_len=20):
    rows = []
    for n in range(10, 41):
        seq = flank * flank_len + core * n + flank * flank_len
        rows.append({"n": n, "domains": [list(t) for t in TransmembraneDomainPredictor.predict_tm_domains(seq)]})
    first = next((r["n"] for r in rows if r["domains"]), None)
    return {"construct": f"{flank}x{flank_len} + {core}xn + {flank}x{flank_len}", "rows": rows,
            "shortest_detected": first}


def examples_tm():
    out = []
    text = repo_data.EXAMPLE_FASTA.read_text()
    desc = {ln[1:].split(None, 1)[0]: ln[1:].split(None, 1)[1] for ln in text.splitlines() if ln.startswith(">")}
    for sid, seq in read_fasta(str(repo_data.EXAMPLE_FASTA)).items():
        tm = TransmembraneDomainPredictor.predict_tm_domains(seq)
        out.append({"id": sid, "description": desc.get(sid, ""), "length": len(seq), "n_tm": len(tm),
                    "domains": [list(t) for t in tm]})
    return out


def phylo_reachability():
    species = list(PLANT_SPECIES)
    counts = {}
    for bits in product([False, True], repeat=len(species)):
        p = classify_distribution("G", dict(zip(species, bits))).pattern
        counts[p] = counts.get(p, 0) + 1
    branches = ["Chlamydomonas_loss", "Conserved", "Angiosperm_specific", "Land_plant_innovation",
                "Bryophyte_origin", "Sporadic", "Unknown"]
    return {"combinations": 2 ** len(species), "counts": counts,
            "never_returned": [b for b in branches if b not in counts]}


def benchmarks():
    return {"ef_hand": ef_hand(), "tm_min_length": tm_min_length(), "examples_tm": examples_tm(),
            "phylo": phylo_reachability(), "ligands": repo_data.ligands()}
