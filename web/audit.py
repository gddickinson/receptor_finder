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
    """Scan all of AT2G23456, then check each known 12-residue loop start is among the hits. Decoy: each loop with its
    position-1 Asp mutated to Ala. Also records which calmodulin (UniProt P0DP23) residues AT2G23456 reproduces."""
    seq = read_fasta(str(repo_data.EXAMPLE_FASTA))["AT2G23456"]
    pat = MotifScanner.MOTIFS["EF_hand"]
    hits = MotifScanner.scan_sequence(seq, pat)
    loops = []
    for lp in EF_LOOPS:
        pos = seq.find(lp)
        mutant = seq[:pos] + "A" + seq[pos + 1:]
        loops.append({"loop": lp, "next": seq[pos + 12], "position": pos, "hit": pos in hits,
                      "d1a_mutant_hit": pos in MotifScanner.scan_sequence(mutant, pat)})
    return {"pattern": pat, "source": "canonical 12-residue EF-hand loop (THIRD_PARTY.md)", "loops": loops,
            "detected": sum(x["hit"] for x in loops), "decoys_detected": sum(x["d1a_mutant_hit"] for x in loops),
            "hits_in_AT2G23456": hits, "calmodulin_match": calmodulin_match(seq)}


def _efhand_bench():
    import json
    return json.loads((repo_data.ROOT / "tests" / "fixtures" / "uniprot_efhand_benchmark.json").read_text())


def calmodulin_match(seq):
    """Longest identical stretch between AT2G23456 and human calmodulin CALM1 (UniProt P0DP23), 1-based inclusive."""
    from difflib import SequenceMatcher
    cam = next(e for e in _efhand_bench()["entries"] if e["accession"] == "P0DP23")["sequence"]
    m = max(SequenceMatcher(None, seq, cam, autojunk=False).get_matching_blocks(), key=lambda b: b.size)
    return {"accession": "P0DP23", "length": m.size, "example_range": [m.a + 1, m.a + m.size],
            "uniprot_range": [m.b + 1, m.b + m.size]}


def ef_hand_benchmark():
    """EF_hand on UniProt-annotated EF-hand proteins (tests/fixtures/uniprot_efhand_benchmark.json, fetched once).
    An annotated canonical loop starts at s when residues s, s+2, s+4 and s+11 are all Ca2+ ligands in UniProt."""
    bench = _efhand_bench()
    pat = MotifScanner.MOTIFS["EF_hand"]
    rows = []
    for e in bench["entries"]:
        ca = set(e["ca_binding_residues"])
        loops = sorted(s for s in ca if {s, s + 2, s + 4, s + 11} <= ca)
        hits = [h + 1 for h in MotifScanner.scan_sequence(e["sequence"], pat)]
        rows.append({"accession": e["accession"], "id": e["id"], "name": e["name"], "length": len(e["sequence"]),
                     "ef_hand_domains": len(e["ef_hand_domains"]), "annotated_loops": loops, "hits": hits,
                     "found": len(set(hits) & set(loops)), "extra": len(set(hits) - set(loops))})
    return {"source": bench["source"], "retrieved": bench["retrieved"], "rows": rows,
            "found": sum(r["found"] for r in rows), "annotated": sum(len(r["annotated_loops"]) for r in rows),
            "extra": sum(r["extra"] for r in rows), "ef_hand_domains": sum(r["ef_hand_domains"] for r in rows)}


def tm_min_length(flank="K", core="L", flank_len=20):
    rows = []
    for n in range(10, 41):
        seq = flank * flank_len + core * n + flank * flank_len
        rows.append({"n": n, "domains": [list(t) for t in TransmembraneDomainPredictor.predict_tm_domains(seq)]})
    first = next((r["n"] for r in rows if r["domains"]), None)
    return {"construct": f"{flank}x{flank_len} + {core}xn + {flank}x{flank_len}", "rows": rows,
            "shortest_detected": first}


def tm_benchmark():
    """UniProt-annotated TM proteins and soluble decoys (tests/fixtures/uniprot_tm_benchmark.json, fetched once)."""
    import json
    bench = json.loads((repo_data.ROOT / "tests" / "fixtures" / "uniprot_tm_benchmark.json").read_text())
    rows = []
    for e in bench["entries"]:
        pred = TransmembraneDomainPredictor.predict_tm_domains(e["sequence"])
        ann = e["transmembrane"]
        hit = lambda s, t: any(a < t and s - 1 < b for a, b in pred)  # noqa: E731
        rows.append({"accession": e["accession"], "id": e["id"], "name": e["name"], "length": len(e["sequence"]),
                     "annotated": ann, "predicted": [list(p) for p in pred],
                     "found": sum(hit(s, t) for s, t in ann),
                     "extra": sum(not any(a < t and s - 1 < b for s, t in ann) for a, b in pred)})
    return {"source": bench["source"], "retrieved": bench["retrieved"], "rows": rows,
            "found": sum(r["found"] for r in rows), "annotated": sum(len(r["annotated"]) for r in rows),
            "extra": sum(r["extra"] for r in rows)}


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
    return {"ef_hand": ef_hand(), "ef_hand_benchmark": ef_hand_benchmark(), "tm_min_length": tm_min_length(), "tm_benchmark": tm_benchmark(),
            "examples_tm": examples_tm(),
            "phylo": phylo_reachability(), "ligands": repo_data.ligands()}
