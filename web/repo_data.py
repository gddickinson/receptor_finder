"""Read-only facts for the web interface, taken from the repo's own code and committed files (nothing invented here).

The Python modules live in src/ and import each other as top-level modules (as in tests/), so src/ goes on sys.path.
"""
import hashlib
import inspect
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from docking import LigandLibrary  # noqa: E402
from phylogenetics import PLANT_SPECIES  # noqa: E402
from scoring import EvolutionaryScorer, LocalizationScorer  # noqa: E402
from sequence_analysis import MotifScanner, SequenceDatabase, TransmembraneDomainPredictor  # noqa: E402
from utils import read_fasta  # noqa: E402

EXAMPLE_FASTA = ROOT / "data" / "example_sequences.fasta"


def _defaults(fn):
    return {k: p.default for k, p in inspect.signature(fn).parameters.items() if p.default is not inspect.Parameter.empty}


def project_description() -> str:
    for line in (ROOT / "pyproject.toml").read_text().splitlines():
        if line.startswith("description"):
            return line.split("=", 1)[1].strip().strip('"')
    return ""


def motifs():
    return {
        "motifs": MotifScanner.MOTIFS,
        "kyte_doolittle": TransmembraneDomainPredictor.HYDROPHOBICITY,
        "tm_defaults": _defaults(TransmembraneDomainPredictor.predict_tm_domains),
        "filter_defaults": _defaults(SequenceDatabase.filter_database),
    }


def scoring_tables():
    return {
        "pattern_scores": EvolutionaryScorer.PATTERN_SCORES,
        "localization_scores": LocalizationScorer.LOCALIZATION_SCORES,
        "signal_type_scores": LocalizationScorer.SIGNAL_TYPE_SCORES,
        "plant_species": PLANT_SPECIES,
    }


# Phosphate groups each ligand has by definition (InsP3 = inositol 1,4,5-trisphosphate, InsP4 = a tetrakisphosphate,
# cADPR / NAD+ / 8-NH2-cADPR = one pyrophosphate, NAADP = pyrophosphate + 2'-phosphate). Counted against the SMILES.
EXPECTED_PHOSPHATES = {"InsP3": 3, "InsP4": 4, "cADPR": 2, "NAADP": 3, "NAD": 2, "8-NH2-cADPR": 2}
CYTOSINE = "N6C=CC(=NC6=O)N"


def ligands():
    out = []
    for name, smi in LigandLibrary.LIGANDS.items():
        p = smi.count("P")
        notes = []
        if p != EXPECTED_PHOSPHATES[name]:
            notes.append(f"{p} phosphorus atoms in the SMILES; {name} has {EXPECTED_PHOSPHATES[name]} phosphates")
        if name.startswith("Ins") and "O1" in smi:
            notes.append("ring closes through an oxygen (a pyranose sugar); myo-inositol is a carbocycle")
        if CYTOSINE in smi:
            notes.append("contains a cytosine base (" + CYTOSINE + "); cADPR has adenine and no cytosine")
        out.append({"name": name, "smiles": smi, "phosphorus": p, "expected_phosphates": EXPECTED_PHOSPHATES[name],
                    "notes": notes})
    return out


def examples():
    seqs = read_fasta(str(EXAMPLE_FASTA))
    text = EXAMPLE_FASTA.read_text()
    descs = {ln[1:].split(None, 1)[0]: (ln[1:].split(None, 1) + [""])[1] for ln in text.splitlines() if ln.startswith(">")}
    return {"source": "data/example_sequences.fasta", "fasta": text,
            "records": [{"id": k, "description": descs.get(k, ""), "length": len(v)} for k, v in seqs.items()]}


# What each committed data/output file is, from reading the code that writes it (see the commit history and the
# placeholder notes in src/). A test checks every committed file under data/ and outputs/ has an entry.
FILE_STATUS = {
    "data/example_sequences.fasta": ("example input",
        "Written by create_example_database() in sequence_analysis.py, which calls them 'hypothetical candidates'. "
        "AT1G12345 and AT3G34567 are low-complexity designed sequences; most of AT2G23456 reproduces the calmodulin "
        "EF-hand region followed by a designed tail. The IDs are placeholders, not Arabidopsis loci."),
    "outputs/sequences_for_prediction.fasta": ("example input",
        "Byte-for-byte the same three sequences as data/example_sequences.fasta."),
    "outputs/structures/AT1G12345.pdb": ("placeholder", "Mock structure from AlphaFold2Interface.predict_structure: "
        "CA atoms on a straight line (3.8 A steps) and simulated pLDDT (sine + random noise). Not a prediction."),
    "outputs/structures/AT2G23456.pdb": ("placeholder", "As above: mock structure, simulated pLDDT."),
    "outputs/structures/AT3G34567.pdb": ("placeholder", "As above: mock structure, simulated pLDDT."),
    "outputs/score_distribution.png": ("figure, untraceable",
        "Drawn by visualization.py; the scores behind it are not committed, so it is not shown here."),
    "outputs/score_heatmap.png": ("figure, untraceable", "As above."),
    "outputs/docking_comparison.png": ("figure, untraceable",
        "As above. Docking in this repo is a placeholder (random affinities), so no docking result is real."),
    "outputs/structure_quality.png": ("figure, untraceable", "As above; pLDDT values in this repo are simulated."),
}


def committed_files():
    out = []
    for rel, (kind, note) in FILE_STATUS.items():
        f = ROOT / rel
        if f.exists():
            data = f.read_bytes()
            out.append({"path": rel, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()[:16],
                        "kind": kind, "note": note})
    return out
