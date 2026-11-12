"""FastAPI web app for Receptor Finder: serves the front end in static/ and the JSON API under /api.

Run with `python run_web.py` (http://127.0.0.1:8206, local only). `scripts/build_static.py` exports every path in
STATIC_ENDPOINTS to JSON for the GitHub Pages site, so keep that list in step with the GET routes below. The POST routes
run the repo's Python code; on the static site the same answers come from the JS ports (seq_engine.js, score_engine.js),
which are node-tested against this Python. No route calls an external database or service.
"""
from pathlib import Path
from typing import Dict, List, Optional

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from web import audit, repo_data
from sequence_analysis import MotifScanner, SequenceDatabase, TransmembraneDomainPredictor  # via repo_data's path
import scoring as sc

STATIC = Path(__file__).resolve().parent / "static"
MAX_FASTA = 2_000_000  # characters; the local app is for small sets, not whole proteomes in one request

app = FastAPI(title="Receptor Finder")
app.mount("/static", StaticFiles(directory=STATIC), name="static")


@app.get("/", include_in_schema=False)
def index():
    return FileResponse(STATIC / "index.html")


@app.get("/api/info")
def info():
    m = repo_data.motifs()
    return {"title": "Receptor Finder", "description": repo_data.project_description(),
            "n_motifs": len(m["motifs"]), "n_ligands": len(repo_data.ligands()),
            "n_examples": len(repo_data.examples()["records"]), "n_species": len(repo_data.PLANT_SPECIES),
            "n_score_fields": len(sc.CandidateScore.__dataclass_fields__) - 1,
            "tm_defaults": m["tm_defaults"], "filter_defaults": m["filter_defaults"]}


@app.get("/api/motifs")
def motifs():
    return repo_data.motifs()


@app.get("/api/scoring")
def scoring_tables():
    return repo_data.scoring_tables()


@app.get("/api/examples")
def examples():
    return repo_data.examples()


@app.get("/api/files")
def files():
    return repo_data.committed_files()


@app.get("/api/audit")
def audit_results():
    return audit.benchmarks()


class AnalyseIn(BaseModel):
    fasta: str
    window: int = 19
    threshold: float = 1.6
    min_length: int = 200
    max_length: int = 2000
    min_tm: int = 2
    max_tm: int = 10


@app.post("/api/analyse")
def analyse(body: AnalyseIn):
    """Same shape as seq_engine.analyse(): TM domains, motif hits and filter verdicts per sequence (Python)."""
    if len(body.fasta) > MAX_FASTA:
        raise HTTPException(413, "FASTA too large for one request")
    if not 5 <= body.window <= 60:
        raise HTTPException(422, "window must be 5-60")
    tmp = SequenceDatabase()
    _load_text(tmp, body.fasta)
    out = []
    for s in tmp.sequences:
        tm = TransmembraneDomainPredictor.predict_tm_domains(s.sequence, body.window, body.threshold)
        length_ok = body.min_length <= s.length <= body.max_length
        tm_ok = body.min_tm <= len(tm) <= body.max_tm
        out.append({"id": s.sequence_id, "description": s.description, "length": s.length,
                    "tm": [list(t) for t in tm], "n_tm": len(tm),
                    "motifs": {k: MotifScanner.scan_sequence(s.sequence, p) for k, p in MotifScanner.MOTIFS.items()},
                    "length_ok": length_ok, "tm_ok": tm_ok, "passes": length_ok and tm_ok})
    return {"engine": "python", "results": out}


def _load_text(db, text):
    """SequenceDatabase.load_from_fasta reads a path; feed it the posted text through a temporary file."""
    import os
    import tempfile
    fd, path = tempfile.mkstemp(suffix=".fasta")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(text)
        db.load_from_fasta(path)
    finally:
        os.unlink(path)


class ScoreIn(BaseModel):
    affinity: float = -8.0
    mean_plddt: float = 80.0
    tm_region_plddt: Optional[float] = None
    pocket_raw: float = 8.0
    pattern: str = "Unknown"
    copies: int = 1
    ancestral: int = 1
    correlation: float = 0.0
    tissues: List[str] = []
    expected_tissues: List[str] = []
    localization: str = "cytosol"
    has_signal: bool = False
    signal_type: Optional[str] = None
    n_tm: int = 0
    pore_loop: bool = False
    selectivity_filter: bool = False
    topology_match: float = 0.0


@app.post("/api/score")
def score(b: ScoreIn):
    """Build a CandidateScore from raw evidence with the scorers in src/scoring.py."""
    c = sc.CandidateScore(
        protein_id="query",
        pocket_score=sc.StructuralScorer.score_pocket(b.pocket_raw),
        docking_score=sc.StructuralScorer.score_docking(b.affinity),
        structure_quality=sc.StructuralScorer.score_structure_quality(b.mean_plddt, b.tm_region_plddt),
        phylo_pattern_score=sc.EvolutionaryScorer.score_phylogenetic_pattern(b.pattern),
        expansion_score=sc.EvolutionaryScorer.score_gene_expansion(b.copies, b.ancestral),
        coexpression_score=sc.ExpressionScorer.score_coexpression(b.correlation),
        tissue_specificity=sc.ExpressionScorer.score_tissue_specificity(b.tissues, b.expected_tissues),
        localization_score=sc.LocalizationScorer.score_localization(b.localization),
        signal_peptide_score=sc.LocalizationScorer.score_signal_peptide(b.has_signal, b.signal_type),
        tm_domain_score=sc.TopologyScorer.score_tm_domains(b.n_tm),
        pore_architecture=sc.TopologyScorer.score_pore_architecture(b.pore_loop, b.selectivity_filter),
        topology_match=max(0.0, min(5.0, b.topology_match)),
    )
    fields: Dict[str, float] = {k: float(getattr(c, k)) for k in sc.CandidateScore.__dataclass_fields__
                                if k != "protein_id"}
    return {"engine": "python", "fields": fields, "total": c.total_score, "priority": c.priority_class}


# GET paths (without the /api prefix) written to the static site as api/<path>.json.
STATIC_ENDPOINTS = ["/info", "/motifs", "/scoring", "/examples", "/files", "/audit"]
