"""Web interface tests: every exported endpoint answers, the static build is complete, and the JS engine tests pass.

Run: python -m pytest -q tests/test_web.py
"""
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from web.app import STATIC, STATIC_ENDPOINTS, app  # noqa: E402

client = TestClient(app)


def test_index_served():
    r = client.get("/")
    assert r.status_code == 200 and "static/main.js" in r.text


@pytest.mark.parametrize("path", STATIC_ENDPOINTS)
def test_endpoint(path):
    r = client.get("/api" + path)
    assert r.status_code == 200, r.text
    assert r.json() is not None


def test_static_build(tmp_path):
    sys.path.insert(0, str(ROOT / "scripts"))
    import build_static
    build_static.build(tmp_path / "site")
    assert (tmp_path / "site" / "index.html").read_text().count("APP_STATIC") == 1
    for path in STATIC_ENDPOINTS:
        assert (tmp_path / "site" / "api" / (path.strip("/") + ".json")).exists()


def test_static_assets_referenced_exist():
    html = (STATIC / "index.html").read_text()
    for name in ("style.css", "main.js"):
        assert f"static/{name}" in html and (STATIC / name).exists()


@pytest.mark.skipif(not shutil.which("node"), reason="node not installed")
def test_js_engine():
    tests = sorted(STATIC.glob("*.test.mjs"))
    for t in tests:
        r = subprocess.run(["node", str(t)], capture_output=True, text=True)
        assert r.returncode == 0, r.stdout + r.stderr


# --- project-specific API tests -----------------------------------------------------------------------------------
EXAMPLE = (ROOT / "data" / "example_sequences.fasta").read_text()


def test_analyse_matches_python_core():
    sys.path.insert(0, str(ROOT / "src"))
    from sequence_analysis import TransmembraneDomainPredictor
    r = client.post("/api/analyse", json={"fasta": EXAMPLE}).json()
    assert r["engine"] == "python" and [x["id"] for x in r["results"]] == ["AT1G12345", "AT2G23456", "AT3G34567"]
    from utils import read_fasta
    seqs = read_fasta(str(ROOT / "data" / "example_sequences.fasta"))
    for x in r["results"]:
        assert x["tm"] == [list(t) for t in TransmembraneDomainPredictor.predict_tm_domains(seqs[x["id"]])]
    assert not any(x["passes"] for x in r["results"])


def test_analyse_rejects_bad_window():
    assert client.post("/api/analyse", json={"fasta": EXAMPLE, "window": 1}).status_code == 422


def test_score_uses_repo_scorers():
    r = client.post("/api/score", json={"affinity": -8.0, "localization": "ER", "n_tm": 6,
                                        "pore_loop": True, "selectivity_filter": True}).json()
    f = r["fields"]
    assert f["docking_score"] == 8.5 and f["localization_score"] == 9.0
    assert f["tm_domain_score"] == 5.0 and f["pore_architecture"] == 5.0
    assert abs(r["total"] - sum(f.values())) < 1e-9 and r["priority"] in ("HIGH", "MEDIUM", "LOW")


def test_every_committed_data_file_is_described():
    from web.repo_data import FILE_STATUS
    out = subprocess.run(["git", "ls-files", "data", "outputs"], cwd=ROOT, capture_output=True, text=True).stdout.split()
    if out:  # inside a git checkout
        assert sorted(out) == sorted(FILE_STATUS), "describe new data/output files in web/repo_data.py FILE_STATUS"


def test_audit_numbers():
    a = client.get("/api/audit").json()
    seq = EXAMPLE.split(">AT2G23456")[1].split(">")[0].split("\n", 1)[1].replace("\n", "")
    assert all(l["position"] == seq.find(l["loop"]) >= 0 for l in a["ef_hand"]["loops"])
    assert a["ef_hand"]["detected"] == 4 and a["ef_hand"]["decoys_detected"] == 0
    assert a["tm_min_length"]["shortest_detected"] == 23
    assert set(a["phylo"]["never_returned"]) == {"Angiosperm_specific", "Land_plant_innovation", "Bryophyte_origin", "Unknown"}
    lig = {x["name"]: x for x in a["ligands"]}
    assert lig["InsP3"]["phosphorus"] == 4 and lig["NAD"]["notes"] == []


def test_no_network_in_web_code():
    """The web app never calls an external service: no HTTP client in Python, only relative fetch() in JS."""
    import re
    for f in (ROOT / "web").rglob("*.py"):
        assert not re.search(r"^\s*(import|from)\s+(requests|urllib|httpx|aiohttp|http\.client)", f.read_text(), re.M), f
    for f in (ROOT / "web" / "static").glob("*.js"):
        assert not re.search(r"fetch\(\s*['\"`]https?:", f.read_text()), f
