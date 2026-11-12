"""The JS-vs-Python fixture (tests/fixtures/engine_fixture.json) must match what the current Python code produces.

If this fails after a change to src/, rerun scripts/make_engine_fixture.py and the node tests. Numbers are compared with
a tolerance because libm differs in the last digit between macOS and Linux.
"""
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import make_engine_fixture  # noqa: E402


def _same(a, b, path="$"):
    if isinstance(a, bool) or isinstance(b, bool) or a is None or b is None or isinstance(a, str):
        assert a == b, path
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)):
        assert math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12), f"{path}: {a} vs {b}"
    elif isinstance(a, dict):
        assert a.keys() == b.keys(), path
        for k in a:
            _same(a[k], b[k], f"{path}.{k}")
    else:
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            _same(x, y, f"{path}[{i}]")


def test_fixture_is_current():
    saved = json.loads(make_engine_fixture.OUT.read_text())
    fresh = json.loads(json.dumps(make_engine_fixture.build()))
    _same(fresh, saved)
