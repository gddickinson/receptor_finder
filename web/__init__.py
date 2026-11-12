"""Web interface for Receptor Finder (FastAPI app + static front end).

Contents:
- app.py        FastAPI routes under /api and STATIC_ENDPOINTS for the GitHub Pages export
- repo_data.py  read-only facts from the repo: info, method, motifs, scoring tables, ligands, committed files
- audit.py      benchmarks of the shipped motif / TM / phylogeny rules, computed by the repo's own Python code
- static/       the ES-module front end (pages, JS engines seq_engine.js / score_engine.js and their node tests)
"""
