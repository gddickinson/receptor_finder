# Session log

## 2026-11-12 -- web interface, evidence audit, bug fixes (Claude)

**Web interface** (`web/`, `run_web.py`, port 8206; static site via `scripts/build_static.py` + `.github/workflows/pages.yml`).
Pages: Home, Method, Sequence filter, Scoring, Audit, About. The sequence filter and scoring run through the Python code
in the local app and through JS ports (`seq_engine.js`, `score_engine.js`) on the static site; the ports are checked
against the Python on 28 sequences and 529 scoring cases (`tests/fixtures/engine_fixture.json`, hydropathy profiles
bitwise equal by reproducing numpy's pairwise summation). No route or test calls an external service.

**What the repo holds (audit).** No proteome search results or candidate lists. `data/example_sequences.fasta` holds three
hypothetical examples (placeholder IDs); `outputs/structures/*.pdb` are mock structures (straight CA trace, simulated
pLDDT); the four `outputs/*.png` figures were drawn from data that is not committed, so the site does not show them.
Nothing here looks like an unpublished finding.

**Bugs fixed** (each with a regression test in `tests/test_regressions.py`, separate commits):
- TM domain ends ran one residue past the domain (C-terminal domains ended at len + 1).
- `CandidateRanker.rank_candidates` silently returned an unsorted list for fields not ending in `_score`.
- Placeholder PDBs were labelled "AlphaFold2 PREDICTION" (writer and the three committed files relabelled).
- `score_localization('ER')` scored 0 instead of 9 (case mismatch).
- `score_docking` slopes were reversed: -8.5 kcal/mol scored 10.75/10, weaker binders scored higher.

**Needs George's decision** (benchmarks on the Audit page; code unchanged):
- `MotifScanner.MOTIFS['EF_hand']` (`D.{2}[DNS]...[DNS].{2}[EQ]`) matches 0 of the 4 calmodulin EF-hand loops in
  AT2G23456. A canonical alternative is PROSITE PS00018.
- TM prediction needs 15 consecutive windows above threshold: the shortest poly-Leu stretch called is 23 residues, so
  real TM helices (~20 residues) are often missed; the "multi-pass" example gets 1 segment.
- `classify_distribution` can only return Chlamydomonas_loss, Conserved or Sporadic (the other branches are
  unreachable); a family absent everywhere is "Sporadic".
- `LigandLibrary` SMILES: InsP3/InsP4 are pyranose sugars with 4/5 phosphates, cADPR and 8-NH2-cADPR contain cytosine,
  NAADP has a nitro group and 2 phosphates.
- `src/scoring.py` (569 lines) and `src/gui_app.py` (1178) are over the 500-line budget (pre-existing).
