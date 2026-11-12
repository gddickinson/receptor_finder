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

**Follow-up the same day (lead, under George's "fix any other problems" rule).** The five open items were fixed, each
in its own "Fix ..." commit with tests and sources in THIRD_PARTY.md:
- EF_hand is now PROSITE PS00018, copied from the entry. It finds all 4 calmodulin loops and rejects their D1A mutants.
- TM prediction uses the standard Kyte & Doolittle (1982) rule (window 19, mean >= 1.6;
  doi:10.1016/0022-2836(82)90515-0). `min_windows` and `merge_gap` are configurable. On UniProt-annotated proteins it
  finds 13/16 helices (the old rule found 1/16) and nothing in calmodulin or myoglobin. Its one extra call is the
  glycophorin A signal peptide.
- `classify_distribution` follows docs/in_silico_calcium_projects.md section 3.1 and can return all 7 patterns. The
  old "Chlamydomonas_loss" went to families absent from the alga, the opposite of the documented definition.
- Ligand SMILES are PubChem isomeric SMILES (CIDs 439456, 107758, 123847, 123952, 5892), fetched once on 2026-11-12.
  8-NH2-cADPR is derived from cADPR, because PubChem's only record for it (CID 135600320) is not cyclic.
- scoring.py was split into scoring.py + scorers.py, and gui_app.py into gui_app.py + 3 mixins (no behaviour change).

**Still for George:**
- Confirm the derived 8-NH2-cADPR structure.
- PROSITE is CC BY-NC-ND 4.0 while this repo is MIT; only the attributed pattern string is reproduced.
- Pattern, TM and docking-score results all changed with these fixes. Nothing downstream is committed, so no stored
  result is affected.

## 2026-11-12 (later) -- 8-NH2-cADPR, PROSITE licence, pocket specificity (Claude)
Both items left "for George" are now settled (details and sources in THIRD_PARTY.md):
- **8-NH2-cADPR confirmed.** The SMILES (unchanged) is the neutral form of ChEMBL CHEMBL3350988 / PubChem CID
  118719244 (Zhang et al. 2008, doi:10.1021/jm7010386). The InChIKey is identical once the zwitterion is neutralised,
  and so are all 8 stereocentres. Walseth & Lee 1993 (doi:10.1016/0167-4889(93)90199-Y) describe it as cADPR with
  the 8-proton replaced. ChEBI has no entry, and the Guide to Pharmacology needs an account.
  `tests/test_ligand_structures.py` checks the formula with H, the 5 rings, the macrocycle, the 8-amino carbon, and
  that removing the amino group gives cADPR exactly.
- **PROSITE removed.** PROSITE as a whole, patterns included, is CC BY-NC-ND 4.0 with a fee for commercial users,
  which is incompatible with MIT. EF_hand is now our own regex for the canonical 12-residue loop,
  `D.[DENQST].[DENQST]G.{5}[DE]`, built from the literature (InterPro IPR002048, Mazumder 2014, Lin 2021).
  - It finds all 4 calmodulin loops and 13/13 UniProt-annotated loops in 6 EF-hand proteins
    (tests/fixtures/uniprot_efhand_benchmark.json). PROSITE found 12/13.
  - It hits none of the 5 pseudo or non-binding EF-hands.
  - It is less specific on random sequence (about 23 vs 1 hits per Mb).
  - The PROSITE string is still in the git history (commit e756dc6).
- **`PocketCharacterizer.assess_specificity`** marked every cADPR call specific without a test. It now needs a NAD+
  score more than 2 points lower, the same rule as InsP3 vs InsP4 (optional `nad_score`).
- Unchanged by decision: `TopologyScorer.score_tm_domains` rewards 4-8 TM (channel-like) while the filter keeps 2-10.
  That is a wide filter followed by a narrower score, not a contradiction.
- Tests: pytest 93 passed; node 529 scoring cases and 39 sequences match Python; screenshots --static and --base all ok
  on 6 routes, with 0 px overflow at 390 px.
