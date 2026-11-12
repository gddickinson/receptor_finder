# Receptor Finder -- Interface Map

## Top-level Files
- `launch_gui.py` -- Entry point; launches `ReceptorFinderGUI` Tkinter app
- `requirements.txt` -- Pinned pip dependencies
- `pyproject.toml` -- Modern packaging configuration
- `README.md` -- Project documentation
- `THIRD_PARTY.md` -- sources, retrieval dates and licences of PROSITE, UniProt and PubChem data used

## Package: `src/`

### Core Pipeline Modules
- `sequence_analysis.py`
  - `ProteinSequence` (dataclass): sequence + annotations
  - `TransmembraneDomainPredictor`: Kyte-Doolittle TM prediction (window 19, threshold 1.6; configurable)
  - `SequenceFilter`: filter by length, TM domains, motifs
  - `MotifScanner`: regex motif scanning (EF_hand = PROSITE PS00018)
  - `SequenceDatabase`: FASTA I/O and database operations

- `structure_prediction.py`
  - `PredictedStructure` (dataclass): PDB + pLDDT quality metrics
  - `AlphaFold2Interface`: predict / batch_predict (placeholder for ColabFold)
  - `StructureAnalyzer`: pLDDT parsing, region quality, TM detection
  - `StructureDatabase`: manage predicted structures

- `binding_pocket.py`
  - `BindingPocket` (dataclass): pocket geometry + chemistry
  - `PocketDetector`: fpocket wrapper / placeholder detection
  - `PocketScorer`: score_insp3_binding, score_cadpr_binding
  - `PocketCharacterizer`: pharmacophore, specificity assessment
  - `PocketDatabase`: collection management

- `docking.py`
  - `DockingResult` (dataclass): affinity + interactions
  - `LigandLibrary`: PubChem/ChEMBL SMILES for InsP3, InsP4, cADPR, NAADP, NAD+, 8-NH2-cADPR (`LIGAND_SOURCES`)
  - `AutoDockVinaInterface`: dock_ligand, prepare_receptor
  - `DockingAnalyzer`: compare_ligands, assess_specificity
  - `DockingDatabase`: result collection

- `scoring.py`
  - `CandidateScore` (dataclass): 100-point multi-criteria score
  - `CandidateRanker`: rank, filter, export to DataFrame/CSV
  - re-exports the scorers below
- `scorers.py`
  - `StructuralScorer`, `EvolutionaryScorer` (`PATTERN_SCORES`), `ExpressionScorer`,
    `LocalizationScorer` (`LOCALIZATION_SCORES`, `SIGNAL_TYPE_SCORES`), `TopologyScorer`: individual scoring functions

- `phylogenetics.py`
  - `GeneFamily` (dataclass): copy numbers across species
  - `PhylogeneticPattern` (dataclass): distribution classification
  - `classify_distribution()`: pattern detection (Chlamydomonas_loss, etc.)
  - `build_copy_number_table()`: species-by-family DataFrame

- `visualization.py`
  - `CandidatePlotter`: score distributions, heatmaps, docking comparisons

### Support Modules
- `utils.py`
  - `read_fasta()`, `write_fasta()`: FASTA I/O helpers
  - `parse_pdb_residues()`, `pdb_sequence()`: PDB parsing
  - `validate_fasta()`, `validate_pdb()`, `validate_pdbqt()`: file validation

- `check_dependencies.py` -- Verifies required packages are installed
- `gui_app.py` -- `ReceptorFinderGUI`: 6-tab Tkinter interface (state + `main()`), built from three mixins:
  - `gui_tabs.py` -- `GuiTabsMixin`: menu bar and tab layouts
  - `gui_actions_analysis.py` -- `GuiAnalysisMixin`: sequence, structure and pocket actions
  - `gui_actions_results.py` -- `GuiResultsMixin`: docking, scoring, plots, help dialogs

## Web interface: `web/`
- `run_web.py` (repo root) -- starts the FastAPI app on http://127.0.0.1:8206
- `web/app.py` -- routes under /api (`/info`, `/motifs`, `/scoring`, `/examples`, `/files`, `/audit`; POST `/analyse`,
  `/score` run the src/ code) and `STATIC_ENDPOINTS` for the static export. No route calls an external service.
- `web/repo_data.py` -- read-only facts from src/ and committed files (motifs, defaults, scoring tables, ligand checks,
  `FILE_STATUS` describing every committed data/output file)
- `web/audit.py` -- benchmarks of the shipped rules (EF-hand loops, shortest TM run, phylogeny reachability)
- `web/static/` -- ES-module front end: `main.js` (router), pages `home.js`, `method.js`, `sequence.js` (+ `seq_view.js`),
  `scoring.js`, `audit.js`, `about.js`; `static_api.js` answers api() on the static site
  - `seq_engine.js` -- JS port of sequence_analysis.py (KD profile in numpy's summation order, TM domains, motifs, FASTA, filter)
  - `score_engine.js` -- JS port of the scoring.py scorers, PocketScorer, classify_distribution, docking specificity
  - `*.test.mjs` -- node tests replaying `tests/fixtures/engine_fixture.json`
- `scripts/make_engine_fixture.py` -- writes the Python-vs-JS fixture; `scripts/build_static.py` -- exports the static
  site; `scripts/screenshots.py` -- browser smoke test / README images
- `.github/workflows/pages.yml` -- tests, builds and publishes the static site on push to main

## `tests/`
- `test_scoring.py` -- Score ranges, priority thresholds, ranking
- `test_sequence_analysis.py` -- TM prediction, motif scanning, FASTA I/O
- `test_utils.py` -- FASTA/PDB helpers and validation
- `test_regressions.py` -- regression tests for bugs fixed 2026-11-12 (TM ends, ranking field, placeholder PDB label,
  ER localization, docking score direction)
- `test_engine_fixture.py` -- the JS-vs-Python fixture is current
- `test_ligand_structures.py` -- 8-NH2-cADPR is cyclic, has the 8-amino group and the right formula
- `smiles_graph.py` -- minimal SMILES reader used by the ligand tests (formula with H, ring count, ring bonds)
- `test_web.py` -- API, static build, JS engine node tests, audit numbers, no network calls

## Pipeline Flow
1. `sequence_analysis` -- Filter proteome for TM proteins with relevant motifs
2. `structure_prediction` -- Predict structures via AlphaFold2
3. `binding_pocket` -- Detect and score pockets for InsP3/cADPR binding
4. `docking` -- Dock ligands and assess specificity
5. `scoring` -- Combine all evidence into 100-point ranking
6. `visualization` -- Plot results
