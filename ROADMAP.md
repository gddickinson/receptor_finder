# Receptor Finder -- Roadmap

## Current State
A computational toolkit for discovering candidate InsP3 and cADPR receptors in plant genomes through convergent evolution analysis. Modules: `sequence_analysis.py` (proteome filtering, TM prediction, motif scanning), `structure_prediction.py` (AlphaFold2 interface), `binding_pocket.py` (pocket detection and scoring), `docking.py` (AutoDock Vina interface), `scoring.py` (multi-criteria 100-point ranking with dataclasses), `visualization.py`, `gui_app.py`, and `check_dependencies.py`. The `scoring.py` is well-designed with typed dataclasses and clear score categories. Has a GUI launcher (`launch_gui.py`).

## Short-term Improvements
- [x] Implement `src/phylogenetics.py` -- gene family expansion analysis, taxonomic distribution patterns
- [x] Flesh out `src/utils.py` with shared helpers (FASTA I/O, PDB parsing, file validation)
- [x] Add unit tests for `scoring.py` (verify score ranges, priority classification thresholds)
- [x] Add unit tests for `sequence_analysis.py` (TM prediction, motif scanning on known sequences)
- [ ] Add input validation in `docking.py` for missing PDB/PDBQT files and malformed ligand structures
- [ ] Add timeout handling for AlphaFold2 and Vina subprocess calls in `structure_prediction.py` and `docking.py`
- [x] Create `INTERFACE.md` for project navigation

## Feature Enhancements
- [ ] Add batch processing pipeline script that chains all 5 phases (filter -> predict -> pockets -> dock -> rank)
- [ ] Implement expression data integration in `scoring.py` -- pull from TAIR or Phytozome databases
- [ ] Add interactive candidate comparison view in `gui_app.py` with side-by-side pocket visualizations
- [ ] Implement P2Rank as alternative to fpocket in `binding_pocket.py`
- [ ] Add specificity testing: dock InsP3 vs InsP4 and require selectivity in `docking.py`
- [ ] Generate HTML report with embedded plots for each analysis run
- [ ] Add support for multiple species beyond Arabidopsis (rice, maize, Marchantia)

## Long-term Vision
- [ ] Integrate molecular dynamics simulation for top candidates (GROMACS or OpenMM)
- [ ] Build a pre-computed results database for common plant proteomes
- [ ] Add web interface for submitting custom proteomes and viewing results
- [ ] Implement machine learning re-ranking using features from all pipeline stages
- [ ] Cloud deployment for structure prediction (ColabFold API integration)
- [ ] Integrate with experimental tracking: link predictions to CRISPR knockout phenotypes

## Technical Debt
- [ ] `src/gui_app.py` likely handles too many concerns -- split into separate panels/dialogs
- [ ] `visualization.py` is listed as "future" in README but exists -- verify completeness
- [ ] `check_dependencies.py` should run automatically on import, not as a separate script
- [ ] External tool paths (fpocket, Vina, colabfold_batch) are likely hardcoded -- make configurable
- [x] No `pyproject.toml` or modern packaging -- add for pip installability
- [x] Missing `.gitignore` for `structures/`, `outputs/`, and large FASTA files
