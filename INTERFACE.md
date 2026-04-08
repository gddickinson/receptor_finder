# Receptor Finder -- Interface Map

## Top-level Files
- `launch_gui.py` -- Entry point; launches `ReceptorFinderGUI` Tkinter app
- `requirements.txt` -- Pinned pip dependencies
- `pyproject.toml` -- Modern packaging configuration
- `README.md` -- Project documentation

## Package: `src/`

### Core Pipeline Modules
- `sequence_analysis.py`
  - `ProteinSequence` (dataclass): sequence + annotations
  - `TransmembraneDomainPredictor`: hydrophobicity-based TM prediction
  - `SequenceFilter`: filter by length, TM domains, motifs
  - `MotifScanner`: regex motif scanning for Ca2+ channel signatures
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
  - `LigandLibrary`: SMILES for InsP3, cADPR, InsP4, etc.
  - `AutoDockVinaInterface`: dock_ligand, prepare_receptor
  - `DockingAnalyzer`: compare_ligands, assess_specificity
  - `DockingDatabase`: result collection

- `scoring.py`
  - `CandidateScore` (dataclass): 100-point multi-criteria score
  - `StructuralScorer`, `EvolutionaryScorer`, `ExpressionScorer`,
    `LocalizationScorer`, `TopologyScorer`: individual scoring functions
  - `CandidateRanker`: rank, filter, export to DataFrame/CSV

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
- `gui_app.py` -- `ReceptorFinderGUI`: 5-tab Tkinter interface

## `tests/`
- `test_scoring.py` -- Score ranges, priority thresholds, ranking
- `test_sequence_analysis.py` -- TM prediction, motif scanning, FASTA I/O
- `test_utils.py` -- FASTA/PDB helpers and validation

## Pipeline Flow
1. `sequence_analysis` -- Filter proteome for TM proteins with relevant motifs
2. `structure_prediction` -- Predict structures via AlphaFold2
3. `binding_pocket` -- Detect and score pockets for InsP3/cADPR binding
4. `docking` -- Dock ligands and assess specificity
5. `scoring` -- Combine all evidence into 100-point ranking
6. `visualization` -- Plot results
