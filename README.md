# Receptor Finder

A computational toolkit for discovering candidate InsP₃ and cADPR receptors in plant genomes through convergent evolution analysis.

## Overview

This project implements **Project #1** from the computational plant biology research proposals: genome-wide discovery of candidate plant Ca²⁺ receptors using structural prediction, binding pocket analysis, molecular docking, and phylogenetic methods.

## Scientific Background

Plants lack canonical InsP₃ receptors (InsP₃Rs) and ryanodine receptors (RyRs) found in animals, yet show robust functional evidence for InsP₃-induced and cADPR-induced Ca²⁺ release. This suggests **convergent evolution** - plants evolved novel protein families to perform analogous functions.

This toolkit systematically searches for these "missing" receptors by:
1. Identifying structural features required for ligand binding
2. Screening plant proteomes for proteins with convergent motifs
3. Predicting 3D structures and binding pockets
4. Simulating ligand docking
5. Ranking candidates using multi-criteria scoring

## Features

### Core Capabilities

- **Sequence Analysis**: Protein filtering, TM domain prediction, motif scanning
- **Structure Prediction**: AlphaFold2 interface, quality assessment
- **Binding Pocket Detection**: Pocket identification and characterization
- **Molecular Docking**: AutoDock Vina interface, multi-ligand comparison
- **Multi-Criteria Scoring**: Integrated ranking system (100-point scale)
- **Phylogenetic Analysis**: Distribution patterns, gene family expansion
- **Comprehensive Visualization**: Heatmaps, score distributions, candidate comparison

### Ligand Library

Built-in ligands for docking:
- InsP₃ (inositol 1,4,5-trisphosphate)
- InsP₄ (for specificity testing)
- cADPR (cyclic ADP-ribose)
- NAADP
- NAD⁺ (control)
- 8-NH₂-cADPR (antagonist)

## Installation

### Requirements

- Python 3.7+
- 500 MB disk space
- 2 GB RAM (recommended)

### Setup

```bash
cd receptor_finder
pip install -r requirements.txt
```

### Dependencies

Core packages (required):
- numpy, scipy, pandas
- matplotlib, seaborn
- biopython

External tools (optional but recommended):
- **AlphaFold2/ColabFold**: Structure prediction
- **fpocket**: Binding pocket detection
- **AutoDock Vina**: Molecular docking
- **HMMER**: HMM-based searches
- **BLAST+**: Sequence homology

## Quick Start

### Test Modules

```bash
# Test sequence analysis
python src/sequence_analysis.py

# Test structure prediction
python src/structure_prediction.py

# Test binding pocket analysis
python src/binding_pocket.py

# Test molecular docking
python src/docking.py

# Test scoring system
python src/scoring.py
```

### Example Workflow

```python
from src.sequence_analysis import SequenceDatabase, SequenceFilter
from src.structure_prediction import AlphaFold2Interface
from src.binding_pocket import PocketDetector, PocketScorer
from src.docking import AutoDockVinaInterface
from src.scoring import CandidateRanker, CandidateScore

# 1. Load and filter sequences
db = SequenceDatabase()
db.load_from_fasta("arabidopsis_proteome.fasta", species="Arabidopsis")
filtered = db.filter_database(min_tm=4, max_tm=8)

# 2. Predict structures
af2 = AlphaFold2Interface()
for seq in filtered.sequences:
    structure = af2.predict_structure(seq.sequence, seq.sequence_id)

# 3. Detect binding pockets
detector = PocketDetector()
scorer = PocketScorer()
for seq in filtered.sequences:
    pockets = detector.detect_pockets(f"{seq.sequence_id}.pdb", seq.sequence_id)
    for pocket in pockets:
        insp3_score = scorer.score_insp3_binding(pocket)

# 4. Rank candidates
ranker = CandidateRanker()
# ... add candidates with scores
ranked = ranker.rank_candidates()
```

## Project Structure

```
receptor_finder/
├── src/                          # Source code
│   ├── sequence_analysis.py      # Sequence filtering and analysis
│   ├── structure_prediction.py   # AlphaFold2 interface
│   ├── binding_pocket.py         # Pocket detection and scoring
│   ├── docking.py                # Molecular docking
│   ├── scoring.py                # Multi-criteria ranking
│   ├── phylogenetics.py          # (Future: phylogenetic analysis)
│   ├── visualization.py          # (Future: plotting tools)
│   └── utils.py                  # Utility functions
├── data/                         # Data files
│   └── example_sequences.fasta   # Example protein sequences
├── structures/                   # Predicted structures
├── outputs/                      # Analysis results
├── requirements.txt              # Python dependencies
└── README.md                     # This file
```

## Workflow Overview

### Phase 1: Sequence Analysis and Filtering

**Input**: Plant proteome (FASTA format)

**Steps**:
1. Filter by length (200-2000 aa)
2. Predict transmembrane domains (require 2-10 TM)
3. Scan for binding motifs (Arg/Lys clusters, aromatic residues)
4. Predict subcellular localization (tonoplast/ER preferred)

**Output**: ~500 candidate proteins

### Phase 2: Structure Prediction

**Tool**: AlphaFold2 (via ColabFold)

**For each candidate**:
1. Generate 5 models
2. Select highest pLDDT
3. Assess quality (mean pLDDT, TM region quality)
4. Filter: Keep only pLDDT > 70 in TM regions

**Output**: High-quality structures (~300)

### Phase 3: Binding Pocket Analysis

**Tool**: fpocket or P2Rank

**For each structure**:
1. Detect pockets (volume > 200 ų)
2. Characterize chemistry (charge, hydrophobicity)
3. Count key residues (basic, aromatic, acidic)
4. Score for InsP₃ binding (0-13 points)
5. Score for cADPR binding (0-12 points)

**Output**: Pockets with scores

### Phase 4: Molecular Docking

**Tool**: AutoDock Vina

**For top pockets**:
1. Prepare receptor (PDBQT format)
2. Prepare ligands (InsP₃, InsP₄, cADPR)
3. Dock with exhaustiveness=32
4. Extract binding affinity
5. Test specificity (InsP₃ vs InsP₄)

**Criteria**: Affinity < -7 kcal/mol AND InsP₃ > InsP₄

**Output**: Docking scores

### Phase 5: Multi-Criteria Ranking

**Scoring system** (100 points total):

- **Structural** (30 pts)
  - Pocket score: 0-13
  - Docking affinity: 0-10
  - Structure quality: 0-7

- **Evolutionary** (20 pts)
  - Phylogenetic pattern: 0-10
  - Gene expansion: 0-10

- **Expression** (20 pts)
  - Coexpression with Ca²⁺ genes: 0-10
  - Tissue specificity: 0-10

- **Localization** (15 pts)
  - Predicted localization: 0-10
  - Signal peptides: 0-5

- **Topology** (15 pts)
  - TM domain count: 0-5
  - Pore architecture: 0-5
  - Similarity to channels: 0-5

**Priority classes**:
- HIGH: ≥60 points (top 20-50 candidates)
- MEDIUM: 40-59 points
- LOW: <40 points

**Output**: Ranked candidate list

## Scoring Examples

### High-Priority Candidate

```
Protein: AT1G12345
Total Score: 78/100 (HIGH PRIORITY)

Structural (24/30):
  - Pocket score: 11/13 (Arg/Lys cluster, aromatic residues)
  - Docking: 9/10 (−8.7 kcal/mol for InsP₃)
  - Quality: 4/7 (pLDDT 75)

Evolutionary (18/20):
  - Phylo pattern: 10/10 (Chlamydomonas loss)
  - Expansion: 8/10 (5 copies vs 1 ancestral)

Expression (16/20):
  - Coexpression: 8/10 (r=0.72 with TPC1, CPK5)
  - Tissue: 8/10 (guard cells, root hairs)

Localization (13/15):
  - Location: 10/10 (tonoplast)
  - Signal: 3/5 (ER signal peptide)

Topology (7/15):
  - TM domains: 5/5 (6 TM helices)
  - Pore: 2/5 (no clear re-entrant loop)
```

### Medium-Priority Candidate

```
Protein: AT2G23456
Total Score: 52/100 (MEDIUM PRIORITY)

Structural (18/30):
  - Pocket score: 8/13
  - Docking: 7/10 (−7.8 kcal/mol)
  - Quality: 3/7 (pLDDT 68)

Evolutionary (10/20):
  - Phylo pattern: 6/10 (Angiosperm-specific)
  - Expansion: 4/10 (2 copies)

... (lower scores in other categories)
```

## Output Files

The toolkit generates:

- **sequences/filtered_candidates.fasta**: Filtered protein sequences
- **structures/*.pdb**: Predicted 3D structures
- **pockets/pocket_analysis.csv**: Binding pocket characteristics
- **docking/docking_results.csv**: Ligand binding affinities
- **rankings/candidate_rankings.csv**: Final ranked list

## Integration with External Tools

### AlphaFold2 (ColabFold)

```bash
# Batch prediction
colabfold_batch sequences.fasta structures/ --num-models 1 --num-recycle 3
```

### fpocket

```bash
# Pocket detection
fpocket -f protein.pdb
```

### AutoDock Vina

```bash
# Prepare receptor
prepare_receptor4.py -r protein.pdb -o protein.pdbqt

# Dock ligand
vina --receptor protein.pdbqt --ligand InsP3.pdbqt \
     --center_x 10 --center_y 10 --center_z 10 \
     --size_x 20 --size_y 20 --size_z 20 \
     --out docking_result.pdbqt --exhaustiveness 32
```

## Expected Results

Running the complete pipeline on Arabidopsis proteome (~30,000 proteins):

1. **Initial filtering**: ~500 candidates (2% of proteome)
2. **After structure prediction**: ~300 high-quality (60%)
3. **After pocket analysis**: ~100 with suitable pockets (33%)
4. **After docking**: ~50 with good affinity (50%)
5. **High-priority**: 20-50 candidates for experimental validation

## Experimental Validation Roadmap

For high-priority candidates:

1. **CRISPR knockout**: Loss-of-function phenotypes
2. **Patch-clamp**: Electrophysiology on isolated vacuoles
3. **Ca²⁺ imaging**: InsP₃/cADPR-induced Ca²⁺ release
4. **Binding assays**: [³H]InsP₃ or [³²P]cADPR binding
5. **Complementation**: Rescue of phenotypes

## Computational Resources

- **Sequence analysis**: Minutes (laptop)
- **Structure prediction**: Hours-Days (GPU recommended)
  - AlphaFold2: ~3 min/protein (GPU)
  - 500 proteins: ~25 GPU-hours
- **Docking**: Hours (CPU sufficient)
  - Vina: ~5-10 min/ligand
  - 50 proteins × 3 ligands: ~4-8 hours
- **Total**: 1-3 days on modest computing cluster

## Limitations and Considerations

1. **Structure prediction accuracy**: AlphaFold2 performs well on individual domains but may struggle with multi-domain proteins

2. **Docking reliability**: Scoring functions are approximate; binding affinities are predictions, not measurements

3. **Phylogenetic inference**: Requires comprehensive taxon sampling

4. **False positives**: High computational scores don't guarantee function - experimental validation is essential

5. **Novel folds**: True convergent evolution may produce folds AlphaFold2 hasn't seen

## License

MIT License - Free for academic and commercial use

## Contact

For questions, issues, or collaboration:
- Open an issue on GitHub
- Check documentation
- Review example workflows

## Acknowledgments

Methods based on:
- AlphaFold2 (Jumper et al. 2021)
- fpocket (Le Guilloux et al. 2009)
- AutoDock Vina (Trott & Olson 2010)
- Plant Ca²⁺ signaling literature

## Future Development

Planned features:
- Automated phylogenetic analysis
- Expression data integration
- Batch processing pipeline
- Web interface
- Pre-computed results database

---

**This is a research tool for computational biology!**
