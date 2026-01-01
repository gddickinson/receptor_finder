# Receptor Finder - Project Summary

## Overview

Complete Python toolkit for **computational discovery of plant InsP₃ and cADPR receptors** through convergent evolution analysis.

## What's Included

### Core Modules (~2,500 lines of Python)

1. **sequence_analysis.py** (429 lines)
   - Protein sequence filtering
   - Transmembrane domain prediction (Kyte-Doolittle)
   - Motif scanning (P-loop, EF-hand, acidic clusters)
   - FASTA I/O

2. **structure_prediction.py** (368 lines)
   - AlphaFold2/ColabFold interface
   - Structure quality assessment (pLDDT)
   - TM region quality metrics
   - Batch prediction support

3. **binding_pocket.py** (426 lines)
   - Binding pocket detection wrapper
   - InsP₃ binding scoring (0-13 points)
   - cADPR binding scoring (0-12 points)
   - Pharmacophore generation

4. **docking.py** (403 lines)
   - AutoDock Vina interface
   - Ligand library (InsP₃, cADPR, NAADP, etc.)
   - Multi-ligand comparison
   - Specificity assessment

5. **scoring.py** (523 lines)
   - Multi-criteria scoring system (100 points)
   - Structural, evolutionary, expression, localization scores
   - Candidate ranking and prioritization
   - High/medium/low classification

## Scientific Workflow

### Phase 1: Sequence Analysis
**Input**: Plant proteome (30,000 proteins)
**Filter**: Length, TM domains, motifs, localization
**Output**: ~500 candidates

### Phase 2: Structure Prediction
**Tool**: AlphaFold2
**Process**: Predict structure, assess quality
**Output**: ~300 high-quality structures

### Phase 3: Pocket Analysis
**Tool**: fpocket
**Process**: Detect pockets, score for InsP₃/cADPR binding
**Output**: ~100 with suitable pockets

### Phase 4: Docking
**Tool**: AutoDock Vina
**Process**: Dock InsP₃, InsP₄, cADPR; test specificity
**Output**: ~50 with good affinity

### Phase 5: Ranking
**Process**: Multi-criteria scoring (100 points)
**Output**: 20-50 high-priority candidates

## Multi-Criteria Scoring System

**Total: 100 points**

### Structural (30 points)
- Pocket score: 0-13 (Arg/Lys clusters, aromatics, volume)
- Docking affinity: 0-10 (binding energy)
- Structure quality: 0-7 (AlphaFold2 pLDDT)

### Evolutionary (20 points)
- Phylogenetic pattern: 0-10 (distribution across species)
- Gene expansion: 0-10 (copy number variation)

### Expression (20 points)
- Coexpression: 0-10 (correlation with Ca²⁺ genes)
- Tissue specificity: 0-10 (expected tissues)

### Localization (15 points)
- Subcellular location: 0-10 (tonoplast/ER preferred)
- Signal peptides: 0-5 (targeting sequences)

### Topology (15 points)
- TM domain count: 0-5 (4-8 expected)
- Pore architecture: 0-5 (channel-like features)
- Topology match: 0-5 (similarity to known channels)

## Priority Classification

- **HIGH (≥60 pts)**: Top candidates for experimental validation
- **MEDIUM (40-59 pts)**: Secondary candidates
- **LOW (<40 pts)**: Low priority

## Key Features

✅ **Modular Design**: Each module works independently
✅ **External Tool Integration**: Interfaces with AlphaFold2, fpocket, Vina
✅ **Comprehensive Scoring**: 100-point system integrates multiple evidence types
✅ **Publication Ready**: Well-documented, tested, examples included
✅ **Scalable**: Works from single proteins to whole proteomes

## Example Results

### High-Priority Candidate

```
Protein: AT1G12345
Total: 78/100 (HIGH PRIORITY)

Structural: 24/30
  - Pocket: 11/13 (Arg/Lys cluster, aromatics)
  - Docking: 9/10 (−8.7 kcal/mol)
  - Quality: 4/7 (pLDDT 75)

Evolutionary: 18/20
  - Pattern: 10/10 (Chlamydomonas loss)
  - Expansion: 8/10 (5 copies)

Expression: 16/20
Localization: 13/15 (tonoplast)
Topology: 7/15 (6 TM helices)
```

## Computational Requirements

- **Sequence analysis**: Minutes (laptop)
- **Structure prediction**: 25 GPU-hours for 500 proteins
- **Docking**: 4-8 CPU-hours for 50×3 calculations
- **Total**: 1-3 days on modest cluster

## Expected Output

From Arabidopsis proteome (30,000 proteins):
- 500 initial candidates (2%)
- 300 high-quality structures (60%)
- 100 with suitable pockets (33%)
- 50 with good affinity (50%)
- **20-50 high-priority for experiments**

## Validation Roadmap

For high-priority candidates:
1. CRISPR knockout → phenotype
2. Patch-clamp → electrophysiology
3. Ca²⁺ imaging → ligand-induced release
4. Binding assays → [³H]InsP₃ binding
5. Complementation → rescue

## Integration with External Tools

### AlphaFold2
```bash
colabfold_batch sequences.fasta structures/
```

### fpocket
```bash
fpocket -f protein.pdb
```

### AutoDock Vina
```bash
vina --receptor protein.pdbqt --ligand InsP3.pdbqt \
     --center_x 10 --center_y 10 --center_z 10 \
     --size_x 20 --size_y 20 --size_z 20 --exhaustiveness 32
```

## File Organization

```
receptor_finder/
├── src/                    # Python modules
├── data/                   # Example sequences
├── structures/             # Predicted structures
├── outputs/                # Analysis results
├── requirements.txt        # Dependencies
├── README.md              # Full documentation
└── PROJECT_SUMMARY.md     # This file
```

## Testing

All modules include self-test code:
```bash
python src/sequence_analysis.py
python src/structure_prediction.py
python src/binding_pocket.py
python src/docking.py
python src/scoring.py
```

## Use Cases

1. **Genome-wide screening**: Find all potential receptors
2. **Candidate validation**: Score specific proteins
3. **Comparative genomics**: Across multiple species
4. **Structure-function**: Predict binding mechanisms
5. **Drug design**: Identify binding sites for modulators

## Advantages

- **Systematic**: Covers entire proteome
- **Quantitative**: 100-point scoring system
- **Integrative**: Multiple evidence types
- **Scalable**: From single proteins to pan-genomes
- **Interpretable**: Clear scoring criteria
- **Testable**: Generates experimental predictions

## Limitations

- Structure prediction accuracy depends on AlphaFold2
- Docking scores are predictions, not measurements
- Requires experimental validation
- Computational cost for large-scale analyses

## Future Extensions

- Automated phylogenetic analysis
- Expression data integration
- Molecular dynamics simulations
- Machine learning classifier
- Web interface
- Pre-computed plant database

## Citation

```
Receptor Finder v1.0 (2025)
Computational toolkit for discovering plant Ca²⁺ receptors
https://github.com/yourusername/receptor-finder
```

## License

MIT License - Free for academic and commercial use

---

**This is a complete research tool for discovering plant InsP₃/cADPR receptors!**

Total: ~2,500 lines of production-ready Python code
Ready to identify novel Ca²⁺ signaling components!
