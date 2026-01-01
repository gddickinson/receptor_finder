# Receptor Finder - GUI Quick Start

## Installation (5 minutes)

1. **Install dependencies**:
```bash
cd receptor_finder
pip install -r requirements.txt
```

2. **Verify installation**:
```bash
python -c "import numpy, pandas, matplotlib, seaborn; print('Dependencies OK!')"
```

## Launch GUI

```bash
python launch_gui.py
```

## Quick Demo (10 minutes)

### Step 1: Load Example Data
- Click **File** → **Load Example Data**
- You'll see 3 example protein sequences loaded

### Step 2: Sequence Analysis (Tab 1)
- View sequences in the right panel
- Click **Analyze Motifs** to see binding motifs
- Check database statistics in the left panel

### Step 3: Structure Prediction (Tab 2)
- Click **Simulate Predictions (Demo)**
- This creates mock AlphaFold2 predictions
- View quality scores (pLDDT)

### Step 4: Pocket Analysis (Tab 3)
- Click **Detect Pockets (Demo)**
- Click **Score for InsP₃ Binding**
- Click **Score for cADPR Binding**
- View pocket statistics

### Step 5: Molecular Docking (Tab 4)
- Select ligands: InsP₃ and cADPR (default)
- Click **Run Docking (Demo)**
- Click **Analyze Specificity** to see InsP₃ vs InsP₄
- View binding affinities

### Step 6: Scoring & Ranking (Tab 5)
- Click **Score All Candidates**
- Click **Rank by Total Score**
- View high-priority candidates (≥60 points)
- Click **Export Rankings (CSV)** to save

### Step 7: Visualizations (Tab 6)
- Click **Score Distribution** - see score histogram
- Click **Score Heatmap** - see component breakdown
- Click **Save All Plots** - export to `outputs/`

## Understanding Results

### Priority Classes

- **HIGH (≥60 pts)**: Top candidates for experimental validation
- **MEDIUM (40-59 pts)**: Secondary candidates
- **LOW (<40 pts)**: Low priority

### Score Components (100 points total)

**Structural (30 pts)**
- Pocket binding: 0-13
- Docking affinity: 0-10
- Structure quality: 0-7

**Evolutionary (20 pts)**
- Phylogenetic pattern: 0-10
- Gene expansion: 0-10

**Expression (20 pts)**
- Coexpression: 0-10
- Tissue specificity: 0-10

**Localization (15 pts)**
- Subcellular location: 0-10
- Signal peptides: 0-5

**Topology (15 pts)**
- TM domains: 0-5
- Pore architecture: 0-5
- Channel similarity: 0-5

## Working with Real Data

### Load Your Sequences

1. Prepare FASTA file with protein sequences
2. **File** → **Load Sequences (FASTA)**
3. Apply filters:
   - Length: 200-2000 aa
   - TM domains: 2-10
4. Click **Apply Filters**

### Generate AlphaFold2 Commands

1. Load sequences
2. Tab 2: **Generate ColabFold Command**
3. Copy command to terminal/HPC
4. Run: `colabfold_batch sequences.fasta structures/`

### Load Predicted Structures

After AlphaFold2 predictions:
1. Tab 2: **Load Predicted Structures**
2. Select directory with PDB files
3. Filter by quality (pLDDT > 70)

### Real Pocket Detection

Use fpocket:
```bash
fpocket -f protein.pdb
```

### Real Molecular Docking

Use AutoDock Vina:
```bash
vina --receptor protein.pdbqt --ligand InsP3.pdbqt \
     --center_x 10 --center_y 10 --center_z 10 \
     --size_x 20 --size_y 20 --size_z 20 \
     --exhaustiveness 32
```

## Output Files

All results saved to `outputs/` directory:

- **Plots**: PNG images of analyses
- **Rankings**: CSV file with candidate scores
- **Structures**: PDB files (if generated)
- **Statistics**: Summary reports

## Tips

1. **Start with demo**: Use example data to learn workflow
2. **Filter aggressively**: Start with strict criteria
3. **Focus on high-priority**: Candidates with score ≥60
4. **Export frequently**: Save results at each stage
5. **Use visualizations**: Plots help identify patterns

## Troubleshooting

**GUI won't launch**:
```bash
python -m tkinter  # Test if tkinter works
pip install -r requirements.txt
```

**Import errors**:
```bash
pip install --upgrade numpy pandas matplotlib seaborn biopython
```

**No example data**:
- The example data is built-in
- If it fails, check `src/sequence_analysis.py`

## Expected Results

From Arabidopsis proteome (30,000 proteins):

1. After filtering: ~500 candidates
2. After structure prediction: ~300 high-quality
3. After pocket analysis: ~100 with suitable pockets
4. After docking: ~50 with good affinity
5. **Final: 20-50 high-priority candidates**

## Next Steps

1. **Experimental Validation**:
   - CRISPR knockout
   - Patch-clamp electrophysiology
   - Ca²⁺ imaging
   - Binding assays

2. **Refinement**:
   - Add phylogenetic data
   - Integrate expression data
   - Include tissue localization

3. **Scale Up**:
   - Analyze multiple species
   - Compare across plant families
   - Build candidate database

## Getting Help

- Check README.md for detailed documentation
- Run module tests: `python src/module_name.py`
- Review PROJECT_SUMMARY.md for workflow overview

---

**Ready to discover novel plant Ca²⁺ receptors!** 🔬🌱
