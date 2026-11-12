"""
Docking, scoring, plotting and help actions of the Receptor Finder GUI.

Mixin for gui_app.ReceptorFinderGUI (split out of gui_app.py to keep files
under 500 lines). Methods use the state set up in ReceptorFinderGUI.__init__.
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent))

from sequence_analysis import (SequenceDatabase, ProteinSequence, SequenceFilter,
                               TransmembraneDomainPredictor, MotifScanner)
from structure_prediction import AlphaFold2Interface, StructureDatabase, StructureAnalyzer
from binding_pocket import PocketDetector, PocketScorer, PocketDatabase
from docking import AutoDockVinaInterface, DockingDatabase, LigandLibrary, DockingAnalyzer
from scoring import (CandidateScore, CandidateRanker, StructuralScorer,
                    EvolutionaryScorer, LocalizationScorer, TopologyScorer)
from visualization import CandidatePlotter


class GuiResultsMixin:
    """Docking, scoring, plotting and help actions of the Receptor Finder GUI."""

    def run_docking(self):
        """Run molecular docking."""
        if len(self.pocket_db) == 0:
            messagebox.showwarning("Warning", "Detect pockets first")
            return
        
        # Get selected ligands
        selected_ligands = [lig for lig, var in self.ligand_vars.items() if var.get()]
        
        if not selected_ligands:
            messagebox.showwarning("Warning", "Select at least one ligand")
            return
        
        self.status_var.set(f"Running docking with {len(selected_ligands)} ligands...")
        self.root.update()
        
        vina = AutoDockVinaInterface()
        
        # Dock to first 5 proteins
        proteins = list(set(p.protein_id for p in list(self.pocket_db.pockets.values())[:5]))
        
        for protein_id in proteins:
            for ligand in selected_ligands:
                results = vina.dock_ligand(
                    f"{protein_id}.pdbqt",
                    f"{ligand}.pdbqt",
                    center=(10.0, 10.0, 10.0)
                )
                for result in results:
                    self.docking_db.add_result(result)
        
        self.status_var.set(f"Completed {len(self.docking_db)} docking calculations")
        self.update_docking_display()
        
        messagebox.showinfo("Complete", f"Docking complete: {len(self.docking_db)} results (demo mode)")
    
    def analyze_specificity(self):
        """Analyze InsP3 vs InsP4 specificity."""
        if len(self.docking_db) == 0:
            messagebox.showwarning("Warning", "Run docking first")
            return
        
        analyzer = DockingAnalyzer()
        
        # Find proteins with both InsP3 and InsP4 results
        results_text = "SPECIFICITY ANALYSIS:\n\n"
        
        proteins = set(r.protein_id for r in self.docking_db.results)
        
        for protein_id in list(proteins)[:5]:
            insp3_results = [r for r in self.docking_db.results 
                           if r.protein_id == protein_id and r.ligand_id == 'InsP3']
            insp4_results = [r for r in self.docking_db.results
                           if r.protein_id == protein_id and r.ligand_id == 'InsP4']
            
            if insp3_results and insp4_results:
                insp3_best = min(insp3_results, key=lambda x: x.binding_affinity)
                insp4_best = min(insp4_results, key=lambda x: x.binding_affinity)
                
                assessment = analyzer.assess_specificity(
                    insp3_best.binding_affinity,
                    insp4_best.binding_affinity
                )
                
                results_text += f"{protein_id}:\n"
                results_text += f"  InsP₃: {assessment['insp3_affinity']:.2f} kcal/mol\n"
                results_text += f"  InsP₄: {assessment['insp4_affinity']:.2f} kcal/mol\n"
                results_text += f"  Δ: {assessment['delta_affinity']:.2f} kcal/mol\n"
                results_text += f"  Specific: {assessment['is_specific']}\n\n"
        
        # Show in dialog
        dialog = tk.Toplevel(self.root)
        dialog.title("Specificity Analysis")
        dialog.geometry("500x600")
        
        text_widget = scrolledtext.ScrolledText(dialog, width=60, height=35)
        text_widget.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        text_widget.insert(tk.END, results_text)
        text_widget.config(state=tk.DISABLED)
    
    def compare_ligands(self):
        """Compare docking results across ligands."""
        if len(self.docking_db) == 0:
            messagebox.showwarning("Warning", "Run docking first")
            return
        
        analyzer = DockingAnalyzer()
        
        # Group by ligand
        results_dict = {}
        for result in self.docking_db.results:
            if result.ligand_id not in results_dict:
                results_dict[result.ligand_id] = []
            results_dict[result.ligand_id].append(result)
        
        comparison = analyzer.compare_ligands(results_dict)
        
        # Display
        dialog = tk.Toplevel(self.root)
        dialog.title("Ligand Comparison")
        dialog.geometry("600x400")
        
        text_widget = scrolledtext.ScrolledText(dialog, width=70, height=20)
        text_widget.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        text_widget.insert(tk.END, comparison.to_string())
        text_widget.config(state=tk.DISABLED)
    
    def update_docking_display(self):
        """Update docking display."""
        self.docking_display_text.delete(1.0, tk.END)
        
        if len(self.docking_db) == 0:
            self.docking_display_text.insert(tk.END, "No docking results. Run docking first.")
            return
        
        df = self.docking_db.to_dataframe()
        self.docking_display_text.insert(tk.END, df.to_string())
        
        # Stats
        self.docking_stats_text.delete(1.0, tk.END)
        stats = f"""DOCKING STATISTICS:

Total calculations: {len(self.docking_db)}

Mean affinity: {df['binding_affinity'].mean():.2f} kcal/mol

Best affinity: {df['binding_affinity'].min():.2f} kcal/mol

Good binders (<-7): {df['is_good_binding'].sum()}

Ligands tested: {df['ligand_id'].nunique()}
"""
        self.docking_stats_text.insert(tk.END, stats)
    
    def score_all_candidates(self):
        """Score all candidates with multi-criteria system."""
        if len(self.sequence_db) == 0:
            messagebox.showwarning("Warning", "No candidates to score")
            return
        
        self.status_var.set("Scoring candidates...")
        self.root.update()
        
        # Clear existing scores
        self.ranker = CandidateRanker()
        
        # Score each protein
        structural_scorer = StructuralScorer()
        
        for seq in self.sequence_db.sequences:
            # Get data for this protein
            pocket = self.pocket_db.get_pockets_by_protein(seq.sequence_id)
            structure = self.structure_db.get_structure(seq.sequence_id)
            docking = self.docking_db.get_results_by_protein(seq.sequence_id)
            
            # Calculate scores
            pocket_score = pocket[0].insp3_score if pocket and pocket[0].insp3_score else 5.0
            docking_score = structural_scorer.score_docking(
                min([r.binding_affinity for r in docking])
            ) if docking else 3.0
            quality_score = structural_scorer.score_structure_quality(
                structure.mean_plddt, structure.tm_region_quality
            ) if structure else 3.0
            
            # Create candidate score
            candidate = CandidateScore(
                protein_id=seq.sequence_id,
                pocket_score=pocket_score,
                docking_score=docking_score,
                structure_quality=quality_score,
                phylo_pattern_score=6.0,  # Would need phylogenetic data
                expansion_score=4.0,
                coexpression_score=5.0,  # Would need expression data
                tissue_specificity=6.0,
                localization_score=8.0,  # From prediction
                signal_peptide_score=2.0,
                tm_domain_score=5.0 if seq.n_transmembrane and 4 <= seq.n_transmembrane <= 8 else 2.0,
                pore_architecture=2.0,
                topology_match=2.0
            )
            
            self.ranker.add_candidate(candidate)
        
        self.status_var.set(f"Scored {len(self.ranker.candidates)} candidates")
        self.update_scoring_display()
        
        messagebox.showinfo("Complete", f"Scored {len(self.ranker.candidates)} candidates")
    
    def rank_candidates(self):
        """Rank candidates by score."""
        if len(self.ranker.candidates) == 0:
            messagebox.showwarning("Warning", "Score candidates first")
            return
        
        ranked = self.ranker.rank_candidates()
        self.update_scoring_display()
        messagebox.showinfo("Complete", f"Ranked {len(ranked)} candidates")
    
    def filter_by_priority(self):
        """Filter candidates by priority."""
        if len(self.ranker.candidates) == 0:
            messagebox.showwarning("Warning", "Score candidates first")
            return
        
        priority = self.priority_var.get()
        
        if priority == "HIGH":
            filtered = self.ranker.get_high_priority(threshold=60)
        elif priority == "MEDIUM":
            filtered = [c for c in self.ranker.candidates if c.total_score >= 40]
        else:
            filtered = self.ranker.candidates
        
        messagebox.showinfo("Filtering", f"Found {len(filtered)} {priority.lower()} priority candidates")
    
    def export_rankings(self):
        """Export rankings to CSV."""
        if len(self.ranker.candidates) == 0:
            messagebox.showwarning("Warning", "Score candidates first")
            return
        
        filename = filedialog.asksaveasfilename(
            title="Export Rankings",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if filename:
            try:
                self.ranker.save_rankings(filename)
                messagebox.showinfo("Success", f"Saved rankings to {Path(filename).name}")
            except Exception as e:
                messagebox.showerror("Error", str(e))
    
    def update_scoring_display(self):
        """Update scoring display."""
        self.scoring_display_text.delete(1.0, tk.END)
        
        if len(self.ranker.candidates) == 0:
            self.scoring_display_text.insert(tk.END, "No candidates scored. Click 'Score All Candidates'.")
            return
        
        df = self.ranker.to_dataframe()
        df = df.sort_values('total_score', ascending=False)
        self.scoring_display_text.insert(tk.END, df.to_string())
        
        # Stats
        self.scoring_stats_text.delete(1.0, tk.END)
        stats = f"""SCORING STATISTICS:

Total candidates: {len(self.ranker.candidates)}

HIGH priority (≥60): {len(self.ranker.get_high_priority(60))}

MEDIUM (40-59): {len([c for c in self.ranker.candidates if 40 <= c.total_score < 60])}

LOW (<40): {len([c for c in self.ranker.candidates if c.total_score < 40])}

Top score: {df['total_score'].max():.1f}
Mean score: {df['total_score'].mean():.1f}
"""
        self.scoring_stats_text.insert(tk.END, stats)
    
    def plot_score_distribution(self):
        """Plot score distribution."""
        if len(self.ranker.candidates) == 0:
            messagebox.showwarning("Warning", "Score candidates first")
            return
        
        df = self.ranker.to_dataframe()
        
        self.vis_figure.clear()
        CandidatePlotter.plot_score_distribution(df, output_path=None)
        self.vis_canvas.draw()
        
        # Also save
        CandidatePlotter.plot_score_distribution(
            df, 
            output_path=str(self.output_dir / "score_distribution.png")
        )
    
    def plot_score_heatmap(self):
        """Plot score heatmap."""
        if len(self.ranker.candidates) == 0:
            messagebox.showwarning("Warning", "Score candidates first")
            return
        
        df = self.ranker.to_dataframe()
        
        CandidatePlotter.plot_score_heatmap(
            df,
            output_path=str(self.output_dir / "score_heatmap.png")
        )
        
        messagebox.showinfo("Saved", f"Heatmap saved to {self.output_dir / 'score_heatmap.png'}")
    
    def plot_docking_comparison(self):
        """Plot docking comparison."""
        if len(self.docking_db) == 0:
            messagebox.showwarning("Warning", "Run docking first")
            return
        
        df = self.docking_db.to_dataframe()
        
        CandidatePlotter.plot_docking_comparison(
            df,
            output_path=str(self.output_dir / "docking_comparison.png")
        )
        
        messagebox.showinfo("Saved", f"Comparison saved to {self.output_dir / 'docking_comparison.png'}")
    
    def plot_structure_quality(self):
        """Plot structure quality."""
        if len(self.structure_db) == 0:
            messagebox.showwarning("Warning", "No structures available")
            return
        
        df = self.structure_db.to_dataframe()
        
        CandidatePlotter.plot_structure_quality(
            df,
            output_path=str(self.output_dir / "structure_quality.png")
        )
        
        messagebox.showinfo("Saved", f"Plot saved to {self.output_dir / 'structure_quality.png'}")
    
    def plot_pocket_characteristics(self):
        """Plot pocket characteristics."""
        if len(self.pocket_db) == 0:
            messagebox.showwarning("Warning", "No pockets available")
            return
        
        df = self.pocket_db.to_dataframe()
        
        CandidatePlotter.plot_pocket_characteristics(
            df,
            output_path=str(self.output_dir / 'pocket_characteristics.png')
        )
        
        messagebox.showinfo("Saved", f"Plot saved to {self.output_dir / 'pocket_characteristics.png'}")
    
    def save_all_plots(self):
        """Save all available plots."""
        saved = []
        
        if len(self.ranker.candidates) > 0:
            df = self.ranker.to_dataframe()
            CandidatePlotter.plot_score_distribution(df, str(self.output_dir / "score_distribution.png"))
            saved.append("score_distribution.png")
        
        if len(self.docking_db) > 0:
            df = self.docking_db.to_dataframe()
            CandidatePlotter.plot_docking_comparison(df, str(self.output_dir / "docking_comparison.png"))
            saved.append("docking_comparison.png")
        
        messagebox.showinfo("Saved", f"Saved {len(saved)} plots to {self.output_dir}/")
    
    def generate_report(self):
        """Generate comprehensive report."""
        messagebox.showinfo("Info", "Report generation feature coming soon!")
    
    def show_about(self):
        """Show about dialog."""
        about_text = """Receptor Finder v1.0

Computational discovery of plant InsP₃ and cADPR receptors
through convergent evolution analysis.

Features:
- Sequence filtering and analysis
- Structure prediction interface (AlphaFold2)
- Binding pocket detection and scoring
- Molecular docking simulations
- Multi-criteria scoring (100-point system)
- Comprehensive visualization

Author: George (2025)
License: MIT"""
        
        messagebox.showinfo("About", about_text)
    
    def show_workflow(self):
        """Show workflow guide."""
        workflow_text = """WORKFLOW GUIDE

1. SEQUENCES
   - Load FASTA file or example data
   - Apply filters (length, TM domains)
   - Analyze motifs

2. STRUCTURES
   - Generate AlphaFold2 commands
   - Or simulate predictions (demo)
   - Filter by quality (pLDDT)

3. POCKETS
   - Detect binding pockets
   - Score for InsP₃/cADPR binding
   - Filter high-scoring pockets

4. DOCKING
   - Select ligands
   - Run docking simulations
   - Analyze specificity

5. SCORING
   - Score all candidates (100 pts)
   - Rank by total score
   - Filter by priority class

6. VISUALIZATIONS
   - Generate plots
   - Export results
   - Create report

Expected Results:
- High priority: 20-50 candidates
- Ready for experimental validation"""
        
        dialog = tk.Toplevel(self.root)
        dialog.title("Workflow Guide")
        dialog.geometry("500x600")
        
        text_widget = scrolledtext.ScrolledText(dialog, width=60, height=35)
        text_widget.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        text_widget.insert(tk.END, workflow_text)
        text_widget.config(state=tk.DISABLED)
