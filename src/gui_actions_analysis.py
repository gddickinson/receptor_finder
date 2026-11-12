"""
Sequence, structure and pocket actions of the Receptor Finder GUI.

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


class GuiAnalysisMixin:
    """Sequence, structure and pocket actions of the Receptor Finder GUI."""

    def load_sequences(self):
        """Load sequences from FASTA file."""
        filename = filedialog.askopenfilename(
            title="Load FASTA File",
            filetypes=[("FASTA files", "*.fasta *.fa *.faa"), ("All files", "*.*")]
        )
        if filename:
            try:
                species = "Unknown"  # Could add dialog to get species name
                self.sequence_db.load_from_fasta(filename, species=species)
                self.status_var.set(f"Loaded {len(self.sequence_db)} sequences from {Path(filename).name}")
                self.update_sequence_display()
                messagebox.showinfo("Success", f"Loaded {len(self.sequence_db)} sequences!")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load sequences:\n{str(e)}")
    
    def load_example_data(self):
        """Load example sequences."""
        from sequence_analysis import create_example_database
        self.sequence_db = create_example_database()
        self.status_var.set(f"Loaded {len(self.sequence_db)} example sequences")
        self.update_sequence_display()
        messagebox.showinfo("Success", "Example data loaded!")
    
    def save_sequences(self):
        """Save filtered sequences to FASTA."""
        if len(self.sequence_db) == 0:
            messagebox.showwarning("Warning", "No sequences to save")
            return
        
        filename = filedialog.asksaveasfilename(
            title="Save Sequences",
            defaultextension=".fasta",
            filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")]
        )
        if filename:
            try:
                self.sequence_db.save_to_fasta(filename)
                self.status_var.set(f"Saved {len(self.sequence_db)} sequences")
                messagebox.showinfo("Success", f"Saved sequences to {Path(filename).name}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save:\n{str(e)}")
    
    def apply_sequence_filters(self):
        """Apply filtering criteria to sequences."""
        if len(self.sequence_db) == 0:
            messagebox.showwarning("Warning", "Load sequences first")
            return
        
        try:
            min_len = int(self.min_length_var.get())
            max_len = int(self.max_length_var.get())
            min_tm = int(self.min_tm_var.get())
            max_tm = int(self.max_tm_var.get())
            
            original_count = len(self.sequence_db)
            self.sequence_db = self.sequence_db.filter_database(
                min_length=min_len,
                max_length=max_len,
                min_tm=min_tm,
                max_tm=max_tm
            )
            
            self.status_var.set(f"Filtered: {len(self.sequence_db)}/{original_count} sequences retained")
            self.update_sequence_display()
            
            messagebox.showinfo("Filtering Complete", 
                              f"Retained {len(self.sequence_db)} of {original_count} sequences")
        except Exception as e:
            messagebox.showerror("Error", f"Filtering failed:\n{str(e)}")
    
    def analyze_motifs(self):
        """Analyze motifs in sequences."""
        if len(self.sequence_db) == 0:
            messagebox.showwarning("Warning", "Load sequences first")
            return
        
        # Analyze first 10 sequences
        results = "MOTIF ANALYSIS (first 10 sequences):\n\n"
        
        for seq in self.sequence_db.sequences[:10]:
            results += f"{seq.sequence_id}:\n"
            motif_counts = MotifScanner.get_motif_counts(seq.sequence)
            for motif, count in motif_counts.items():
                if count > 0:
                    results += f"  {motif}: {count}\n"
            results += "\n"
        
        # Display in dialog
        dialog = tk.Toplevel(self.root)
        dialog.title("Motif Analysis Results")
        dialog.geometry("500x600")
        
        text_widget = scrolledtext.ScrolledText(dialog, width=60, height=35)
        text_widget.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        text_widget.insert(tk.END, results)
        text_widget.config(state=tk.DISABLED)
    
    def update_sequence_display(self):
        """Update sequence display."""
        self.seq_display_text.delete(1.0, tk.END)
        
        if len(self.sequence_db) == 0:
            self.seq_display_text.insert(tk.END, "No sequences loaded. Use File menu to load FASTA file.")
            return
        
        # Show table
        df = self.sequence_db.to_dataframe()
        self.seq_display_text.insert(tk.END, df.to_string())
        
        # Update stats
        self.seq_stats_text.delete(1.0, tk.END)
        stats = f"""DATABASE STATISTICS:

Total sequences: {len(self.sequence_db)}

Length range: {df['length'].min()}-{df['length'].max()} aa

TM domains:
  Mean: {df['n_transmembrane'].mean():.1f}
  Range: {df['n_transmembrane'].min()}-{df['n_transmembrane'].max()}

Species: {df['species'].nunique()}
"""
        self.seq_stats_text.insert(tk.END, stats)
    
    def generate_colabfold_command(self):
        """Generate ColabFold command."""
        if len(self.sequence_db) == 0:
            messagebox.showwarning("Warning", "Load sequences first")
            return
        
        # Save sequences to temp file
        temp_fasta = self.output_dir / "sequences_for_prediction.fasta"
        self.sequence_db.save_to_fasta(str(temp_fasta))
        
        cmd = AlphaFold2Interface.get_colabfold_command(
            str(temp_fasta),
            str(self.output_dir / "structures")
        )
        
        # Show command in dialog
        dialog = tk.Toplevel(self.root)
        dialog.title("ColabFold Command")
        dialog.geometry("600x200")
        
        ttk.Label(dialog, text="Run this command to predict structures:",
                 font=('Arial', 10, 'bold')).pack(pady=10)
        
        text_widget = scrolledtext.ScrolledText(dialog, width=70, height=6)
        text_widget.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        text_widget.insert(tk.END, cmd)
        text_widget.config(state=tk.DISABLED)
        
        ttk.Button(dialog, text="Copy to Clipboard",
                  command=lambda: self.root.clipboard_append(cmd)).pack(pady=5)
    
    def simulate_structure_prediction(self):
        """Simulate structure prediction (demo mode)."""
        if len(self.sequence_db) == 0:
            messagebox.showwarning("Warning", "Load sequences first")
            return
        
        self.status_var.set("Simulating structure predictions...")
        self.root.update()
        
        interface = AlphaFold2Interface(output_dir=str(self.output_dir / "structures"))
        
        for seq in self.sequence_db.sequences[:min(10, len(self.sequence_db))]:
            structure = interface.predict_structure(seq.sequence, seq.sequence_id)
            self.structure_db.add_structure(structure)
        
        self.status_var.set(f"Simulated predictions for {len(self.structure_db)} proteins")
        self.update_structure_display()
        
        messagebox.showinfo("Complete", f"Predicted {len(self.structure_db)} structures (demo mode)")
    
    def load_structures(self):
        """Load predicted structures."""
        messagebox.showinfo("Info", "In production, would load PDB files with pLDDT scores")
    
    def filter_structures(self):
        """Filter structures by quality."""
        if len(self.structure_db) == 0:
            messagebox.showwarning("Warning", "No structures to filter")
            return
        
        try:
            min_plddt = float(self.min_plddt_var.get())
            original_count = len(self.structure_db)
            
            high_quality = self.structure_db.get_high_quality_structures(min_plddt=min_plddt)
            
            messagebox.showinfo("Filtering", 
                              f"Found {len(high_quality)} of {original_count} structures with pLDDT ≥ {min_plddt}")
        except Exception as e:
            messagebox.showerror("Error", str(e))
    
    def update_structure_display(self):
        """Update structure display."""
        self.struct_display_text.delete(1.0, tk.END)
        
        if len(self.structure_db) == 0:
            self.struct_display_text.insert(tk.END, "No structures loaded. Run predictions or load PDB files.")
            return
        
        df = self.structure_db.to_dataframe()
        self.struct_display_text.insert(tk.END, df.to_string())
        
        # Stats
        self.struct_stats_text.delete(1.0, tk.END)
        stats = f"""STRUCTURE STATISTICS:

Total structures: {len(self.structure_db)}

Mean pLDDT: {df['mean_plddt'].mean():.1f}

High quality (>70): {df['is_high_quality'].sum()}

Quality range: {df['mean_plddt'].min():.1f}-{df['mean_plddt'].max():.1f}
"""
        self.struct_stats_text.insert(tk.END, stats)
    
    def detect_pockets(self):
        """Detect binding pockets."""
        if len(self.structure_db) == 0:
            messagebox.showwarning("Warning", "No structures available")
            return
        
        self.status_var.set("Detecting binding pockets...")
        self.root.update()
        
        detector = PocketDetector()
        
        for structure in list(self.structure_db.structures.values())[:min(10, len(self.structure_db))]:
            pockets = detector.detect_pockets(structure.pdb_file, structure.protein_id)
            for pocket in pockets:
                self.pocket_db.add_pocket(pocket)
        
        self.status_var.set(f"Detected {len(self.pocket_db)} pockets")
        self.update_pocket_display()
        
        messagebox.showinfo("Complete", f"Detected {len(self.pocket_db)} pockets (demo mode)")
    
    def score_insp3_binding(self):
        """Score pockets for InsP3 binding."""
        if len(self.pocket_db) == 0:
            messagebox.showwarning("Warning", "Detect pockets first")
            return
        
        scorer = PocketScorer()
        for pocket in self.pocket_db.pockets.values():
            scorer.score_insp3_binding(pocket)
        
        self.update_pocket_display()
        messagebox.showinfo("Complete", "Scored all pockets for InsP₃ binding")
    
    def score_cadpr_binding(self):
        """Score pockets for cADPR binding."""
        if len(self.pocket_db) == 0:
            messagebox.showwarning("Warning", "Detect pockets first")
            return
        
        scorer = PocketScorer()
        for pocket in self.pocket_db.pockets.values():
            scorer.score_cadpr_binding(pocket)
        
        self.update_pocket_display()
        messagebox.showinfo("Complete", "Scored all pockets for cADPR binding")
    
    def filter_pockets(self):
        """Filter high-scoring pockets."""
        if len(self.pocket_db) == 0:
            messagebox.showwarning("Warning", "No pockets to filter")
            return
        
        try:
            min_score = float(self.min_pocket_score_var.get())
            high_scoring = self.pocket_db.get_high_scoring_pockets('InsP3', min_score=min_score)
            
            messagebox.showinfo("Filtering",
                              f"Found {len(high_scoring)} pockets with InsP₃ score ≥ {min_score}")
        except Exception as e:
            messagebox.showerror("Error", str(e))
    
    def update_pocket_display(self):
        """Update pocket display."""
        self.pocket_display_text.delete(1.0, tk.END)
        
        if len(self.pocket_db) == 0:
            self.pocket_display_text.insert(tk.END, "No pockets detected. Run pocket detection first.")
            return
        
        df = self.pocket_db.to_dataframe()
        self.pocket_display_text.insert(tk.END, df.to_string())
        
        # Stats
        self.pocket_stats_text.delete(1.0, tk.END)
        stats = f"""POCKET STATISTICS:

Total pockets: {len(self.pocket_db)}

Mean volume: {df['volume'].mean():.1f} ų

Mean InsP₃ score: {df['insp3_score'].mean():.1f}/13

High-scoring InsP₃ (>7): {(df['insp3_score'] > 7).sum()}

Mean cADPR score: {df['cadpr_score'].mean():.1f}/12
"""
        self.pocket_stats_text.insert(tk.END, stats)
