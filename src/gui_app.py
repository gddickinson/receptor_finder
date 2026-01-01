"""
Main GUI Application
====================

Tabbed GUI interface for Receptor Finder - discovering plant InsP₃/cADPR receptors.
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


class ReceptorFinderGUI:
    """Main GUI application for Receptor Finder."""
    
    def __init__(self, root):
        self.root = root
        self.root.title("Receptor Finder v1.0 - Plant InsP₃/cADPR Receptor Discovery")
        self.root.geometry("1400x900")
        
        # Initialize data structures
        self.sequence_db = SequenceDatabase()
        self.structure_db = StructureDatabase()
        self.pocket_db = PocketDatabase()
        self.docking_db = DockingDatabase()
        self.ranker = CandidateRanker()
        
        # Output directory
        self.output_dir = Path("outputs")
        self.output_dir.mkdir(exist_ok=True)
        
        # Create GUI
        self.create_menu()
        self.create_notebook()
        
        # Status bar
        self.status_var = tk.StringVar(value="Ready - Load sequences to begin")
        status_bar = ttk.Label(root, textvariable=self.status_var, relief=tk.SUNKEN)
        status_bar.pack(side=tk.BOTTOM, fill=tk.X)
    
    def create_menu(self):
        """Create menu bar."""
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)
        
        # File menu
        file_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(label="Load Sequences (FASTA)", command=self.load_sequences)
        file_menu.add_command(label="Save Filtered Sequences", command=self.save_sequences)
        file_menu.add_separator()
        file_menu.add_command(label="Load Example Data", command=self.load_example_data)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.root.quit)
        
        # Help menu
        help_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Help", menu=help_menu)
        help_menu.add_command(label="About", command=self.show_about)
        help_menu.add_command(label="Workflow Guide", command=self.show_workflow)
    
    def create_notebook(self):
        """Create tabbed interface."""
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Create tabs
        self.create_sequence_tab()
        self.create_structure_tab()
        self.create_pocket_tab()
        self.create_docking_tab()
        self.create_scoring_tab()
        self.create_visualization_tab()
    
    def create_sequence_tab(self):
        """Create sequence management tab."""
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="1. Sequences")
        
        # Left panel: Controls
        left_panel = ttk.Frame(tab, width=300)
        left_panel.pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=5)
        left_panel.pack_propagate(False)
        
        ttk.Label(left_panel, text="Sequence Analysis", 
                 font=('Arial', 12, 'bold')).pack(pady=5)
        
        ttk.Button(left_panel, text="Load FASTA File",
                  command=self.load_sequences).pack(pady=5, fill=tk.X)
        
        ttk.Button(left_panel, text="Load Example Sequences",
                  command=self.load_example_data).pack(pady=5, fill=tk.X)
        
        ttk.Separator(left_panel, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)
        
        # Filtering options
        ttk.Label(left_panel, text="Filter Criteria", 
                 font=('Arial', 10, 'bold')).pack(pady=5)
        
        # Length filter
        filter_frame = ttk.Frame(left_panel)
        filter_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(filter_frame, text="Length (aa):").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.min_length_var = tk.StringVar(value="200")
        ttk.Entry(filter_frame, textvariable=self.min_length_var, width=10).grid(row=0, column=1, pady=2)
        ttk.Label(filter_frame, text="to").grid(row=0, column=2, pady=2)
        self.max_length_var = tk.StringVar(value="2000")
        ttk.Entry(filter_frame, textvariable=self.max_length_var, width=10).grid(row=0, column=3, pady=2)
        
        # TM domain filter
        ttk.Label(filter_frame, text="TM domains:").grid(row=1, column=0, sticky=tk.W, pady=2)
        self.min_tm_var = tk.StringVar(value="2")
        ttk.Entry(filter_frame, textvariable=self.min_tm_var, width=10).grid(row=1, column=1, pady=2)
        ttk.Label(filter_frame, text="to").grid(row=1, column=2, pady=2)
        self.max_tm_var = tk.StringVar(value="10")
        ttk.Entry(filter_frame, textvariable=self.max_tm_var, width=10).grid(row=1, column=3, pady=2)
        
        ttk.Button(left_panel, text="Apply Filters",
                  command=self.apply_sequence_filters).pack(pady=10, fill=tk.X)
        
        ttk.Button(left_panel, text="Analyze Motifs",
                  command=self.analyze_motifs).pack(pady=5, fill=tk.X)
        
        ttk.Separator(left_panel, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)
        
        # Statistics
        ttk.Label(left_panel, text="Database Statistics:", 
                 font=('Arial', 10, 'bold')).pack(pady=5)
        
        self.seq_stats_text = scrolledtext.ScrolledText(left_panel, width=35, height=10)
        self.seq_stats_text.pack(fill=tk.BOTH, expand=True)
        
        # Right panel: Display
        right_panel = ttk.Frame(tab)
        right_panel.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        ttk.Label(right_panel, text="Sequence Database", 
                 font=('Arial', 12, 'bold')).pack(pady=5)
        
        self.seq_display_text = scrolledtext.ScrolledText(right_panel, width=80, height=35)
        self.seq_display_text.pack(fill=tk.BOTH, expand=True)
    
    def create_structure_tab(self):
        """Create structure prediction tab."""
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="2. Structures")
        
        # Left panel
        left_panel = ttk.Frame(tab, width=300)
        left_panel.pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=5)
        left_panel.pack_propagate(False)
        
        ttk.Label(left_panel, text="Structure Prediction", 
                 font=('Arial', 12, 'bold')).pack(pady=5)
        
        info_text = ("AlphaFold2 structure prediction interface.\n\n"
                    "This module provides commands for running\n"
                    "ColabFold or AlphaFold2.\n\n"
                    "For actual predictions, use external tools.")
        
        ttk.Label(left_panel, text=info_text, justify=tk.LEFT,
                 wraplength=280).pack(pady=10)
        
        ttk.Button(left_panel, text="Generate ColabFold Command",
                  command=self.generate_colabfold_command).pack(pady=5, fill=tk.X)
        
        ttk.Button(left_panel, text="Simulate Predictions (Demo)",
                  command=self.simulate_structure_prediction).pack(pady=5, fill=tk.X)
        
        ttk.Button(left_panel, text="Load Predicted Structures",
                  command=self.load_structures).pack(pady=5, fill=tk.X)
        
        ttk.Separator(left_panel, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)
        
        ttk.Label(left_panel, text="Quality Filters:", 
                 font=('Arial', 10, 'bold')).pack(pady=5)
        
        filter_frame = ttk.Frame(left_panel)
        filter_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(filter_frame, text="Min pLDDT:").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.min_plddt_var = tk.StringVar(value="70")
        ttk.Entry(filter_frame, textvariable=self.min_plddt_var, width=10).grid(row=0, column=1, pady=2)
        
        ttk.Button(left_panel, text="Filter by Quality",
                  command=self.filter_structures).pack(pady=5, fill=tk.X)
        
        ttk.Separator(left_panel, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)
        
        self.struct_stats_text = scrolledtext.ScrolledText(left_panel, width=35, height=12)
        self.struct_stats_text.pack(fill=tk.BOTH, expand=True)
        
        # Right panel
        right_panel = ttk.Frame(tab)
        right_panel.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        ttk.Label(right_panel, text="Predicted Structures", 
                 font=('Arial', 12, 'bold')).pack(pady=5)
        
        self.struct_display_text = scrolledtext.ScrolledText(right_panel, width=80, height=35)
        self.struct_display_text.pack(fill=tk.BOTH, expand=True)
    
    def create_pocket_tab(self):
        """Create pocket analysis tab."""
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="3. Pockets")
        
        # Left panel
        left_panel = ttk.Frame(tab, width=300)
        left_panel.pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=5)
        left_panel.pack_propagate(False)
        
        ttk.Label(left_panel, text="Binding Pocket Analysis", 
                 font=('Arial', 12, 'bold')).pack(pady=5)
        
        ttk.Button(left_panel, text="Detect Pockets (Demo)",
                  command=self.detect_pockets).pack(pady=5, fill=tk.X)
        
        ttk.Button(left_panel, text="Score for InsP₃ Binding",
                  command=self.score_insp3_binding).pack(pady=5, fill=tk.X)
        
        ttk.Button(left_panel, text="Score for cADPR Binding",
                  command=self.score_cadpr_binding).pack(pady=5, fill=tk.X)
        
        ttk.Separator(left_panel, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)
        
        ttk.Label(left_panel, text="Filter Options:", 
                 font=('Arial', 10, 'bold')).pack(pady=5)
        
        filter_frame = ttk.Frame(left_panel)
        filter_frame.pack(fill=tk.X, pady=5)
        
        ttk.Label(filter_frame, text="Min Score:").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.min_pocket_score_var = tk.StringVar(value="7.0")
        ttk.Entry(filter_frame, textvariable=self.min_pocket_score_var, width=10).grid(row=0, column=1, pady=2)
        
        ttk.Button(left_panel, text="Filter High-Scoring Pockets",
                  command=self.filter_pockets).pack(pady=5, fill=tk.X)
        
        ttk.Separator(left_panel, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)
        
        self.pocket_stats_text = scrolledtext.ScrolledText(left_panel, width=35, height=15)
        self.pocket_stats_text.pack(fill=tk.BOTH, expand=True)
        
        # Right panel
        right_panel = ttk.Frame(tab)
        right_panel.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        ttk.Label(right_panel, text="Detected Binding Pockets", 
                 font=('Arial', 12, 'bold')).pack(pady=5)
        
        self.pocket_display_text = scrolledtext.ScrolledText(right_panel, width=80, height=35)
        self.pocket_display_text.pack(fill=tk.BOTH, expand=True)
    
    def create_docking_tab(self):
        """Create molecular docking tab."""
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="4. Docking")
        
        # Left panel
        left_panel = ttk.Frame(tab, width=300)
        left_panel.pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=5)
        left_panel.pack_propagate(False)
        
        ttk.Label(left_panel, text="Molecular Docking", 
                 font=('Arial', 12, 'bold')).pack(pady=5)
        
        ttk.Label(left_panel, text="Select Ligands:", 
                 font=('Arial', 10, 'bold')).pack(pady=5)
        
        # Ligand checkboxes
        self.ligand_vars = {}
        ligands = ['InsP3', 'InsP4', 'cADPR', 'NAADP']
        for ligand in ligands:
            var = tk.BooleanVar(value=(ligand in ['InsP3', 'cADPR']))
            self.ligand_vars[ligand] = var
            ttk.Checkbutton(left_panel, text=ligand, variable=var).pack(anchor=tk.W)
        
        ttk.Separator(left_panel, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)
        
        ttk.Button(left_panel, text="Run Docking (Demo)",
                  command=self.run_docking).pack(pady=5, fill=tk.X)
        
        ttk.Button(left_panel, text="Analyze Specificity",
                  command=self.analyze_specificity).pack(pady=5, fill=tk.X)
        
        ttk.Button(left_panel, text="Compare Ligands",
                  command=self.compare_ligands).pack(pady=5, fill=tk.X)
        
        ttk.Separator(left_panel, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)
        
        self.docking_stats_text = scrolledtext.ScrolledText(left_panel, width=35, height=15)
        self.docking_stats_text.pack(fill=tk.BOTH, expand=True)
        
        # Right panel
        right_panel = ttk.Frame(tab)
        right_panel.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        ttk.Label(right_panel, text="Docking Results", 
                 font=('Arial', 12, 'bold')).pack(pady=5)
        
        self.docking_display_text = scrolledtext.ScrolledText(right_panel, width=80, height=35)
        self.docking_display_text.pack(fill=tk.BOTH, expand=True)
    
    def create_scoring_tab(self):
        """Create scoring and ranking tab."""
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="5. Scoring & Ranking")
        
        # Left panel
        left_panel = ttk.Frame(tab, width=300)
        left_panel.pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=5)
        left_panel.pack_propagate(False)
        
        ttk.Label(left_panel, text="Multi-Criteria Scoring", 
                 font=('Arial', 12, 'bold')).pack(pady=5)
        
        ttk.Button(left_panel, text="Score All Candidates",
                  command=self.score_all_candidates).pack(pady=5, fill=tk.X)
        
        ttk.Button(left_panel, text="Rank by Total Score",
                  command=self.rank_candidates).pack(pady=5, fill=tk.X)
        
        ttk.Separator(left_panel, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)
        
        ttk.Label(left_panel, text="Priority Filter:", 
                 font=('Arial', 10, 'bold')).pack(pady=5)
        
        self.priority_var = tk.StringVar(value="ALL")
        ttk.Radiobutton(left_panel, text="All Candidates", 
                       variable=self.priority_var, value="ALL").pack(anchor=tk.W)
        ttk.Radiobutton(left_panel, text="High Priority (≥60)", 
                       variable=self.priority_var, value="HIGH").pack(anchor=tk.W)
        ttk.Radiobutton(left_panel, text="Medium+ (≥40)", 
                       variable=self.priority_var, value="MEDIUM").pack(anchor=tk.W)
        
        ttk.Button(left_panel, text="Apply Filter",
                  command=self.filter_by_priority).pack(pady=5, fill=tk.X)
        
        ttk.Separator(left_panel, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)
        
        ttk.Button(left_panel, text="Export Rankings (CSV)",
                  command=self.export_rankings).pack(pady=5, fill=tk.X)
        
        ttk.Separator(left_panel, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)
        
        self.scoring_stats_text = scrolledtext.ScrolledText(left_panel, width=35, height=12)
        self.scoring_stats_text.pack(fill=tk.BOTH, expand=True)
        
        # Right panel
        right_panel = ttk.Frame(tab)
        right_panel.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        ttk.Label(right_panel, text="Candidate Rankings", 
                 font=('Arial', 12, 'bold')).pack(pady=5)
        
        self.scoring_display_text = scrolledtext.ScrolledText(right_panel, width=80, height=35)
        self.scoring_display_text.pack(fill=tk.BOTH, expand=True)
    
    def create_visualization_tab(self):
        """Create visualization tab."""
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text="6. Visualizations")
        
        # Left panel
        left_panel = ttk.Frame(tab, width=300)
        left_panel.pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=5)
        left_panel.pack_propagate(False)
        
        ttk.Label(left_panel, text="Generate Plots", 
                 font=('Arial', 12, 'bold')).pack(pady=5)
        
        ttk.Button(left_panel, text="Score Distribution",
                  command=self.plot_score_distribution).pack(pady=5, fill=tk.X)
        
        ttk.Button(left_panel, text="Score Heatmap",
                  command=self.plot_score_heatmap).pack(pady=5, fill=tk.X)
        
        ttk.Button(left_panel, text="Docking Comparison",
                  command=self.plot_docking_comparison).pack(pady=5, fill=tk.X)
        
        ttk.Button(left_panel, text="Structure Quality",
                  command=self.plot_structure_quality).pack(pady=5, fill=tk.X)
        
        ttk.Button(left_panel, text="Pocket Characteristics",
                  command=self.plot_pocket_characteristics).pack(pady=5, fill=tk.X)
        
        ttk.Separator(left_panel, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)
        
        ttk.Label(left_panel, text="Export:", 
                 font=('Arial', 10, 'bold')).pack(pady=5)
        
        ttk.Button(left_panel, text="Save All Plots",
                  command=self.save_all_plots).pack(pady=5, fill=tk.X)
        
        ttk.Button(left_panel, text="Generate Report",
                  command=self.generate_report).pack(pady=5, fill=tk.X)
        
        # Right panel
        right_panel = ttk.Frame(tab)
        right_panel.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        ttk.Label(right_panel, text="Visualization", 
                 font=('Arial', 12, 'bold')).pack(pady=5)
        
        self.vis_figure = Figure(figsize=(10, 8))
        self.vis_canvas = FigureCanvasTkAgg(self.vis_figure, right_panel)
        self.vis_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
    
    # Callback methods
    
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


def main():
    """Main entry point."""
    root = tk.Tk()
    app = ReceptorFinderGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
