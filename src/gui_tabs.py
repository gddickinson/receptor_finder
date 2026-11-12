"""
Menu bar and the six notebook tabs of the Receptor Finder GUI.

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


class GuiTabsMixin:
    """Menu bar and the six notebook tabs of the Receptor Finder GUI."""

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
