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
from gui_tabs import GuiTabsMixin
from gui_actions_analysis import GuiAnalysisMixin
from gui_actions_results import GuiResultsMixin


class ReceptorFinderGUI(GuiTabsMixin, GuiAnalysisMixin, GuiResultsMixin):
    """Main GUI application for Receptor Finder (tabs and actions live in the gui_* mixin modules)."""
    
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
    

def main():
    """Main entry point."""
    root = tk.Tk()
    app = ReceptorFinderGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
