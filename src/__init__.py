"""
Receptor Finder
===============

A comprehensive toolkit for discovering candidate InsP₃ and cADPR receptors 
in plant genomes through convergent evolution analysis.

Modules:
--------
- sequence_analysis: Sequence retrieval, filtering, and homolog detection
- structure_prediction: AlphaFold2 structure prediction interface
- binding_pocket: Binding pocket detection and characterization
- docking: Molecular docking simulations
- phylogenetics: Phylogenetic analysis and evolutionary patterns
- scoring: Multi-criteria candidate ranking
- visualization: Plotting and visualization tools
- utils: Utility functions and helpers

Author: George (2025)
License: MIT
"""

__version__ = "1.0.0"
__author__ = "George"

from . import sequence_analysis
from . import structure_prediction
from . import binding_pocket
from . import docking
from . import phylogenetics
from . import scoring
from . import visualization
from . import utils

__all__ = [
    'sequence_analysis',
    'structure_prediction',
    'binding_pocket',
    'docking',
    'phylogenetics',
    'scoring',
    'visualization',
    'utils'
]
