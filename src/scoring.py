"""
Scoring Module
==============

Multi-criteria scoring system for ranking receptor candidates.
Integrates structural, functional, and evolutionary evidence.
"""

import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple
from pathlib import Path

# The individual scorers live in scorers.py; re-exported here so
# "from scoring import StructuralScorer" keeps working.
try:
    from .scorers import (StructuralScorer, EvolutionaryScorer, ExpressionScorer,  # noqa: F401
                          LocalizationScorer, TopologyScorer)
except ImportError:  # imported as a top-level module with src/ on sys.path
    from scorers import (StructuralScorer, EvolutionaryScorer, ExpressionScorer,  # noqa: F401
                         LocalizationScorer, TopologyScorer)


@dataclass
class CandidateScore:
    """Multi-criteria score for a receptor candidate."""
    
    protein_id: str
    
    # Structural scores (max 30 points)
    pocket_score: float = 0.0        # 0-13 for InsP3 or 0-12 for cADPR
    docking_score: float = 0.0       # 0-10 based on affinity
    structure_quality: float = 0.0    # 0-7 based on pLDDT
    
    # Evolutionary scores (max 20 points)
    phylo_pattern_score: float = 0.0  # 0-10 for distribution pattern
    expansion_score: float = 0.0       # 0-10 for copy number variation
    
    # Expression scores (max 20 points)
    coexpression_score: float = 0.0   # 0-10 correlation with Ca2+ genes
    tissue_specificity: float = 0.0    # 0-10 expected tissue expression
    
    # Localization scores (max 15 points)
    localization_score: float = 0.0   # 0-10 for tonoplast/ER
    signal_peptide_score: float = 0.0  # 0-5 for targeting signals
    
    # Topology scores (max 15 points)
    tm_domain_score: float = 0.0      # 0-5 for appropriate TM count
    pore_architecture: float = 0.0     # 0-5 for channel-like structure
    topology_match: float = 0.0        # 0-5 similarity to known channels
    
    @property
    def total_score(self) -> float:
        """Calculate total score (max 100)."""
        return (self.pocket_score + self.docking_score + self.structure_quality +
                self.phylo_pattern_score + self.expansion_score +
                self.coexpression_score + self.tissue_specificity +
                self.localization_score + self.signal_peptide_score +
                self.tm_domain_score + self.pore_architecture + self.topology_match)
    
    @property
    def priority_class(self) -> str:
        """Classify candidate priority."""
        total = self.total_score
        if total >= 60:
            return "HIGH"
        elif total >= 40:
            return "MEDIUM"
        else:
            return "LOW"


class CandidateRanker:
    """Rank and prioritize receptor candidates."""
    
    def __init__(self):
        self.candidates: List[CandidateScore] = []
    
    def add_candidate(self, candidate: CandidateScore):
        """Add candidate score."""
        self.candidates.append(candidate)
    
    def rank_candidates(self, by_score: str = "total") -> List[CandidateScore]:
        """
        Rank candidates by score.
        
        Parameters
        ----------
        by_score : str
            Which score to rank by ('total', 'pocket', 'docking', etc.)
            
        Returns
        -------
        ranked : List[CandidateScore]
            Ranked candidates (highest first)
        """
        if by_score == "total":
            return sorted(self.candidates, key=lambda x: x.total_score, reverse=True)
        # Accept a field name ('structure_quality', 'tissue_specificity') or
        # the short form of a *_score field ('pocket' -> 'pocket_score').
        fields = CandidateScore.__dataclass_fields__
        if by_score in fields and by_score != "protein_id":
            attr = by_score
        elif f"{by_score}_score" in fields:
            attr = f"{by_score}_score"
        else:
            raise ValueError(f"Unknown score field: {by_score!r}")
        return sorted(self.candidates, key=lambda x: getattr(x, attr), reverse=True)
    
    def get_high_priority(self, threshold: float = 60.0) -> List[CandidateScore]:
        """Get high-priority candidates."""
        return [c for c in self.candidates if c.total_score >= threshold]
    
    def to_dataframe(self) -> pd.DataFrame:
        """Convert to DataFrame."""
        data = []
        for candidate in self.candidates:
            data.append({
                'protein_id': candidate.protein_id,
                'total_score': candidate.total_score,
                'priority_class': candidate.priority_class,
                'pocket_score': candidate.pocket_score,
                'docking_score': candidate.docking_score,
                'structure_quality': candidate.structure_quality,
                'phylo_pattern_score': candidate.phylo_pattern_score,
                'expansion_score': candidate.expansion_score,
                'coexpression_score': candidate.coexpression_score,
                'localization_score': candidate.localization_score,
                'tm_domain_score': candidate.tm_domain_score
            })
        return pd.DataFrame(data)
    
    def save_rankings(self, filepath: str):
        """Save rankings to CSV."""
        df = self.to_dataframe()
        df = df.sort_values('total_score', ascending=False)
        df.to_csv(filepath, index=False)


if __name__ == "__main__":
    # Test code
    print("Testing Scoring Module\n")
    
    # Test 1: Structural scoring
    print("Test 1: Structural Scoring")
    scorer = StructuralScorer()
    
    pocket_score = scorer.score_pocket(10.5, max_pocket_score=13)
    docking_score = scorer.score_docking(-8.5)
    quality_score = scorer.score_structure_quality(85.0, tm_region_quality=78.0)
    
    print(f"  Pocket score: {pocket_score:.1f}/13")
    print(f"  Docking score: {docking_score:.1f}/10")
    print(f"  Quality score: {quality_score:.1f}/7")
    
    # Test 2: Evolutionary scoring
    print("\nTest 2: Evolutionary Scoring")
    evo_scorer = EvolutionaryScorer()
    
    phylo_score = evo_scorer.score_phylogenetic_pattern("Chlamydomonas_loss")
    expansion_score = evo_scorer.score_gene_expansion(copy_number=5, ancestral_copy_number=1)
    
    print(f"  Phylogenetic pattern score: {phylo_score:.1f}/10")
    print(f"  Expansion score: {expansion_score:.1f}/10")
    
    # Test 3: Complete candidate score
    print("\nTest 3: Complete Candidate Scoring")
    
    candidate = CandidateScore(
        protein_id="AT1G12345",
        pocket_score=10.5,
        docking_score=8.5,
        structure_quality=6.5,
        phylo_pattern_score=10.0,
        expansion_score=5.0,
        coexpression_score=7.5,
        tissue_specificity=8.0,
        localization_score=10.0,
        signal_peptide_score=5.0,
        tm_domain_score=5.0,
        pore_architecture=4.0,
        topology_match=4.0
    )
    
    print(f"  Protein: {candidate.protein_id}")
    print(f"  Total score: {candidate.total_score:.1f}/100")
    print(f"  Priority: {candidate.priority_class}")
    
    # Test 4: Candidate ranking
    print("\nTest 4: Candidate Ranking")
    
    ranker = CandidateRanker()
    
    # Create mock candidates
    for i in range(5):
        mock_candidate = CandidateScore(
            protein_id=f"AT{i}G{12345+i}",
            pocket_score=8.0 + i,
            docking_score=7.0 + i*0.5,
            structure_quality=5.0 + i*0.3,
            phylo_pattern_score=6.0 + i*0.8,
            expansion_score=4.0 + i*0.6,
            localization_score=8.0,
            tm_domain_score=5.0
        )
        ranker.add_candidate(mock_candidate)
    
    # Rank
    ranked = ranker.rank_candidates()
    
    print(f"  Total candidates: {len(ranker.candidates)}")
    print(f"  High priority: {len(ranker.get_high_priority(threshold=50))}")
    
    print("\n  Top 3 candidates:")
    for i, cand in enumerate(ranked[:3]):
        print(f"    {i+1}. {cand.protein_id}: {cand.total_score:.1f} pts ({cand.priority_class})")
    
    print("\nTest complete!")
