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


class StructuralScorer:
    """Score structural features."""
    
    @staticmethod
    def score_pocket(pocket_score: float, max_pocket_score: float = 13.0) -> float:
        """
        Normalize pocket score to 0-13 scale.
        
        Parameters
        ----------
        pocket_score : float
            Raw pocket score from PocketScorer
        max_pocket_score : float
            Maximum possible pocket score
            
        Returns
        -------
        normalized_score : float
            Score 0-13
        """
        return min(pocket_score, max_pocket_score)
    
    @staticmethod
    def score_docking(binding_affinity: float) -> float:
        """
        Score docking affinity (0-10 points).
        
        Affinity scale (kcal/mol):
        - < -9.0: Excellent (10 pts)
        - -7.0 to -9.0: Good (7-10 pts)
        - -5.0 to -7.0: Moderate (3-7 pts)
        - > -5.0: Poor (0-3 pts)
        
        Parameters
        ----------
        binding_affinity : float
            Binding affinity in kcal/mol (more negative = better)
            
        Returns
        -------
        score : float
            Score 0-10
        """
        if binding_affinity < -9.0:
            return 10.0
        elif binding_affinity < -7.0:
            # Linear interpolation: -9 to -7 maps to 10 to 7
            return 10.0 + (binding_affinity + 9.0) * (3.0 / 2.0)
        elif binding_affinity < -5.0:
            # Linear interpolation: -7 to -5 maps to 7 to 3
            return 7.0 + (binding_affinity + 7.0) * (4.0 / 2.0)
        else:
            # > -5.0: scale from 3 to 0
            return max(0, 3.0 + (binding_affinity + 5.0) * 0.6)
    
    @staticmethod
    def score_structure_quality(mean_plddt: float,
                                tm_region_quality: Optional[float] = None) -> float:
        """
        Score structure quality (0-7 points).
        
        Based on AlphaFold2 pLDDT:
        - > 90: Very high confidence (7 pts)
        - 70-90: High confidence (5-7 pts)
        - 50-70: Medium confidence (2-5 pts)
        - < 50: Low confidence (0-2 pts)
        
        TM region quality adds up to 2 extra points.
        
        Returns
        -------
        score : float
            Score 0-7
        """
        # Base score from mean pLDDT
        if mean_plddt > 90:
            score = 7.0
        elif mean_plddt > 70:
            score = 5.0 + (mean_plddt - 70) * (2.0 / 20.0)
        elif mean_plddt > 50:
            score = 2.0 + (mean_plddt - 50) * (3.0 / 20.0)
        else:
            score = (mean_plddt / 50.0) * 2.0
        
        # Penalize if TM region quality is poor
        if tm_region_quality is not None and tm_region_quality < 70:
            score *= (tm_region_quality / 70.0)
        
        return min(score, 7.0)


class EvolutionaryScorer:
    """Score evolutionary patterns."""
    
    @staticmethod
    def score_phylogenetic_pattern(pattern: str) -> float:
        """
        Score phylogenetic distribution pattern (0-10 points).
        
        Patterns:
        - "Chlamydomonas_loss": High score (10 pts) - compensatory evolution
        - "Land_plant_innovation": Medium-high (8 pts)
        - "Angiosperm_specific": Medium (6 pts)
        - "Conserved": Lower (4 pts) - less likely to be novel receptor
        - "Sporadic": Low (2 pts)
        
        Parameters
        ----------
        pattern : str
            Phylogenetic pattern classification
            
        Returns
        -------
        score : float
            Score 0-10
        """
        pattern_scores = {
            "Chlamydomonas_loss": 10.0,
            "Land_plant_innovation": 8.0,
            "Angiosperm_specific": 6.0,
            "Bryophyte_origin": 7.0,
            "Conserved": 4.0,
            "Sporadic": 2.0,
            "Unknown": 0.0
        }
        
        return pattern_scores.get(pattern, 0.0)
    
    @staticmethod
    def score_gene_expansion(copy_number: int,
                            ancestral_copy_number: int = 1) -> float:
        """
        Score gene family expansion (0-10 points).
        
        Expansion suggests functional importance.
        
        Parameters
        ----------
        copy_number : int
            Number of copies in focal species
        ancestral_copy_number : int
            Estimated ancestral copy number
            
        Returns
        -------
        score : float
            Score 0-10
        """
        expansion_factor = copy_number / max(ancestral_copy_number, 1)
        
        if expansion_factor >= 10:
            return 10.0
        elif expansion_factor >= 5:
            return 8.0
        elif expansion_factor >= 2:
            return 5.0 + (expansion_factor - 2) * (3.0 / 3.0)
        else:
            return expansion_factor * 2.5


class ExpressionScorer:
    """Score expression patterns."""
    
    @staticmethod
    def score_coexpression(correlation: float) -> float:
        """
        Score coexpression with Ca²⁺ signaling genes (0-10 points).
        
        Based on Pearson correlation coefficient.
        
        Parameters
        ----------
        correlation : float
            Correlation coefficient (-1 to 1)
            
        Returns
        -------
        score : float
            Score 0-10
        """
        # Convert correlation to score
        if correlation > 0.8:
            return 10.0
        elif correlation > 0.6:
            return 7.0 + (correlation - 0.6) * (3.0 / 0.2)
        elif correlation > 0.4:
            return 4.0 + (correlation - 0.4) * (3.0 / 0.2)
        elif correlation > 0.2:
            return 1.0 + (correlation - 0.2) * (3.0 / 0.2)
        else:
            return max(0, correlation * 5.0)
    
    @staticmethod
    def score_tissue_specificity(tissues: List[str],
                                 expected_tissues: List[str]) -> float:
        """
        Score expression in expected tissues (0-10 points).
        
        Expected tissues for Ca²⁺ receptors might include:
        - Guard cells (for stomatal regulation)
        - Root hairs (for Ca²⁺ waves)
        - Pollen tubes (for tip growth)
        
        Parameters
        ----------
        tissues : List[str]
            Tissues where gene is expressed
        expected_tissues : List[str]
            Expected tissues for receptor
            
        Returns
        -------
        score : float
            Score 0-10
        """
        if not expected_tissues:
            return 5.0  # No expectation
        
        overlap = len(set(tissues) & set(expected_tissues))
        fraction = overlap / len(expected_tissues)
        
        return fraction * 10.0


class LocalizationScorer:
    """Score subcellular localization."""
    
    @staticmethod
    def score_localization(predicted_location: str) -> float:
        """
        Score subcellular localization (0-10 points).
        
        For InsP₃/cADPR receptors, expect:
        - Tonoplast/vacuole: High (10 pts)
        - ER: High (9 pts)
        - Plasma membrane: Medium (6 pts)
        - Other: Low (2 pts)
        
        Parameters
        ----------
        predicted_location : str
            Predicted localization
            
        Returns
        -------
        score : float
            Score 0-10
        """
        localization_scores = {
            "tonoplast": 10.0,
            "vacuole": 10.0,
            "ER": 9.0,
            "endoplasmic_reticulum": 9.0,
            "plasma_membrane": 6.0,
            "endomembrane": 7.0,
            "golgi": 5.0,
            "mitochondria": 2.0,
            "chloroplast": 2.0,
            "nucleus": 1.0,
            "cytosol": 1.0
        }
        
        return localization_scores.get(predicted_location.lower(), 0.0)
    
    @staticmethod
    def score_signal_peptide(has_signal: bool,
                            signal_type: Optional[str] = None) -> float:
        """
        Score signal peptide/targeting (0-5 points).
        
        Parameters
        ----------
        has_signal : bool
            Whether signal peptide detected
        signal_type : Optional[str]
            Type of signal (ER, vacuolar, etc.)
            
        Returns
        -------
        score : float
            Score 0-5
        """
        if not has_signal:
            return 0.0
        
        if signal_type:
            type_scores = {
                "ER": 5.0,
                "vacuolar": 5.0,
                "secretory": 4.0,
                "mitochondrial": 1.0,
                "chloroplast": 1.0
            }
            return type_scores.get(signal_type, 2.0)
        else:
            return 2.0


class TopologyScorer:
    """Score membrane topology."""
    
    @staticmethod
    def score_tm_domains(n_tm: int, expected_tm: Tuple[int, int] = (4, 8)) -> float:
        """
        Score number of TM domains (0-5 points).
        
        Expected: 4-8 TM domains for tetrameric Ca²⁺ channel
        
        Parameters
        ----------
        n_tm : int
            Number of predicted TM domains
        expected_tm : Tuple[int, int]
            Expected range (min, max)
            
        Returns
        -------
        score : float
            Score 0-5
        """
        min_tm, max_tm = expected_tm
        
        if min_tm <= n_tm <= max_tm:
            return 5.0
        elif n_tm == min_tm - 1 or n_tm == max_tm + 1:
            return 3.0
        elif n_tm >= 2:
            return 1.0
        else:
            return 0.0
    
    @staticmethod
    def score_pore_architecture(has_pore_loop: bool,
                                has_selectivity_filter: bool) -> float:
        """
        Score channel-like pore architecture (0-5 points).
        
        Parameters
        ----------
        has_pore_loop : bool
            Predicted pore-forming loop/re-entrant helix
        has_selectivity_filter : bool
            Acidic residues suggesting Ca²⁺ selectivity
            
        Returns
        -------
        score : float
            Score 0-5
        """
        score = 0.0
        
        if has_pore_loop:
            score += 3.0
        if has_selectivity_filter:
            score += 2.0
        
        return score


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
        else:
            return sorted(self.candidates, key=lambda x: getattr(x, f"{by_score}_score", 0), 
                        reverse=True)
    
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
