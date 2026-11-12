"""
Scorers
=======

Individual scoring rules for receptor candidates (structural, evolutionary,
expression, localization and topology evidence). Each returns points on the
scale given in its docstring; scoring.CandidateScore adds them up.
Split out of scoring.py to keep files under 500 lines.
"""

from typing import List, Optional, Tuple


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
            return 10.0 - (binding_affinity + 9.0) * (3.0 / 2.0)
        elif binding_affinity < -5.0:
            # Linear interpolation: -7 to -5 maps to 7 to 3
            return 7.0 - (binding_affinity + 7.0) * (4.0 / 2.0)
        else:
            # >= -5.0: scale from 3 (at -5) down to 0 (at 0 kcal/mol)
            return max(0.0, 3.0 - (binding_affinity + 5.0) * 0.6)
    
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

    PATTERN_SCORES = {
        "Chlamydomonas_loss": 10.0,
        "Land_plant_innovation": 8.0,
        "Angiosperm_specific": 6.0,
        "Bryophyte_origin": 7.0,
        "Conserved": 4.0,
        "Sporadic": 2.0,
        "Unknown": 0.0
    }
    
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
        return EvolutionaryScorer.PATTERN_SCORES.get(pattern, 0.0)
    
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

    LOCALIZATION_SCORES = {
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
    SIGNAL_TYPE_SCORES = {
        "ER": 5.0,
        "vacuolar": 5.0,
        "secretory": 4.0,
        "mitochondrial": 1.0,
        "chloroplast": 1.0
    }
    
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
        table = {k.lower(): v for k, v in LocalizationScorer.LOCALIZATION_SCORES.items()}
        return table.get(predicted_location.lower(), 0.0)
    
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
            return LocalizationScorer.SIGNAL_TYPE_SCORES.get(signal_type, 2.0)
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
