"""
Binding Pocket Module
=====================

Detect and characterize ligand binding pockets in protein structures.
Provides tools for pocket identification, scoring, and comparison to known receptors.
"""

import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple
from pathlib import Path


@dataclass
class BindingPocket:
    """Data class for binding pocket."""
    
    pocket_id: str
    protein_id: str
    
    # Geometric properties
    center_coords: Optional[Tuple[float, float, float]] = None
    volume: Optional[float] = None  # Å³
    depth: Optional[float] = None   # Å
    
    # Chemical properties
    charge: float = 0.0
    hydrophobicity: float = 0.0
    
    # Residue composition
    residues: Optional[List[str]] = None
    aromatic_count: int = 0
    basic_count: int = 0
    acidic_count: int = 0
    polar_count: int = 0
    
    # Ligand-specific scores
    insp3_score: Optional[float] = None
    cadpr_score: Optional[float] = None
    
    def __post_init__(self):
        if self.residues:
            self._count_residue_types()
    
    def _count_residue_types(self):
        """Count residue types in pocket."""
        aromatic = {'PHE', 'TRP', 'TYR'}
        basic = {'ARG', 'LYS', 'HIS'}
        acidic = {'ASP', 'GLU'}
        polar = {'SER', 'THR', 'ASN', 'GLN', 'CYS'}
        
        self.aromatic_count = sum(1 for r in self.residues if r in aromatic)
        self.basic_count = sum(1 for r in self.residues if r in basic)
        self.acidic_count = sum(1 for r in self.residues if r in acidic)
        self.polar_count = sum(1 for r in self.residues if r in polar)


class PocketDetector:
    """
    Detect binding pockets in protein structures.
    
    In practice, would use tools like:
    - fpocket
    - P2Rank
    - SiteMap (Schrödinger)
    - DEPTH
    
    This is a simplified implementation for demonstration.
    """
    
    @staticmethod
    def detect_pockets(pdb_file: str, protein_id: str,
                      min_volume: float = 200.0) -> List[BindingPocket]:
        """
        Detect binding pockets in structure.
        
        NOTE: This is a placeholder. Real implementation would:
        - Use fpocket or P2Rank
        - Calculate grid-based pockets
        - Identify surface cavities
        
        Parameters
        ----------
        pdb_file : str
            Path to PDB file
        protein_id : str
            Protein identifier
        min_volume : float
            Minimum pocket volume (Å³)
            
        Returns
        -------
        pockets : List[BindingPocket]
            Detected pockets
        """
        # Placeholder: Generate mock pockets
        pockets = []
        
        for i in range(3):  # Simulate 3 pockets
            pocket = BindingPocket(
                pocket_id=f"{protein_id}_POCKET_{i+1}",
                protein_id=protein_id,
                center_coords=(10.0 + i*5, 10.0, 10.0 + i*5),
                volume=250.0 + i*50,
                depth=8.0 + i*2,
                charge=2.0 - i,
                hydrophobicity=0.3 + i*0.2,
                residues=['ARG', 'LYS', 'ASP', 'PHE', 'TRP'] if i == 0 
                        else ['LEU', 'VAL', 'ILE', 'ALA'] if i == 1
                        else ['SER', 'THR', 'ASN', 'GLN']
            )
            pockets.append(pocket)
        
        return pockets
    
    @staticmethod
    def get_fpocket_command(pdb_file: str) -> str:
        """
        Generate fpocket command (for user reference).
        
        fpocket is a widely-used pocket detection tool.
        
        Example:
        ```bash
        fpocket -f protein.pdb
        ```
        """
        return f"fpocket -f {pdb_file}"


class PocketScorer:
    """Score binding pockets for InsP₃ and cADPR binding potential."""
    
    @staticmethod
    def score_insp3_binding(pocket: BindingPocket) -> float:
        """
        Score pocket for InsP₃ binding potential.
        
        InsP₃ binding requires:
        - Positive electrostatic potential (for 4,5-bisphosphate)
        - Sufficient volume (>200 Å³)
        - Basic residues (Arg, Lys)
        - Aromatic residues (for inositol ring)
        
        Returns score 0-13 (higher is better)
        """
        score = 0.0
        
        # Volume requirement (0-2 points)
        if pocket.volume and pocket.volume > 200:
            score += 2
        
        # Positive charge for phosphates (0-3 points)
        if pocket.charge > 2:
            score += 3
        elif pocket.charge > 0:
            score += 1.5
        
        # Basic residues (0-4 points)
        score += min(pocket.basic_count, 4)
        
        # Aromatic residues for inositol (0-2 points)
        score += min(pocket.aromatic_count, 2)
        
        # Accessibility (0-2 points)
        if pocket.depth and pocket.depth > 5:
            score += 2
        
        pocket.insp3_score = score
        return score
    
    @staticmethod
    def score_cadpr_binding(pocket: BindingPocket) -> float:
        """
        Score pocket for cADPR binding potential.
        
        cADPR binding requires:
        - Larger volume (>300 Å³) for cyclic structure
        - Some positive charge (pyrophosphate)
        - Hydrophobic character (for adenine)
        - Hydrogen bond donors (for ribose)
        
        Returns score 0-12 (higher is better)
        """
        score = 0.0
        
        # Volume requirement (0-2 points)
        if pocket.volume and pocket.volume > 300:
            score += 2
        
        # Moderate charge (0-2 points)
        if pocket.charge > 0:
            score += 2
        
        # Hydrophobic character (0-3 points)
        if pocket.hydrophobicity > 0.4:
            score += 3
        elif pocket.hydrophobicity > 0.2:
            score += 1.5
        
        # Polar residues for H-bonding (0-3 points)
        score += min(pocket.polar_count, 3)
        
        # Aromatic for adenine (0-2 points)
        score += min(pocket.aromatic_count, 2)
        
        pocket.cadpr_score = score
        return score
    
    @staticmethod
    def compare_to_reference(pocket: BindingPocket,
                           reference_residues: List[str]) -> float:
        """
        Compare pocket residues to reference binding site.
        
        Parameters
        ----------
        pocket : BindingPocket
            Query pocket
        reference_residues : List[str]
            Residues from known binding site
            
        Returns
        -------
        similarity : float
            Jaccard similarity (0-1)
        """
        if not pocket.residues or not reference_residues:
            return 0.0
        
        set_pocket = set(pocket.residues)
        set_ref = set(reference_residues)
        
        intersection = len(set_pocket & set_ref)
        union = len(set_pocket | set_ref)
        
        if union == 0:
            return 0.0
        
        return intersection / union


class PocketCharacterizer:
    """Detailed characterization of binding pockets."""
    
    # Known InsP₃ receptor binding site residues (from animal InsP₃R)
    INSP3R_BINDING_SITE = ['ARG', 'ARG', 'LYS', 'ARG', 'PHE', 'TRP', 
                          'ASP', 'SER', 'THR', 'TYR']
    
    # Known ryanodine receptor residues
    RYR_BINDING_SITE = ['ARG', 'LYS', 'GLU', 'ASP', 'PHE', 'LEU',
                       'SER', 'THR', 'HIS']
    
    @staticmethod
    def create_pharmacophore(pocket: BindingPocket) -> Dict:
        """
        Create pharmacophore model from pocket.
        
        Pharmacophore defines spatial arrangement of features
        important for binding.
        
        Returns
        -------
        pharmacophore : Dict
            Pharmacophore features
        """
        features = {
            'positive_charge_centers': [],
            'aromatic_centers': [],
            'hbond_donors': [],
            'hbond_acceptors': [],
            'hydrophobic_centers': []
        }
        
        # Simplified: count feature types
        if pocket.residues:
            for residue in pocket.residues:
                if residue in ['ARG', 'LYS', 'HIS']:
                    features['positive_charge_centers'].append(residue)
                elif residue in ['PHE', 'TRP', 'TYR']:
                    features['aromatic_centers'].append(residue)
                elif residue in ['SER', 'THR', 'TYR', 'CYS']:
                    features['hbond_donors'].append(residue)
                elif residue in ['ASP', 'GLU', 'ASN', 'GLN']:
                    features['hbond_acceptors'].append(residue)
                elif residue in ['LEU', 'VAL', 'ILE', 'PHE', 'TRP']:
                    features['hydrophobic_centers'].append(residue)
        
        return features
    
    @staticmethod
    def assess_specificity(pocket: BindingPocket,
                          insp3_score: float,
                          cadpr_score: float,
                          insp4_score: Optional[float] = None) -> Dict:
        """
        Assess ligand specificity.
        
        Important: receptor should bind InsP₃ but NOT InsP₄
        (or cADPR but not NAD⁺).
        
        Returns
        -------
        assessment : Dict
            Specificity assessment
        """
        assessment = {
            'primary_ligand': None,
            'is_specific': False,
            'specificity_ratio': None
        }
        
        # Determine primary ligand
        if insp3_score > cadpr_score:
            assessment['primary_ligand'] = 'InsP3'
            
            # Check specificity (InsP₃ > InsP₄)
            if insp4_score is not None:
                if insp3_score > insp4_score + 2:  # At least 2 points difference
                    assessment['is_specific'] = True
                    assessment['specificity_ratio'] = insp3_score / (insp4_score + 0.1)
        else:
            assessment['primary_ligand'] = 'cADPR'
            assessment['is_specific'] = True  # Assume specific unless tested otherwise
        
        return assessment


class PocketDatabase:
    """Manage collection of detected pockets."""
    
    def __init__(self):
        self.pockets: Dict[str, BindingPocket] = {}
    
    def add_pocket(self, pocket: BindingPocket):
        """Add pocket to database."""
        self.pockets[pocket.pocket_id] = pocket
    
    def get_pocket(self, pocket_id: str) -> Optional[BindingPocket]:
        """Retrieve pocket by ID."""
        return self.pockets.get(pocket_id)
    
    def get_pockets_by_protein(self, protein_id: str) -> List[BindingPocket]:
        """Get all pockets for a protein."""
        return [p for p in self.pockets.values() if p.protein_id == protein_id]
    
    def get_high_scoring_pockets(self, ligand: str = 'InsP3',
                                min_score: float = 7.0) -> List[BindingPocket]:
        """Get pockets with high binding scores."""
        if ligand == 'InsP3':
            return [p for p in self.pockets.values() 
                   if p.insp3_score and p.insp3_score >= min_score]
        elif ligand == 'cADPR':
            return [p for p in self.pockets.values()
                   if p.cadpr_score and p.cadpr_score >= min_score]
        else:
            return []
    
    def to_dataframe(self) -> pd.DataFrame:
        """Convert to DataFrame."""
        data = []
        for pocket in self.pockets.values():
            data.append({
                'pocket_id': pocket.pocket_id,
                'protein_id': pocket.protein_id,
                'volume': pocket.volume,
                'depth': pocket.depth,
                'charge': pocket.charge,
                'basic_count': pocket.basic_count,
                'aromatic_count': pocket.aromatic_count,
                'insp3_score': pocket.insp3_score,
                'cadpr_score': pocket.cadpr_score
            })
        return pd.DataFrame(data)
    
    def save_summary(self, filepath: str):
        """Save summary to CSV."""
        df = self.to_dataframe()
        df.to_csv(filepath, index=False)
    
    def __len__(self):
        return len(self.pockets)


if __name__ == "__main__":
    # Test code
    print("Testing Binding Pocket Module\n")
    
    # Test 1: Pocket detection
    print("Test 1: Pocket Detection")
    detector = PocketDetector()
    pockets = detector.detect_pockets("test.pdb", "TEST_PROTEIN")
    
    print(f"  Detected {len(pockets)} pockets")
    for pocket in pockets:
        print(f"    {pocket.pocket_id}: volume={pocket.volume:.1f} Ų, "
              f"charge={pocket.charge:.1f}")
    
    # Test 2: InsP₃ binding scoring
    print("\nTest 2: InsP₃ Binding Scoring")
    scorer = PocketScorer()
    
    for pocket in pockets:
        score = scorer.score_insp3_binding(pocket)
        print(f"  {pocket.pocket_id}: InsP₃ score = {score:.1f}/13")
    
    # Test 3: cADPR binding scoring
    print("\nTest 3: cADPR Binding Scoring")
    for pocket in pockets:
        score = scorer.score_cadpr_binding(pocket)
        print(f"  {pocket.pocket_id}: cADPR score = {score:.1f}/12")
    
    # Test 4: Pharmacophore generation
    print("\nTest 4: Pharmacophore Analysis")
    characterizer = PocketCharacterizer()
    pharmacophore = characterizer.create_pharmacophore(pockets[0])
    
    print(f"  Pocket: {pockets[0].pocket_id}")
    for feature, residues in pharmacophore.items():
        if residues:
            print(f"    {feature}: {len(residues)} sites")
    
    # Test 5: Pocket database
    print("\nTest 5: Pocket Database")
    db = PocketDatabase()
    for pocket in pockets:
        db.add_pocket(pocket)
    
    print(f"  Database size: {len(db)}")
    high_scoring = db.get_high_scoring_pockets('InsP3', min_score=5.0)
    print(f"  High-scoring (InsP₃ > 5): {len(high_scoring)}")
    
    # Show summary
    df = db.to_dataframe()
    print("\n  Summary:")
    print(df[['pocket_id', 'volume', 'charge', 'insp3_score', 'cadpr_score']])
    
    print("\nTest complete!")
