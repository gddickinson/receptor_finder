"""
Phylogenetics Module
====================

Gene family expansion analysis and taxonomic distribution patterns
for receptor candidate evaluation.
"""

import numpy as np
import pandas as pd
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple
from pathlib import Path


# Representative plant species for phylogenetic analysis
PLANT_SPECIES = {
    "Arabidopsis_thaliana": {"clade": "Eudicots", "group": "Land_plants"},
    "Oryza_sativa": {"clade": "Monocots", "group": "Land_plants"},
    "Zea_mays": {"clade": "Monocots", "group": "Land_plants"},
    "Marchantia_polymorpha": {"clade": "Bryophytes", "group": "Land_plants"},
    "Physcomitrella_patens": {"clade": "Bryophytes", "group": "Land_plants"},
    "Chlamydomonas_reinhardtii": {"clade": "Chlorophyte", "group": "Algae"},
    "Selaginella_moellendorffii": {"clade": "Lycophytes", "group": "Land_plants"},
    "Amborella_trichopoda": {"clade": "Basal_angiosperm", "group": "Land_plants"},
}


@dataclass
class GeneFamily:
    """Represents a gene family across species."""

    family_id: str
    description: str = ""
    copy_numbers: Dict[str, int] = field(default_factory=dict)

    @property
    def total_copies(self) -> int:
        """Total copies across all species."""
        return sum(self.copy_numbers.values())

    @property
    def species_count(self) -> int:
        """Number of species that have at least one copy."""
        return sum(1 for v in self.copy_numbers.values() if v > 0)

    def expansion_factor(self, species: str, ancestral: int = 1) -> float:
        """Compute expansion factor for a given species."""
        copies = self.copy_numbers.get(species, 0)
        return copies / max(ancestral, 1)


@dataclass
class PhylogeneticPattern:
    """Classification of a gene's phylogenetic distribution."""

    gene_id: str
    pattern: str  # e.g. "Chlamydomonas_loss", "Land_plant_innovation"
    present_in: List[str] = field(default_factory=list)
    absent_in: List[str] = field(default_factory=list)
    confidence: float = 0.0  # 0-1


def classify_distribution(gene_id: str,
                          presence: Dict[str, bool]) -> PhylogeneticPattern:
    """
    Classify a gene's phylogenetic distribution pattern.

    Parameters
    ----------
    gene_id : str
        Gene identifier.
    presence : Dict[str, bool]
        Mapping of species name to presence (True/False).

    Returns
    -------
    PhylogeneticPattern
        The classified pattern.
    """
    present = [sp for sp, val in presence.items() if val]
    absent = [sp for sp, val in presence.items() if not val]

    has_chlamy = presence.get("Chlamydomonas_reinhardtii", False)
    has_land = any(
        presence.get(sp, False)
        for sp, info in PLANT_SPECIES.items()
        if info["group"] == "Land_plants"
    )
    has_bryophyte = any(
        presence.get(sp, False)
        for sp, info in PLANT_SPECIES.items()
        if info["clade"] == "Bryophytes"
    )
    has_angiosperm = any(
        presence.get(sp, False)
        for sp, info in PLANT_SPECIES.items()
        if info["clade"] in ("Eudicots", "Monocots", "Basal_angiosperm")
    )

    # Determine pattern
    if has_land and not has_chlamy:
        pattern = "Chlamydomonas_loss"
        confidence = 0.9
    elif has_land and has_chlamy:
        pattern = "Conserved"
        confidence = 0.8
    elif has_angiosperm and not has_bryophyte:
        pattern = "Angiosperm_specific"
        confidence = 0.7
    elif has_bryophyte and has_angiosperm:
        pattern = "Land_plant_innovation"
        confidence = 0.8
    elif has_bryophyte and not has_angiosperm:
        pattern = "Bryophyte_origin"
        confidence = 0.6
    elif len(present) <= 1:
        pattern = "Sporadic"
        confidence = 0.5
    else:
        pattern = "Unknown"
        confidence = 0.0

    return PhylogeneticPattern(
        gene_id=gene_id,
        pattern=pattern,
        present_in=present,
        absent_in=absent,
        confidence=confidence,
    )


def build_copy_number_table(families: List[GeneFamily]) -> pd.DataFrame:
    """
    Build a species-by-family copy number table.

    Returns DataFrame with species as rows and families as columns.
    """
    all_species = set()
    for fam in families:
        all_species.update(fam.copy_numbers.keys())

    data = {}
    for fam in families:
        data[fam.family_id] = {
            sp: fam.copy_numbers.get(sp, 0) for sp in all_species
        }

    return pd.DataFrame(data)


if __name__ == "__main__":
    print("Testing Phylogenetics Module\n")

    # Test 1: classify distribution
    presence = {
        "Arabidopsis_thaliana": True,
        "Oryza_sativa": True,
        "Marchantia_polymorpha": True,
        "Chlamydomonas_reinhardtii": False,
    }
    pat = classify_distribution("AT1G12345", presence)
    print(f"Pattern: {pat.pattern} (confidence {pat.confidence:.1f})")

    # Test 2: gene family
    fam = GeneFamily(
        family_id="TPC_family",
        description="Two-pore channel family",
        copy_numbers={"Arabidopsis_thaliana": 2, "Oryza_sativa": 3, "Zea_mays": 4},
    )
    print(f"Family {fam.family_id}: {fam.total_copies} total copies")
    print(f"Expansion in maize: {fam.expansion_factor('Zea_mays'):.1f}x")

    print("\nTest complete!")
