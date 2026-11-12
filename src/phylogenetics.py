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

    Patterns follow docs/in_silico_calcium_projects.md (section 3.1) and
    the names scored in scoring.EvolutionaryScorer:

    - Conserved: in the alga and in angiosperms and early land plants
      (bryophytes or the lycophyte).
    - Chlamydomonas_loss: in the alga (which has a canonical InsP3R), absent
      from early land plants, present in angiosperms ("present in algae,
      lost in land plants, re-expanded in angiosperms").
    - Land_plant_innovation: absent from the alga, present in bryophytes
      and angiosperms ("absent in algae, appeared in bryophytes, conserved
      in all land plants").
    - Angiosperm_specific: present only in flowering plants.
    - Bryophyte_origin: in bryophytes, absent from the alga and angiosperms.
    - Sporadic: present somewhere but in none of the patterns above.
    - Unknown: absent from every species sampled.

    Presence/absence alone cannot separate a loss from a gain; the labels
    are hypotheses to test with a gene tree.

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

    def _any(clades):
        return any(presence.get(sp, False) for sp, info in PLANT_SPECIES.items()
                   if info["clade"] in clades)

    has_chlamy = presence.get("Chlamydomonas_reinhardtii", False)
    has_bryophyte = _any(("Bryophytes",))
    has_lycophyte = _any(("Lycophytes",))
    has_angiosperm = _any(("Eudicots", "Monocots", "Basal_angiosperm"))
    has_early_land = has_bryophyte or has_lycophyte

    if not present:
        pattern, confidence = "Unknown", 0.0
    elif has_chlamy and has_angiosperm and has_early_land:
        pattern, confidence = "Conserved", 0.8
    elif has_chlamy and has_angiosperm:
        pattern, confidence = "Chlamydomonas_loss", 0.9
    elif not has_chlamy and has_bryophyte and has_angiosperm:
        pattern, confidence = "Land_plant_innovation", 0.8
    elif has_angiosperm and not has_chlamy and not has_early_land:
        pattern, confidence = "Angiosperm_specific", 0.7
    elif has_bryophyte and not has_chlamy and not has_angiosperm:
        pattern, confidence = "Bryophyte_origin", 0.6
    else:
        pattern, confidence = "Sporadic", 0.5

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
