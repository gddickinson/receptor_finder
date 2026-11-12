"""Structure checks for LigandLibrary SMILES beyond heavy-atom counts (sources in THIRD_PARTY.md).

8-NH2-cADPR (8-amino-cyclic ADP-ribose, Walseth & Lee 1993) must be cADPR with an amino group on adenine C8: cyclic
(the pyrophosphate closes a macrocycle through the adenine N1-ribose), C15H22N6O13P2, and every stereocentre as in
cADPR. PubChem's record under that name (CID 135600320) is a non-cyclic structure; ChEMBL CHEMBL3350988 (from Zhang et
al. 2008, J. Med. Chem. 51:1623) is the cyclic compound, drawn as the zwitterion.
"""
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import smiles_graph as G  # noqa: E402
from docking import LigandLibrary  # noqa: E402

AMINO = LigandLibrary.LIGANDS["8-NH2-cADPR"]
CADPR = LigandLibrary.LIGANDS["cADPR"]
# PubChem CID 135600320 ("8-Aminoadenosine cyclic 3',5'-(hydrogen phosphate) 5'-ribofuranosyl ester"), retrieved
# 2026-11-12: the record the library must NOT use (open-chain pyrophosphate, enol sugars).
PUBCHEM_135600320 = ("C1NC2=C(C(=N)N1[C@H]3C(=C([C@H](O3)CO)O)O)N=C(N2[C@H]4C(=C[C@H](O4)COP(=O)(O)OP(=O)(O)O)O)N")
# ChEMBL CHEMBL3350988 standard InChIKey (zwitterion) and the InChIKey of its neutral form, computed with RDKit 2026.09.2.
CHEMBL_NEUTRAL_INCHIKEY = "DLAYKDRCOMDICA-AIKIOGRUSA-N"


def _c8_amino_carbons(smiles):
    """Carbons bonded to three nitrogens, one a terminal NH2 on a single bond: the 8-amino purine C8 (adenine C2 has
    two N neighbours, C6 carries =NH)."""
    atoms, bonds = G.parse(smiles)
    out = []
    for i, a in enumerate(atoms):
        nb = G.neighbours(smiles, i)
        if a["el"] != "C" or [atoms[n]["el"] for n in nb].count("N") != 3:
            continue
        if any(len(G.neighbours(smiles, n)) == 1 and bonds.get((i, n), bonds.get((n, i))) == 1 for n in nb):
            out.append(i)
    return out


class TestEightAminoCadpr(unittest.TestCase):

    def test_formula(self):
        self.assertEqual(G.formula(AMINO), {"C": 15, "H": 22, "N": 6, "O": 13, "P": 2})
        self.assertEqual(G.formula(CADPR), {"C": 15, "H": 21, "N": 5, "O": 13, "P": 2})

    def test_cyclic_like_cadpr(self):
        self.assertEqual(G.ring_count(AMINO), G.ring_count(CADPR))
        self.assertEqual(G.ring_count(AMINO), 5)  # purine (2) + 2 riboses + the macrocycle
        self.assertEqual(G.ring_count(PUBCHEM_135600320), 4)
        atoms, bonds = G.parse(AMINO)
        p_bonds = [(i, j) for i, j in bonds if "P" in (atoms[i]["el"], atoms[j]["el"])]
        on_ring = [b for b in p_bonds if G.in_ring(AMINO, *b)]
        self.assertEqual(len(on_ring), 4, "both P atoms must sit in the macrocycle (P-O-P plus the two ester O)")

    def test_eight_amino_group(self):
        self.assertEqual(len(_c8_amino_carbons(AMINO)), 1)
        self.assertEqual(_c8_amino_carbons(CADPR), [])

    def test_only_difference_from_cadpr_is_c8_amino(self):
        # Same atom order and stereo marks as cADPR (CID 123847) once the 8-amino branch is removed.
        self.assertEqual(AMINO.replace("N3C(N)=NC4", "N3C=NC4"), CADPR)
        self.assertEqual(AMINO.count("@"), CADPR.count("@"))

    def test_matches_chembl_record(self):
        try:
            from rdkit import Chem
        except ImportError:
            self.skipTest("RDKit not installed")
        self.assertEqual(Chem.MolToInchiKey(Chem.MolFromSmiles(AMINO)), CHEMBL_NEUTRAL_INCHIKEY)

    def test_source_recorded(self):
        self.assertIn("CHEMBL3350988", LigandLibrary.LIGAND_SOURCES["8-NH2-cADPR"])


if __name__ == "__main__":
    unittest.main()
