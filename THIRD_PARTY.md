# Third-party data and methods

Data and definitions copied from outside sources. Each was retrieved once by hand; no code or test in this repository
contacts these services.

## PROSITE PS00018 (EF_HAND_1): EF-hand motif
- Used as `MotifScanner.MOTIFS['EF_hand']` (src/sequence_analysis.py) and in web/static/seq_engine.js.
- Pattern, copied from the entry:
  `D-{W}-[DNS]-{ILVFYW}-[DENSTG]-[DNQGHRK]-{GP}-[LIVMC]-[DENQSTAGC]-x(2)-[DE]-[LIVMFYW]`.
- Source: https://prosite.expasy.org/PS00018 (entry dated 29-MAY-2024 INFO UPDATE), retrieved 2026-11-12.
- Licence: PROSITE is distributed under CC BY-NC-ND 4.0 (SIB), https://prosite.expasy.org/prosite_license.html; only the pattern string is reproduced here, with attribution.

## Kyte & Doolittle hydropathy: TM prediction
- Kyte J, Doolittle RF (1982). A simple method for displaying the hydropathic character of a protein.
  *J. Mol. Biol.* 157(1):105-132. doi:10.1016/0022-2836(82)90515-0 (DOI checked on Crossref, 2026-11-12).
- Used for the hydropathy scale and the rule in `TransmembraneDomainPredictor.predict_tm_domains`: a 19-residue
  window with mean hydropathy >= 1.6 marks a probable membrane-spanning segment.

## UniProtKB: TM benchmark proteins
- `tests/fixtures/uniprot_tm_benchmark.json`: sequences and "Transmembrane" features of P02724 (GLPA_HUMAN),
  P02945 (BACR_HALSA), P29972 (AQP1_HUMAN) and P0AEJ4 (ENVZ_ECOLI), plus soluble decoys P0DP23 (CALM1_HUMAN) and
  P02144 (MYG_HUMAN).
- Source: UniProtKB REST API (https://rest.uniprot.org/uniprotkb/<accession>.json), retrieved 2026-11-12.
- Licence: CC BY 4.0 (https://www.uniprot.org/help/license).

## PubChem: ligand SMILES (`LigandLibrary.LIGANDS`, src/docking.py)
Fetched once on 2026-11-12 with the PUG REST API,
`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/<name>/property/Title,MolecularFormula,SMILES,ConnectivitySMILES/JSON`.
The library uses PubChem's isomeric `SMILES` (stereochemistry matters for docking); the canonical connectivity SMILES
is listed for reference. PubChem data-use policy: https://www.ncbi.nlm.nih.gov/home/about/policies/.

| Ligand | Query name | CID | PubChem title | Formula | Canonical (connectivity) SMILES |
|---|---|---|---|---|---|
| InsP3 | inositol 1,4,5-trisphosphate | 439456 | Ins(1,4,5)P3 | C6H15O15P3 | `C1(C(C(C(C(C1OP(=O)(O)O)O)OP(=O)(O)O)OP(=O)(O)O)O)O` |
| InsP4 | inositol 1,3,4,5-tetrakisphosphate | 107758 | Inositol 1,3,4,5-tetrakisphosphate | C6H16O18P4 | `C1(C(C(C(C(C1OP(=O)(O)O)OP(=O)(O)O)OP(=O)(O)O)O)OP(=O)(O)O)O` |
| cADPR | cyclic ADP-ribose | 123847 | Cyclic ADP-Ribose | C15H21N5O13P2 | `C1C2C(C(C(O2)N3C=NC4=C3N=CN(C4=N)C5C(C(C(O5)COP(=O)(OP(=O)(O1)O)O)O)O)O)O` |
| NAADP | nicotinic acid adenine dinucleotide phosphate | 123952 | Nicotinic Acid Adenine Dinucleotide Phosphate | C21H27N6O18P3 | `C1=CC(=C[N+](=C1)C2C(C(C(O2)COP(=O)([O-])OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)OP(=O)(O)O)O)O)O)C(=O)O` |
| NAD | NAD+ | 5892 | Nadide | C21H27N7O14P2 | `C1=CC(=C[N+](=C1)C2C(C(C(O2)COP(=O)([O-])OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O)C(=O)N` |

### 8-NH2-cADPR (8-amino-cyclic ADP-ribose)
Checked 2026-11-12 against the literature and three databases. The library SMILES is the neutral form of the compound
recorded in ChEMBL; it is cADPR (CID 123847) with adenine C8-H replaced by C8-NH2, and nothing else changed.

- Primary literature: Walseth TF, Lee HC (1993). Synthesis and characterization of antagonists of
  cyclic-ADP-ribose-induced Ca2+ release. *Biochim. Biophys. Acta* 1178(3):235-242. doi:10.1016/0167-4889(93)90199-Y
  (PMID 8395888). The abstract reports that the 1H-NMR spectrum of 8-amino-cADPR "showed characteristics of that of
  cADPR and confirmed the replacement of the 8-proton", i.e. the cyclic cADPR skeleton with an 8-substituent.
- Structure of the parent: Lee HC, Aarhus R, Levitt D (1994). The crystal structure of cyclic ADP-ribose.
  *Nat. Struct. Biol.* 1:143-144. doi:10.1038/nsb0394-143 (cyclisation through adenine N1). ChEBI CHEBI:31445
  ("cyclic ADP-beta-D-ribose") has the same InChIKey as PubChem CID 123847 (BQOHYSXSASDCEA-KEOHHSTQSA-N), so the
  cADPR stereocentres used here (both ribosides beta-D) are ChEBI's too.
- Database record: ChEMBL **CHEMBL3350988** (deposited in PubChem as CID **118719244**), from Zhang B et al. (2008)
  2'-Deoxy cyclic adenosine 5'-diphosphate ribose derivatives: importance of the 2'-hydroxyl motif for the antagonistic
  activity of 8-substituted cADPR derivatives. *J. Med. Chem.* 51:1623-1636. doi:10.1021/jm7010386 (ChEMBL document
  CHEMBL1142301; antagonist IC50 105 nM in sea urchin egg homogenate). ChEMBL draws it as the zwitterion (N1+ with a
  6-amino group, one phosphate O-). Moving those two charges back (N1=C6 to C6=NH, O- to OH) gives the library SMILES
  exactly: both have InChIKey DLAYKDRCOMDICA-AIKIOGRUSA-N (RDKit 2026.09.2), with all 8 sugar stereocentres equal.
- Formula C15H22N6O13P2, monoisotopic mass 556.0720 (cADPR: C15H21N5O13P2, 541.0611; the difference is NH, as for an
  amino group replacing H). Five independent rings: purine (2), two riboses, and the pyrophosphate macrocycle.
- Not used:
  - PubChem CID 135600320, the record found by the names "8-amino-cADPR" and "8-amino-cyclic ADP-ribose" (CAS
    151898-25-8). It has the right formula but a non-cyclic structure: an open-chain pyrophosphate and enol sugars.
    Its CAS-style synonym, like cADPR's own CAS name, says "N-alpha-D-ribofuranosyl"; that is the CAS index name and
    predates the crystal structure, which showed the N1 riboside is beta.
  - CIDs 10076390 and 22489536 (untitled, same formula): sugar stereocentres inverted relative to cADPR.
  - ChEBI has no entry for 8-amino-cADPR (searched by name and InChIKey). The IUPHAR/BPS Guide to Pharmacology
    could not be checked: its search and web services now require a registered account.
- `tests/test_ligand_structures.py` checks the formula with hydrogens, the ring count, that both P atoms lie in the
  macrocycle, the single 8-amino carbon, that removing the amino group gives the cADPR SMILES exactly, and (when RDKit
  is installed) the InChIKey.

The searches by "NAADP" (CID 52942301, cationic form) and "NAD+" (also CIDs 5893, 15938971, ...) returned several
records; the neutral/zwitterionic parent records above were chosen.
