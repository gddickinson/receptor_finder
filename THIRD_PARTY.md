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
is listed for reference. PubChem data is public domain (https://pubchemdocs.ncbi.nlm.nih.gov/downloads).

| Ligand | Query name | CID | PubChem title | Formula | Canonical (connectivity) SMILES |
|---|---|---|---|---|---|
| InsP3 | inositol 1,4,5-trisphosphate | 439456 | Ins(1,4,5)P3 | C6H15O15P3 | `C1(C(C(C(C(C1OP(=O)(O)O)O)OP(=O)(O)O)OP(=O)(O)O)O)O` |
| InsP4 | inositol 1,3,4,5-tetrakisphosphate | 107758 | Inositol 1,3,4,5-tetrakisphosphate | C6H16O18P4 | `C1(C(C(C(C(C1OP(=O)(O)O)OP(=O)(O)O)OP(=O)(O)O)O)OP(=O)(O)O)O` |
| cADPR | cyclic ADP-ribose | 123847 | Cyclic ADP-Ribose | C15H21N5O13P2 | `C1C2C(C(C(O2)N3C=NC4=C3N=CN(C4=N)C5C(C(C(O5)COP(=O)(OP(=O)(O1)O)O)O)O)O)O` |
| NAADP | nicotinic acid adenine dinucleotide phosphate | 123952 | Nicotinic Acid Adenine Dinucleotide Phosphate | C21H27N6O18P3 | `C1=CC(=C[N+](=C1)C2C(C(C(O2)COP(=O)([O-])OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)OP(=O)(O)O)O)O)O)C(=O)O` |
| NAD | NAD+ | 5892 | Nadide | C21H27N7O14P2 | `C1=CC(=C[N+](=C1)C2C(C(C(O2)COP(=O)([O-])OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O)C(=O)N` |

8-NH2-cADPR: PubChem's only record for this compound (CID 135600320, query "8-amino-cADPR", CAS 151898-25-8) has a
non-cyclic structure with enol sugars (its title is "8-Aminoadenosine cyclic 3',5'-(hydrogen phosphate)
5'-ribofuranosyl ester"), so it was not used. The library SMILES is cADPR (CID 123847) with the adenine C8-H replaced
by C8-NH2 (`N3C=NC4` -> `N3C(N)=NC4`). Its formula, C15H22N6O13P2 (checked with RDKit 2026.09.2), matches CID
135600320. **George: please confirm this structure or supply a preferred source.**

The searches by "NAADP" (CID 52942301, cationic form) and "NAD+" (also CIDs 5893, 15938971, ...) returned several
records; the neutral/zwitterionic parent records above were chosen.
