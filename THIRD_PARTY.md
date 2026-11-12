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
