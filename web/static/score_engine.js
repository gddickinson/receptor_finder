// Pure JS port of the scoring rules in src/scoring.py, src/binding_pocket.py (PocketScorer, assess_specificity),
// src/phylogenetics.py (classify_distribution) and src/docking.py (DockingAnalyzer.assess_specificity).
// No DOM. Node-tested against the Python code (score_engine.test.mjs, tests/fixtures/engine_fixture.json).

export const PATTERN_SCORES = { Chlamydomonas_loss: 10, Land_plant_innovation: 8, Angiosperm_specific: 6,
  Bryophyte_origin: 7, Conserved: 4, Sporadic: 2, Unknown: 0 };
export const LOCALIZATION_SCORES = { tonoplast: 10, vacuole: 10, ER: 9, endoplasmic_reticulum: 9, plasma_membrane: 6,
  endomembrane: 7, golgi: 5, mitochondria: 2, chloroplast: 2, nucleus: 1, cytosol: 1 };
export const SIGNAL_TYPE_SCORES = { ER: 5, vacuolar: 5, secretory: 4, mitochondrial: 1, chloroplast: 1 };
export const PLANT_SPECIES = {
  Arabidopsis_thaliana: ['Eudicots', 'Land_plants'], Oryza_sativa: ['Monocots', 'Land_plants'],
  Zea_mays: ['Monocots', 'Land_plants'], Marchantia_polymorpha: ['Bryophytes', 'Land_plants'],
  Physcomitrella_patens: ['Bryophytes', 'Land_plants'], Chlamydomonas_reinhardtii: ['Chlorophyte', 'Algae'],
  Selaginella_moellendorffii: ['Lycophytes', 'Land_plants'], Amborella_trichopoda: ['Basal_angiosperm', 'Land_plants'],
};
// CandidateScore fields: [name, max points, category].
export const FIELDS = [
  ['pocket_score', 13, 'Structural'], ['docking_score', 10, 'Structural'], ['structure_quality', 7, 'Structural'],
  ['phylo_pattern_score', 10, 'Evolutionary'], ['expansion_score', 10, 'Evolutionary'],
  ['coexpression_score', 10, 'Expression'], ['tissue_specificity', 10, 'Expression'],
  ['localization_score', 10, 'Localization'], ['signal_peptide_score', 5, 'Localization'],
  ['tm_domain_score', 5, 'Topology'], ['pore_architecture', 5, 'Topology'], ['topology_match', 5, 'Topology'],
];

export const scorePocket = (raw, max = 13) => Math.min(raw, max);

export function scoreDocking(a) {
  if (a < -9) return 10;
  if (a < -7) return 10 - (a + 9) * (3 / 2);
  if (a < -5) return 7 - (a + 7) * (4 / 2);
  return Math.max(0, 3 - (a + 5) * 0.6);
}

export function scoreStructureQuality(mean, tm = null) {
  let s;
  if (mean > 90) s = 7;
  else if (mean > 70) s = 5 + (mean - 70) * (2 / 20);
  else if (mean > 50) s = 2 + (mean - 50) * (3 / 20);
  else s = (mean / 50) * 2;
  if (tm !== null && tm !== undefined && tm < 70) s *= tm / 70;
  return Math.min(s, 7);
}

export const scorePhyloPattern = p => PATTERN_SCORES[p] ?? 0;

export function scoreGeneExpansion(copies, ancestral = 1) {
  const f = copies / Math.max(ancestral, 1);
  if (f >= 10) return 10;
  if (f >= 5) return 8;
  if (f >= 2) return 5 + (f - 2) * (3 / 3);
  return f * 2.5;
}

export function scoreCoexpression(r) {
  if (r > 0.8) return 10;
  if (r > 0.6) return 7 + (r - 0.6) * (3 / 0.2);
  if (r > 0.4) return 4 + (r - 0.4) * (3 / 0.2);
  if (r > 0.2) return 1 + (r - 0.2) * (3 / 0.2);
  return Math.max(0, r * 5);
}

export function scoreTissueSpecificity(tissues, expected) {
  if (!expected.length) return 5;
  const a = new Set(tissues), overlap = [...new Set(expected)].filter(t => a.has(t)).length;
  return (overlap / expected.length) * 10;
}

export function scoreLocalization(loc) {
  const key = String(loc).toLowerCase();
  for (const [k, v] of Object.entries(LOCALIZATION_SCORES)) if (k.toLowerCase() === key) return v;
  return 0;
}

export function scoreSignalPeptide(has, type = null) {
  if (!has) return 0;
  if (type) return SIGNAL_TYPE_SCORES[type] ?? 2;
  return 2;
}

export function scoreTmDomains(n, [lo, hi] = [4, 8]) {
  if (lo <= n && n <= hi) return 5;
  if (n === lo - 1 || n === hi + 1) return 3;
  return n >= 2 ? 1 : 0;
}

export const scorePoreArchitecture = (loop, filter) => (loop ? 3 : 0) + (filter ? 2 : 0);

export function totalScore(c) { return FIELDS.reduce((s, [k]) => s + (+c[k] || 0), 0); }
export function priorityClass(total) { return total >= 60 ? 'HIGH' : total >= 40 ? 'MEDIUM' : 'LOW'; }

// --- binding_pocket.py ---------------------------------------------------------------------------------------------
const AROMATIC = new Set(['PHE', 'TRP', 'TYR']), BASIC = new Set(['ARG', 'LYS', 'HIS']),
  ACIDIC = new Set(['ASP', 'GLU']), POLAR = new Set(['SER', 'THR', 'ASN', 'GLN', 'CYS']);
export function residueCounts(residues = []) {
  const n = set => residues.filter(r => set.has(r)).length;
  return { basic: n(BASIC), aromatic: n(AROMATIC), acidic: n(ACIDIC), polar: n(POLAR) };
}

export function scoreInsp3(p) {
  const c = residueCounts(p.residues || []);
  let s = 0;
  if (p.volume && p.volume > 200) s += 2;
  if (p.charge > 2) s += 3; else if (p.charge > 0) s += 1.5;
  s += Math.min(c.basic, 4) + Math.min(c.aromatic, 2);
  if (p.depth && p.depth > 5) s += 2;
  return s;
}

export function scoreCadpr(p) {
  const c = residueCounts(p.residues || []);
  let s = 0;
  if (p.volume && p.volume > 300) s += 2;
  if (p.charge > 0) s += 2;
  if (p.hydrophobicity > 0.4) s += 3; else if (p.hydrophobicity > 0.2) s += 1.5;
  return s + Math.min(c.polar, 3) + Math.min(c.aromatic, 2);
}

// PocketCharacterizer.assess_specificity (pocket-score units).
export function pocketSpecificity(insp3, cadpr, insp4 = null) {
  const a = { primary_ligand: null, is_specific: false, specificity_ratio: null };
  if (insp3 > cadpr) {
    a.primary_ligand = 'InsP3';
    if (insp4 !== null && insp4 !== undefined && insp3 > insp4 + 2) {
      a.is_specific = true; a.specificity_ratio = insp3 / (insp4 + 0.1);
    }
  } else { a.primary_ligand = 'cADPR'; a.is_specific = true; }
  return a;
}

// --- phylogenetics.py ------------------------------------------------------------------------------------------------
export function classifyDistribution(presence) {
  const present = Object.keys(presence).filter(k => presence[k]);
  const any = test => Object.entries(PLANT_SPECIES).some(([sp, [clade, group]]) => test(clade, group) && presence[sp]);
  const chlamy = !!presence.Chlamydomonas_reinhardtii, land = any((c, g) => g === 'Land_plants'),
    bryo = any(c => c === 'Bryophytes'), angio = any(c => ['Eudicots', 'Monocots', 'Basal_angiosperm'].includes(c));
  if (land && !chlamy) return { pattern: 'Chlamydomonas_loss', confidence: 0.9 };
  if (land && chlamy) return { pattern: 'Conserved', confidence: 0.8 };
  if (angio && !bryo) return { pattern: 'Angiosperm_specific', confidence: 0.7 };
  if (bryo && angio) return { pattern: 'Land_plant_innovation', confidence: 0.8 };
  if (bryo && !angio) return { pattern: 'Bryophyte_origin', confidence: 0.6 };
  if (present.length <= 1) return { pattern: 'Sporadic', confidence: 0.5 };
  return { pattern: 'Unknown', confidence: 0 };
}

// --- docking.py ------------------------------------------------------------------------------------------------------
export function dockingSpecificity(insp3, insp4) {
  const delta = insp4 - insp3;
  return { insp3_affinity: insp3, insp4_affinity: insp4, delta_affinity: delta, is_specific: delta > 2,
    specificity_ratio: Math.abs(insp3 / (insp4 + 0.1)) };
}

// Same composition as POST /api/score in web/app.py: raw evidence -> CandidateScore fields, total and priority.
export const SCORE_DEFAULTS = { affinity: -8, mean_plddt: 80, tm_region_plddt: null, pocket_raw: 8, pattern: 'Unknown',
  copies: 1, ancestral: 1, correlation: 0, tissues: [], expected_tissues: [], localization: 'cytosol',
  has_signal: false, signal_type: null, n_tm: 0, pore_loop: false, selectivity_filter: false, topology_match: 0 };
export function computeScore(input) {
  const b = { ...SCORE_DEFAULTS, ...input };
  const fields = {
    pocket_score: scorePocket(b.pocket_raw), docking_score: scoreDocking(b.affinity),
    structure_quality: scoreStructureQuality(b.mean_plddt, b.tm_region_plddt),
    phylo_pattern_score: scorePhyloPattern(b.pattern), expansion_score: scoreGeneExpansion(b.copies, b.ancestral),
    coexpression_score: scoreCoexpression(b.correlation),
    tissue_specificity: scoreTissueSpecificity(b.tissues, b.expected_tissues),
    localization_score: scoreLocalization(b.localization), signal_peptide_score: scoreSignalPeptide(b.has_signal, b.signal_type),
    tm_domain_score: scoreTmDomains(b.n_tm), pore_architecture: scorePoreArchitecture(b.pore_loop, b.selectivity_filter),
    topology_match: Math.max(0, Math.min(5, b.topology_match)),
  };
  const total = totalScore(fields);
  return { engine: 'javascript', fields, total, priority: priorityClass(total) };
}
