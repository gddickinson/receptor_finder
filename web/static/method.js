// Method page: the five phases from the README, with what the code actually does at each step (from src/).
import { route, app, api, esc, list } from './core.js';
import { FIELDS } from './score_engine.js';

const STATUS = { real: ['ok', 'implemented'], rule: ['info', 'rules only'], mock: ['no', 'placeholder'], none: ['warn', 'not in code'] };
const badge = s => `<span class="badge ${STATUS[s][0]}">${STATUS[s][1]}</span>`;

// Steps and tools as written in README.md ("Workflow Overview"); the status of each comes from reading src/.
const PHASES = [
  { n: 1, name: 'Sequence analysis and filtering', file: 'sequence_analysis.py', input: 'Plant proteome (FASTA)', steps: [
    ['Filter by length (200-2000 aa)', 'real', 'SequenceFilter.filter_by_length'],
    ['Predict transmembrane domains (keep 2-10)', 'real', 'Kyte-Doolittle sliding window, TransmembraneDomainPredictor'],
    ['Scan for binding motifs (Arg/Lys clusters, aromatics)', 'real', 'MotifScanner: 6 regular expressions; not used by the default filter'],
    ['Predict subcellular localization (tonoplast/ER preferred)', 'none', 'no predictor in the code; the score takes a localization you supply'],
    ['BLAST / HMMER homology search', 'none', 'listed as optional external tools; no wrapper in the code'],
  ] },
  { n: 2, name: 'Structure prediction', file: 'structure_prediction.py', input: 'Filtered sequences', steps: [
    ['AlphaFold2 / ColabFold prediction', 'mock', 'predict_structure writes a straight CA trace with simulated pLDDT; get_colabfold_command prints the real command'],
    ['Assess quality (mean pLDDT, TM-region pLDDT)', 'real', 'StructureAnalyzer.parse_pdb_plddt / assess_structure_quality, on any real PDB'],
    ['TM regions from the structure', 'mock', 'identify_tm_regions_from_structure returns a 20-residue region every 100 residues'],
  ] },
  { n: 3, name: 'Binding pocket analysis', file: 'binding_pocket.py', input: 'Structures', steps: [
    ['Detect pockets (fpocket / P2Rank)', 'mock', 'detect_pockets returns the same three invented pockets for any file; get_fpocket_command prints the real command'],
    ['Score pockets for InsP₃ (0-13) and cADPR (0-12)', 'rule', 'PocketScorer: volume, charge, residue-class counts, depth; runs on the Scoring page'],
    ['Specificity (InsP₃ over InsP₄)', 'rule', 'PocketCharacterizer.assess_specificity, needs an InsP₄ pocket score'],
  ] },
  { n: 4, name: 'Molecular docking', file: 'docking.py', input: 'Top pockets', steps: [
    ['Dock InsP₃, InsP₄, cADPR with AutoDock Vina', 'mock', 'dock_ligand returns random affinities; get_vina_command prints the real command'],
    ['Ligand library (SMILES)', 'rule', 'six SMILES in LigandLibrary; the Audit page lists problems found in them'],
    ['Specificity: InsP₃ at least 2 kcal/mol better than InsP₄', 'rule', 'DockingAnalyzer.assess_specificity'],
  ] },
  { n: 5, name: 'Multi-criteria ranking', file: 'scoring.py', input: 'All evidence', steps: [
    ['100-point score in five categories, HIGH ≥ 60, MEDIUM 40-59', 'real', 'CandidateScore, CandidateRanker; calculator on the Scoring page'],
    ['Expression evidence (co-expression, tissues)', 'rule', 'scoring rules only; no expression data source in the code'],
    ['Phylogenetic pattern', 'rule', 'phylogenetics.classify_distribution over 8 species; see the Audit page'],
  ] },
];

route('method', async () => {
  const info = await api('/info');
  const cats = [...new Set(FIELDS.map(f => f[2]))];
  app().innerHTML = `<p class="kicker">Method</p><h1>Five phases, from proteome to ranked candidates</h1>
    <p class="lead">The pipeline as described in the README, with what the code does at each step today.
      “Placeholder” steps return mock values and need the external tool (ColabFold, fpocket, AutoDock Vina) to give
      real results; the README's expected funnel (~500 → ~300 → ~100 → ~50 → 20-50 candidates for an Arabidopsis
      proteome) is a projection, not a computed result.</p>
    <div class="pipe">${list(PHASES, p => `<div class="panel"><h3>Phase ${p.n} · ${esc(p.name)}</h3>
      <p class="small mute">src/${esc(p.file)} · input: ${esc(p.input)}</p>
      <div class="scroll"><table><tbody>${list(p.steps, ([s, st, how]) => `<tr><td>${esc(s)}</td><td>${badge(st)}</td>
        <td class="small mute">${esc(how)}</td></tr>`)}</tbody></table></div></div>`)}</div>
    <h2>The 100-point rubric</h2>
    <p>From <code>CandidateScore</code> in src/scoring.py. Priority: HIGH ≥ 60, MEDIUM 40-59, LOW &lt; 40.</p>
    <div class="scroll"><table><thead><tr><th>Category</th><th>Criterion</th><th class="num">Max</th></tr></thead><tbody>
      ${list(cats, c => list(FIELDS.filter(f => f[2] === c), (f, i) => `<tr><td>${i ? '' : esc(c)}</td><td class="mono">${esc(f[0])}</td><td class="num">${f[1]}</td></tr>`))}
      <tr><td></td><td><b>Total</b></td><td class="num"><b>${FIELDS.reduce((s, f) => s + f[1], 0)}</b></td></tr></tbody></table></div>
    <h2>Default filter</h2>
    <p>Length ${info.filter_defaults.min_length}-${info.filter_defaults.max_length} aa, then
      ${info.filter_defaults.min_tm}-${info.filter_defaults.max_tm} predicted TM segments (window ${info.tm_defaults.window},
      threshold ${info.tm_defaults.threshold}; a segment needs 15 consecutive windows above threshold, and segments
      within 5 residues are merged). <a href="#/sequence">Try it</a>.</p>
    <h2>Validation the README proposes</h2>
    <p class="mute">CRISPR knockouts, patch-clamp of isolated vacuoles, Ca²⁺ imaging of InsP₃/cADPR-induced release,
      [³H]InsP₃ or [³²P]cADPR binding assays and complementation. None of this is in the repo.</p>`;
});
