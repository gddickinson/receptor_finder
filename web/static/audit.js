// Audit page: what the committed files are, and benchmarks of the shipped rules computed by the repo's Python code
// (GET /audit, exported as JSON for the static site). It names what each instrument cannot detect.
import { route, app, api, esc, list } from './core.js';

const KIND = { 'example input': 'info', placeholder: 'no', 'figure, untraceable': 'warn' };

route('audit', async () => {
  const [files, a] = await Promise.all([api('/files'), api('/audit')]);
  const ef = a.ef_hand, tm = a.tm_min_length, ph = a.phylo;
  app().innerHTML = `<p class="kicker">Evidence audit</p><h1>What is committed, and what the rules can detect</h1>
    <p class="lead">Every number on this page is computed by the repo's own Python code (web/audit.py) from the committed
      files. The repository holds no proteome search, no ranked candidate list and no real structure, pocket or
      docking result, so this site shows none.</p>

    <h2>Committed data and outputs</h2>
    <div class="scroll"><table><thead><tr><th>File</th><th>Status</th><th class="num">Bytes</th><th>What it is</th></tr></thead><tbody>
      ${list(files, f => `<tr><td class="mono small">${esc(f.path)}<br><span class="mute">sha256 ${esc(f.sha256)}…</span></td>
        <td><span class="badge ${KIND[f.kind] || 'info'}">${esc(f.kind)}</span></td><td class="num">${f.bytes.toLocaleString('en')}</td>
        <td class="small">${esc(f.note)}</td></tr>`)}</tbody></table></div>

    <h2>1 · EF-hand motif against calmodulin's loops</h2>
    <p>Positive set: the four 12-residue Ca²⁺-binding loops in example sequence AT2G23456, whose first ~130 residues
      reproduce calmodulin's EF-hand region. The shipped pattern <code>${esc(ef.pattern)}</code> places its conserved
      positions at 1, 4, 8 and 11; the canonical loop has its Ca²⁺ ligands at 1, 3, 5, 7, 9 and 12.</p>
    <div class="scroll"><table><thead><tr><th>Loop</th><th class="num">Position in AT2G23456</th><th>Shipped EF_hand pattern</th></tr></thead><tbody>
      ${list(ef.loops, l => `<tr><td class="mono">${esc(l.loop)}</td><td class="num">${l.position}</td>
        <td>${l.hits.length ? '<span class="badge ok">match</span>' : '<span class="badge no">no match</span>'}</td></tr>`)}</tbody></table></div>
    <p class="callout"><b>${ef.detected} of ${ef.loops.length}</b> calmodulin EF-hand loops detected, and
      ${ef.hits_in_AT2G23456.length} EF_hand hits in the whole of AT2G23456 (described in the FASTA as a “protein with
      Ca²⁺ binding motifs”). The pattern cannot be used to find EF hands as written; a proposed correction is listed
      in the project report for George to decide.</p>

    <h2>2 · Shortest hydrophobic stretch called a TM segment</h2>
    <p>Construct: <code>${esc(tm.construct)}</code> (poly-Leu, the most hydrophobic common helix, between Lys flanks), run
      through predict_tm_domains with its defaults for n = ${tm.rows[0].n}-${tm.rows.at(-1).n}.</p>
    <p class="callout"><b>Shortest detected: ${tm.shortest_detected ?? 'none'} residues.</b> A segment needs 15
      consecutive above-threshold windows, so even a pure-leucine 20-residue helix (a typical TM length) is not
      called. Natural TM helices, which mix in less hydrophobic residues, need longer runs still. The default filter
      (≥ 2 segments) will therefore miss many genuine multi-pass membrane proteins.</p>
    <div class="scroll"><table><thead><tr><th class="num">n (Leu)</th>${list(tm.rows.filter(r => r.n % 2 === 0 || r.n === tm.shortest_detected), r => `<th class="num">${r.n}</th>`)}</tr></thead>
      <tbody><tr><td>called</td>${list(tm.rows.filter(r => r.n % 2 === 0 || r.n === tm.shortest_detected), r => `<td class="num">${r.domains.length ? '✓' : '·'}</td>`)}</tr></tbody></table></div>

    <h2>3 · The example sequences under the default filter</h2>
    <div class="scroll"><table><thead><tr><th>ID</th><th>Description in the FASTA</th><th class="num">Length</th><th class="num">TM segments</th></tr></thead><tbody>
      ${list(a.examples_tm, e => `<tr><td class="mono">${esc(e.id)}</td><td>${esc(e.description)}</td><td class="num">${e.length}</td>
        <td class="num">${e.n_tm}${e.domains.length ? ` <span class="mute small">(${e.domains.map(d => d.join('-')).join(', ')})</span>` : ''}</td></tr>`)}</tbody></table></div>
    <p class="small mute">None of the three reaches the default minimum of 2 TM segments, including the one described as
      “multi-pass”, so the GUI demo's filter step passes nothing with default settings.</p>

    <h2>4 · Phylogenetic patterns the classifier can return</h2>
    <p>classify_distribution run on all ${ph.combinations} presence/absence combinations of its 8 species:</p>
    <p>${list(Object.entries(ph.counts), ([k, v]) => `<span class="tag a">${esc(k)}: ${v}</span>`)}</p>
    <p class="callout">Never returned: ${list(ph.never_returned, k => `<span class="tag b">${esc(k)}</span>`)}. Any land
      plant presence is caught by the first two branches (Chlamydomonas_loss / Conserved), so the angiosperm- and
      bryophyte-specific branches below them cannot be reached, and a family absent from every species is called
      “Sporadic”. The right branch order is a scientific choice, left for George.</p>

    <h2>5 · Ligand library (docking.py)</h2>
    <p>Phosphorus atoms counted in each SMILES string against the phosphate groups the named ligand has.</p>
    <div class="scroll"><table><thead><tr><th>Ligand</th><th class="num">P in SMILES</th><th class="num">Expected</th><th>Notes</th></tr></thead><tbody>
      ${list(a.ligands, l => `<tr><td class="mono">${esc(l.name)}</td><td class="num">${l.phosphorus}</td><td class="num">${l.expected_phosphates}</td>
        <td class="small">${l.notes.length ? list(l.notes, x => `<span class="err">${esc(x)}</span>`, '<br>') : '<span class="ok">no problem found by these checks</span>'}</td></tr>`)}</tbody></table></div>
    <p class="small mute">Docking itself is a placeholder in this repo, so these SMILES have not yet affected any
      result; they would once real docking is wired in.</p>`;
});
