// Audit page: what the committed files are, and benchmarks of the shipped rules computed by the repo's Python code
// (GET /audit, exported as JSON for the static site). It names what each instrument cannot detect.
import { route, app, api, esc, list } from './core.js';

const KIND = { 'example input': 'info', placeholder: 'no', 'figure, untraceable': 'warn' };

route('audit', async () => {
  const [files, a] = await Promise.all([api('/files'), api('/audit')]);
  const ef = a.ef_hand, tm = a.tm_min_length, tb = a.tm_benchmark, ph = a.phylo;
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
      reproduce calmodulin's EF-hand region. Decoys: the same loops with the position-1 Asp (a Ca²⁺ ligand) changed to
      Ala. The motif is ${esc(ef.source)}, <code class="wrapcode">${esc(ef.pattern)}</code>: 13 positions, the loop plus
      the hydrophobic residue after it.</p>
    <div class="scroll"><table><thead><tr><th>Loop (+ next)</th><th class="num">Position</th><th>Loop</th><th>D1A decoy</th></tr></thead><tbody>
      ${list(ef.loops, l => `<tr><td class="mono">${esc(l.loop)}<span class="mute">${esc(l.next)}</span></td><td class="num">${l.position}</td>
        <td>${l.hit ? '<span class="badge ok">match</span>' : '<span class="badge no">no match</span>'}</td>
        <td>${l.d1a_mutant_hit ? '<span class="badge no">match</span>' : '<span class="badge ok">rejected</span>'}</td></tr>`)}</tbody></table></div>
    <p class="callout"><b>${ef.detected} of ${ef.loops.length}</b> loops found and ${ef.decoys_detected} of
      ${ef.loops.length} decoys matched. What the pattern cannot do: it reads sequence only, so it cannot tell a
      Ca²⁺-binding EF hand from a degenerate one that has lost binding, and it misses non-canonical loops (for
      example the 14-residue S100 pseudo-EF-hand). The earlier pattern, <code>D.{2}[DNS]...[DNS].{2}[EQ]</code>,
      found none of these loops.</p>

    <h2>2 · TM prediction against UniProt annotations</h2>
    <p>The standard Kyte &amp; Doolittle rule (19-residue window, mean hydropathy ≥ 1.6; J. Mol. Biol. 1982,
      157:105-132) on UniProt-annotated membrane proteins and two soluble decoys (${esc(tb.source.split(',')[0])},
      retrieved ${esc(tb.retrieved)}). An annotated helix counts as found if a predicted segment overlaps it.</p>
    <div class="scroll"><table><thead><tr><th>Protein</th><th class="num">Annotated TM</th><th class="num">Found</th><th class="num">Extra</th><th>Predicted segments</th></tr></thead><tbody>
      ${list(tb.rows, r => `<tr><td><span class="mono">${esc(r.accession)}</span> ${esc(r.name)}</td><td class="num">${r.annotated.length}</td>
        <td class="num">${r.found}</td><td class="num">${r.extra}</td><td class="small mono">${esc(r.predicted.map(p => p.join('-')).join(', ') || '-')}</td></tr>`)}</tbody></table></div>
    <p class="callout"><b>${tb.found} of ${tb.annotated}</b> annotated helices found, ${tb.extra} extra segment(s), no
      segment in either soluble decoy. What the rule cannot do: it misses less hydrophobic helices (2 of 7 in
      bacteriorhodopsin, 1 of 6 in aquaporin-1), it calls hydrophobic signal peptides TM segments (the extra
      segment at the glycophorin A N-terminus overlaps its UniProt signal peptide, residues 1-19).</p>
    <p class="small mute">Shortest poly-Leu stretch called (<code>${esc(tm.construct)}</code>): <b>${tm.shortest_detected}</b>
      residues. The earlier default needed 15 consecutive windows above threshold (23 residues of pure Leu) and found
      none of these helices except one in EnvZ; it is still available as <code>min_windows=15</code>.</p>

    <h2>3 · The example sequences under the default filter</h2>
    <div class="scroll"><table><thead><tr><th>ID</th><th>Description in the FASTA</th><th class="num">Length</th><th class="num">TM segments</th></tr></thead><tbody>
      ${list(a.examples_tm, e => `<tr><td class="mono">${esc(e.id)}</td><td>${esc(e.description)}</td><td class="num">${e.length}</td>
        <td class="num">${e.n_tm}${e.domains.length ? ` <span class="mute small">(${e.domains.map(d => d.join('-')).join(', ')})</span>` : ''}</td></tr>`)}</tbody></table></div>
    <p class="small mute">With the default filter (2-10 segments, 200-2000 aa):
      ${list(a.examples_tm, e => `${esc(e.id)} ${e.n_tm >= 2 && e.n_tm <= 10 && e.length >= 200 && e.length <= 2000 ? 'passes' : 'fails'}`, ', ')}.
      These are designed sequences, so their counts say nothing about real proteins.</p>

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
