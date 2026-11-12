// Audit page: what the committed files are, and benchmarks of the shipped rules computed by the repo's Python code
// (GET /audit, exported as JSON for the static site). It names what each instrument cannot detect.
import { route, app, api, esc, list } from './core.js';

const KIND = { 'example input': 'info', placeholder: 'no', 'figure, untraceable': 'warn' };

route('audit', async () => {
  const [files, a] = await Promise.all([api('/files'), api('/audit')]);
  const ef = a.ef_hand, eb = a.ef_hand_benchmark, cm = a.ef_hand.calmodulin_match, tm = a.tm_min_length, tb = a.tm_benchmark, ph = a.phylo;
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
    <p>Positive set: the four 12-residue Ca²⁺-binding loops in example sequence AT2G23456, whose residues
      ${cm.example_range.join('-')} are identical to human calmodulin (UniProt ${esc(cm.accession)}, residues
      ${cm.uniprot_range.join('-')}). Decoys: the same loops with the position-1 Asp (a Ca²⁺ ligand) changed to Ala.
      The motif is the ${esc(ef.source)}, <code class="wrapcode">${esc(ef.pattern)}</code>: Asp at 1, side-chain oxygen
      ligands at 3 and 5, Gly at 6, Glu or Asp at 12.</p>
    <div class="scroll"><table><thead><tr><th>Loop</th><th class="num">Position</th><th>Loop</th><th>D1A decoy</th></tr></thead><tbody>
      ${list(ef.loops, l => `<tr><td class="mono">${esc(l.loop)}</td><td class="num">${l.position}</td>
        <td>${l.hit ? '<span class="badge ok">match</span>' : '<span class="badge no">no match</span>'}</td>
        <td>${l.d1a_mutant_hit ? '<span class="badge no">match</span>' : '<span class="badge ok">rejected</span>'}</td></tr>`)}</tbody></table></div>
    <p>On UniProt-annotated EF-hand proteins (${esc(eb.source.split(',')[0])}, retrieved ${esc(eb.retrieved)}), a loop
      counts as annotated when UniProt lists Ca²⁺ ligands at its positions 1, 3, 5 and 12:</p>
    <div class="scroll"><table><thead><tr><th>Protein</th><th class="num">EF-hand domains</th><th>Annotated loops</th><th>Motif hits</th></tr></thead><tbody>
      ${list(eb.rows, r => `<tr><td><span class="mono">${esc(r.accession)}</span> ${esc(r.name)}</td><td class="num">${r.ef_hand_domains}</td>
        <td class="small mono">${esc(r.annotated_loops.join(', ') || '-')}</td><td class="small mono">${esc(r.hits.join(', ') || '-')}</td></tr>`)}</tbody></table></div>
    <p class="callout"><b>${ef.detected} of ${ef.loops.length}</b> example loops found and ${ef.decoys_detected} of
      ${ef.loops.length} decoys matched; <b>${eb.found} of ${eb.annotated}</b> UniProt-annotated loops found with
      ${eb.extra} other hit(s), and no hit in the ${eb.ef_hand_domains - eb.annotated} EF-hand domains that do not bind
      Ca²⁺ this way (S100 pseudo-EF-hands, troponin C site I, recoverin EF1 and EF4). What the motif cannot do: it
      reads sequence only, so a loop that fits it may still have lost Ca²⁺ binding, and it misses non-canonical loops
      (the 14-residue S100 pseudo-EF-hand). The original pattern, <code>D.{2}[DNS]...[DNS].{2}[EQ]</code>, found none
      of the calmodulin loops.</p>

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
    <p class="callout">${ph.never_returned.length ? `Never returned: ${list(ph.never_returned, k => `<span class="tag b">${esc(k)}</span>`)}.`
      : `All ${Object.keys(ph.counts).length} documented patterns are reachable.`} The patterns follow the definitions in
      docs/in_silico_calcium_projects.md (section 3.1). What the classifier cannot do: presence/absence across 8 species
      cannot separate a loss from a gain or detect re-expansion (copy number); those need a gene tree. An earlier
      version could return only Chlamydomonas_loss, Conserved and Sporadic, and gave “Chlamydomonas_loss” to families
      absent from the alga, the opposite of the documented definition.</p>

    <h2>5 · Ligand library (docking.py)</h2>
    <p>SMILES from PubChem (isomeric SMILES via the PUG REST API, retrieved 2026-11-12; see THIRD_PARTY.md). Checked here:
      phosphorus atoms against the phosphate groups each ligand has, a carbocyclic inositol ring, and no cytosine base.</p>
    <div class="scroll"><table><thead><tr><th>Ligand</th><th>Source</th><th class="num">P in SMILES</th><th class="num">Expected</th><th>Checks</th></tr></thead><tbody>
      ${list(a.ligands, l => `<tr><td class="mono">${esc(l.name)}</td><td class="small">${esc(l.source)}</td><td class="num">${l.phosphorus}</td><td class="num">${l.expected_phosphates}</td>
        <td class="small">${l.notes.length ? list(l.notes, x => `<span class="err">${esc(x)}</span>`, '<br>') : '<span class="ok">pass</span>'}</td></tr>`)}</tbody></table></div>
    <p class="small mute">8-NH2-cADPR is the ChEMBL record CHEMBL3350988 (PubChem CID 118719244), written in its neutral
      form: cADPR (CID 123847) with an amino group on adenine C8, C15H22N6O13P2 (Walseth &amp; Lee 1993). PubChem's
      record found by the name 8-amino-cADPR (CID 135600320) is not cyclic and is not used. The earlier library had InsP3/InsP4 as pyranose sugars with 4/5 phosphates,
      cytosine in cADPR and 8-NH2-cADPR, and a nitro group in NAADP. Docking itself is still a placeholder.</p>`;
});
