// Sequence filter page: paste FASTA, run TM prediction + motif scan + the default filter (Python locally, JS on the
// static site), then inspect one sequence's hydropathy profile and residue map.
import { route, app, api, esc, list, $, isStatic, download } from './core.js';
import { parseFasta, MOTIFS } from './seq_engine.js';
import { profileSvg, residueMap } from './seq_view.js';

let state = { fasta: '', results: null, engine: '', pick: 0, motif: 'Arg_cluster' };
const yes = (ok, a = 'pass', b = 'fail') => `<span class="badge ${ok ? 'ok' : 'no'}">${ok ? a : b}</span>`;

route('sequence', async () => {
  const [info, ex] = await Promise.all([api('/info'), api('/examples')]);
  const f = info.filter_defaults, t = info.tm_defaults;
  if (!state.fasta) state.fasta = ex.fasta;
  app().innerHTML = `<p class="kicker">Phase 1 · sequence analysis</p><h1>Sequence filter</h1>
    <p class="lead">Kyte-Doolittle TM prediction, the six motif patterns and the length / TM-count filter from
      src/sequence_analysis.py, ${isStatic() ? 'ported to JavaScript and run in your browser' : 'run by the Python code'}.
      Your sequences stay ${isStatic() ? 'in this browser tab' : 'on this computer'}.</p>
    <div class="side"><div class="panel">
      <label for="fa">FASTA (pre-filled with the ${ex.records.length} hypothetical examples in ${esc(ex.source)})</label>
      <textarea id="fa" rows="9" spellcheck="false">${esc(state.fasta)}</textarea>
      <div class="row">
        <div><label for="w">Window</label><input id="w" type="number" min="5" max="60" value="${t.window}"></div>
        <div><label for="th">Threshold</label><input id="th" type="number" step="0.1" value="${t.threshold}"></div></div>
      <div class="row">
        <div><label for="mw">Min windows</label><input id="mw" type="number" min="1" value="${t.min_windows}"></div>
        <div><label for="mg">Merge gap</label><input id="mg" type="number" min="0" value="${t.merge_gap}"></div></div>
      <div class="row">
        <div><label for="minl">Min length</label><input id="minl" type="number" value="${f.min_length}"></div>
        <div><label for="maxl">Max length</label><input id="maxl" type="number" value="${f.max_length}"></div></div>
      <div class="row">
        <div><label for="mint">Min TM</label><input id="mint" type="number" value="${f.min_tm}"></div>
        <div><label for="maxt">Max TM</label><input id="maxt" type="number" value="${f.max_tm}"></div></div>
      <button id="run">Analyse</button><button id="ex" class="alt">Reset to examples</button>
      <p class="small mute">Defaults are the Python defaults: the Kyte &amp; Doolittle (1982) rule, a 19-residue window
        with mean hydropathy ≥ 1.6. The example IDs are placeholders, not Arabidopsis loci.</p>
    </div><div id="out"></div></div>`;
  const run = async () => {
    state.fasta = $('#fa').value;
    const recs = parseFasta(state.fasta);
    if (!recs.length) { $('#out').innerHTML = '<p class="err">No FASTA records found (headers start with &gt;).</p>'; return; }
    const body = { fasta: state.fasta, window: +$('#w').value, threshold: +$('#th').value,
      min_windows: +$('#mw').value, merge_gap: +$('#mg').value, min_length: +$('#minl').value,
      max_length: +$('#maxl').value, min_tm: +$('#mint').value, max_tm: +$('#maxt').value };
    $('#run').disabled = true;
    try {
      const r = await api('/analyse', { body });
      state = { ...state, results: r.results, engine: r.engine, body, seqs: Object.fromEntries(recs.map(x => [x.id, x.sequence])) };
      state.pick = Math.min(state.pick, r.results.length - 1);
      show();
    } catch (e) { $('#out').innerHTML = `<p class="err">${esc(e.message)}</p>`; }
    $('#run').disabled = false;
  };
  $('#run').onclick = run;
  $('#ex').onclick = () => { $('#fa').value = ex.fasta; run(); };
  run();
});

function show() {
  const { results: rs, body } = state, mk = Object.keys(MOTIFS);
  const passN = rs.filter(r => r.passes).length;
  const pick = rs[state.pick], seq = state.seqs[pick.id] || '';
  $('#out').innerHTML = `<div class="panel"><h3>${rs.length} sequences · ${passN} pass the filter</h3>
    <p class="small mute">Engine: ${esc(state.engine)}. Length ${body.min_length}-${body.max_length}, TM ${body.min_tm}-${body.max_tm}.</p>
    <div class="scroll"><table><thead><tr><th>ID</th><th class="num">Length</th><th class="num">TM</th><th>Length</th><th>TM</th>
      ${list(mk, m => `<th class="num" title="${esc(MOTIFS[m])}">${esc(m.replace('_', ' '))}</th>`)}</tr></thead><tbody>
      ${list(rs, (r, i) => `<tr class="click${i === state.pick ? ' on' : ''}" data-i="${i}" tabindex="0"><td class="mono">${esc(r.id)}</td>
        <td class="num">${r.length}</td><td class="num">${r.n_tm}</td><td>${yes(r.length_ok)}</td><td>${yes(r.tm_ok)}</td>
        ${list(mk, m => `<td class="num">${r.motifs[m].length}</td>`)}</tr>`)}</tbody></table></div>
    <button class="alt" id="csv">Download table (CSV)</button></div>
    <div class="panel"><h3>${esc(pick.id)} ${pick.description ? '· ' + esc(pick.description) : ''}</h3>
      ${profileSvg(seq, pick.tm, body.window, body.threshold)}
      <p class="small">TM segments: ${pick.tm.length ? list(pick.tm, ([a, b]) => `<span class="tag a">${a}-${b}</span>`) : '<span class="mute">none</span>'}</p>
      <label for="mf">Underline matches of</label>
      <select id="mf">${list(mk, m => `<option ${m === state.motif ? 'selected' : ''}>${esc(m)}</option>`)}</select>
      <p class="small mute mono wrapcode">${esc(MOTIFS[state.motif])} · ${pick.motifs[state.motif].length} non-overlapping matches at
        ${esc(pick.motifs[state.motif].slice(0, 30).join(', ') || '-')}${pick.motifs[state.motif].length > 30 ? ' ...' : ''}</p>
      ${residueMap(seq, pick.tm, MOTIFS[state.motif])}
      ${state.motif === 'EF_hand' ? '<p class="small mute">EF_hand is the canonical 12-residue Ca²⁺-binding loop (Asp at 1, O-ligands at 3 and 5, Gly at 6, Glu/Asp at 12); benchmark on the <a href="#/audit">audit</a> page.</p>' : ''}
    </div>`;
  document.querySelectorAll('tr.click').forEach(tr => {
    const go = () => { state.pick = +tr.dataset.i; show(); };
    tr.onclick = go; tr.onkeydown = e => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); go(); } };
  });
  $('#mf').onchange = e => { state.motif = e.target.value; show(); };
  $('#csv').onclick = () => download('receptor_finder_filter.csv', toCsv(rs, mk), 'text/csv');
}

function toCsv(rs, mk) {
  const q = s => `"${String(s ?? '').replace(/"/g, '""')}"`;
  const head = ['id', 'length', 'n_tm', 'tm_segments', 'length_ok', 'tm_ok', 'passes', ...mk];
  return [head.join(','), ...rs.map(r => [q(r.id), r.length, r.n_tm, q(r.tm.map(t => t.join('-')).join(';')),
    r.length_ok, r.tm_ok, r.passes, ...mk.map(m => r.motifs[m].length)].join(','))].join('\n') + '\n';
}
