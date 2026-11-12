// Scoring page: the 100-point candidate score (POST /score: Python locally, score_engine.js on the static site), plus
// the pocket rules, the phylogenetic pattern classifier and the docking specificity rule (score_engine.js, all tested
// against the Python). Inputs are evidence you supply: nothing here is a result for a real protein.
import { route, app, api, esc, list, $, $$ } from './core.js';
import * as S from './score_engine.js';

const num = (id, label, v, a = {}) => `<div><label for="${id}">${label}</label><input id="${id}" type="number" value="${v}"
  ${a.step ? `step="${a.step}"` : ''} ${a.min !== undefined ? `min="${a.min}"` : ''} ${a.max !== undefined ? `max="${a.max}"` : ''}></div>`;
const sel = (id, label, opts, v) => `<div><label for="${id}">${label}</label><select id="${id}">${list(opts, o =>
  `<option value="${esc(o[0])}" ${o[0] === v ? 'selected' : ''}>${esc(o[1])}</option>`)}</select></div>`;
const chk = (id, label, v = false) => `<label class="check" for="${id}"><input id="${id}" type="checkbox" ${v ? 'checked' : ''}>${label}</label>`;
const fmt = v => (Math.round(v * 100) / 100).toString();

route('scoring', async () => {
  const tabs = await api('/scoring');
  const pats = Object.keys(tabs.pattern_scores), locs = Object.keys(tabs.localization_scores);
  const sig = [['', 'no type given'], ...Object.keys(tabs.signal_type_scores).map(k => [k, k])];
  app().innerHTML = `<p class="kicker">Phase 5 · multi-criteria ranking</p><h1>Score a candidate</h1>
    <p class="lead">Enter the evidence for one protein; each scorer in src/scoring.py turns it into points. Use it to see
      how the rubric weighs the evidence. The values filled in are neutral defaults, not data for any protein.</p>
    <div class="two"><div class="panel" id="inp">
      <h3>Structural</h3><div class="row">${num('aff', 'Docking affinity (kcal/mol)', -8, { step: 0.1 })}
        ${num('pl', 'Mean pLDDT', 80, { min: 0, max: 100 })}</div>
      <div class="row">${num('tmpl', 'TM-region pLDDT (blank = none)', '', { min: 0, max: 100 })}
        ${num('pk', 'Pocket score (InsP₃ 0-13 or cADPR 0-12)', 8, { step: 0.5, min: 0 })}</div>
      <h3>Evolutionary</h3><div class="row">${sel('pat', 'Phylogenetic pattern', pats.map(p => [p, p.replace(/_/g, ' ')]), 'Unknown')}</div>
      <div class="row">${num('cp', 'Copies in focal species', 1, { min: 0 })}${num('anc', 'Ancestral copies', 1, { min: 0 })}</div>
      <h3>Expression</h3><div class="row">${num('r', 'Co-expression r with Ca²⁺ genes', 0, { step: 0.05, min: -1, max: 1 })}</div>
      <div class="row"><div><label for="tis">Tissues expressed (comma list)</label><input id="tis" value=""></div>
        <div><label for="etis">Expected tissues (blank = no expectation)</label><input id="etis" value=""></div></div>
      <h3>Localization</h3><div class="row">${sel('loc', 'Predicted localization', locs.map(l => [l, l.replace(/_/g, ' ')]), 'cytosol')}
        ${sel('st', 'Signal type', sig, '')}</div>${chk('sp', 'Signal peptide detected')}
      <h3>Topology</h3><div class="row">${num('ntm', 'Predicted TM segments', 0, { min: 0 })}${num('tmm', 'Similarity to channels (0-5)', 0, { step: 0.5, min: 0, max: 5 })}</div>
      ${chk('pore', 'Pore loop / re-entrant helix')}${chk('sf', 'Acidic selectivity filter')}
    </div><div><div class="panel" id="res"></div></div></div>
    <h2>Pocket rules (Phase 3)</h2><div class="two"><div class="panel">
      <div class="row">${num('pv', 'Volume (Å³)', 250)}${num('pd', 'Depth (Å)', 8, { step: 0.5 })}</div>
      <div class="row">${num('pc', 'Net charge', 2, { step: 0.5 })}${num('ph', 'Hydrophobicity (0-1)', 0.3, { step: 0.05 })}</div>
      <div><label for="pr">Pocket residues (three-letter codes)</label><input id="pr" value="ARG LYS ASP PHE TRP"></div>
      <div class="row">${num('pi4', 'InsP₄ pocket score (blank = untested)', '', { step: 0.5 })}</div>
      <p class="small mute">The defaults are the first pocket that the placeholder PocketDetector invents for every structure.</p>
    </div><div class="panel" id="pres"></div></div>
    <h2>Phylogenetic pattern (phylogenetics.py)</h2><div class="two"><div class="panel">
      ${list(Object.entries(tabs.plant_species), ([sp, v]) => chk('sp-' + sp, `<i>${esc(sp.replace(/_/g, ' '))}</i> <span class="mute small">${esc(v.clade)}</span>`, sp !== 'Chlamydomonas_reinhardtii'))}
    </div><div class="panel" id="phres"></div></div>
    <h2>Docking specificity (docking.py)</h2><div class="two"><div class="panel"><div class="row">
      ${num('d3', 'InsP₃ affinity (kcal/mol)', -9, { step: 0.1 })}${num('d4', 'InsP₄ affinity (kcal/mol)', -6.5, { step: 0.1 })}</div>
      <p class="small mute">Docking in this repo is a placeholder; enter affinities from a real Vina run.</p></div>
      <div class="panel" id="dres"></div></div>`;
  let t = 0;
  const upd = () => { clearTimeout(t); t = setTimeout(update, 120); };
  $$('#app input, #app select').forEach(el => el.addEventListener('input', upd));
  update();
});

const val = id => $('#' + id).value, n = id => +val(id), opt = id => (val(id) === '' ? null : +val(id));
const csv = id => val(id).split(',').map(s => s.trim()).filter(Boolean);

async function update() {
  const body = { affinity: n('aff'), mean_plddt: n('pl'), tm_region_plddt: opt('tmpl'), pocket_raw: n('pk'), pattern: val('pat'),
    copies: n('cp'), ancestral: n('anc'), correlation: n('r'), tissues: csv('tis'), expected_tissues: csv('etis'),
    localization: val('loc'), has_signal: $('#sp').checked, signal_type: val('st') || null, n_tm: n('ntm'),
    pore_loop: $('#pore').checked, selectivity_filter: $('#sf').checked, topology_match: n('tmm') };
  try {
    const r = await api('/score', { body });
    const cls = { HIGH: 'ok', MEDIUM: 'warn', LOW: 'no' }[r.priority];
    $('#res').innerHTML = `<h3>Total</h3><div class="total">${fmt(r.total)}<span class="mute small"> / 100</span></div>
      <p><span class="badge ${cls}">${r.priority}</span> <span class="small mute">engine: ${esc(r.engine)}</span></p>
      <table><tbody>${list(S.FIELDS, ([k, mx, cat]) => `<tr><td class="small">${esc(k.replace(/_/g, ' '))}<br><span class="mute">${cat}</span></td>
        <td style="width:40%"><div class="bar"><i style="width:${Math.max(0, Math.min(100, r.fields[k] / mx * 100))}%"></i></div></td>
        <td class="num">${fmt(r.fields[k])} / ${mx}</td></tr>`)}</tbody></table>`;
  } catch (e) { $('#res').innerHTML = `<p class="err">${esc(e.message)}</p>`; }

  const pocket = { volume: n('pv'), depth: n('pd'), charge: n('pc'), hydrophobicity: n('ph'),
    residues: val('pr').toUpperCase().split(/[\s,]+/).filter(Boolean) };
  const i3 = S.scoreInsp3(pocket), cd = S.scoreCadpr(pocket), c = S.residueCounts(pocket.residues);
  const spec = S.pocketSpecificity(i3, cd, opt('pi4'));
  $('#pres').innerHTML = `<h3>Pocket scores</h3><p><b class="mono">InsP₃ ${fmt(i3)} / 13</b> · <b class="mono">cADPR ${fmt(cd)} / 12</b></p>
    <p class="small mute">Residue classes: basic ${c.basic}, aromatic ${c.aromatic}, acidic ${c.acidic}, polar ${c.polar}.</p>
    <p>Primary ligand: <b>${spec.primary_ligand}</b> · specific: ${spec.is_specific ? 'yes' : 'no'}${spec.specificity_ratio !== null ? ` (ratio ${fmt(spec.specificity_ratio)})` : ''}</p>
    <p class="small mute">A cADPR call is marked specific without any test (assess_specificity assumes it); an InsP₃
      call needs an InsP₄ pocket score more than 2 points lower.</p>`;

  const presence = Object.fromEntries($$('#app input[id^="sp-"]').map(el => [el.id.slice(3), el.checked]));
  const ph = S.classifyDistribution(presence);
  $('#phres').innerHTML = `<h3>Pattern</h3><p><b class="mono">${esc(ph.pattern)}</b> · confidence ${ph.confidence} ·
    ${S.scorePhyloPattern(ph.pattern)} / 10 points</p><p class="small mute">Only three patterns are reachable with the
    current branch order: see the <a href="#/audit">audit</a>.</p>`;

  const d = S.dockingSpecificity(n('d3'), n('d4'));
  $('#dres').innerHTML = `<h3>Specificity</h3><p>ΔG(InsP₄) − ΔG(InsP₃) = <b class="mono">${fmt(d.delta_affinity)}</b> kcal/mol ·
    ${d.is_specific ? '<span class="badge ok">specific</span>' : '<span class="badge no">not specific</span>'}</p>
    <p class="small mute">Rule: InsP₃ must dock more than 2 kcal/mol better than InsP₄. Docking score for the InsP₃
      affinity: ${fmt(S.scoreDocking(n('d3')))} / 10.</p>`;
}
