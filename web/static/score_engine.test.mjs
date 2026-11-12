// node web/static/score_engine.test.mjs: the JS scoring engine must reproduce the Python scorers on every case in
// tests/fixtures/engine_fixture.json (regenerate with scripts/make_engine_fixture.py).
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import * as S from './score_engine.js';

const fx = JSON.parse(readFileSync(new URL('../../tests/fixtures/engine_fixture.json', import.meta.url)));
const close = (a, b, msg) => assert.ok(Math.abs(a - b) <= 1e-9 * Math.max(1, Math.abs(b)), `${msg}: ${a} vs ${b}`);
const sc = fx.scores;
let n = 0;
const check = (rows, fn, name) => rows.forEach(r => { close(fn(...r.slice(0, -1)), r.at(-1), `${name} ${JSON.stringify(r)}`); n++; });
check(sc.docking, S.scoreDocking, 'docking');
check(sc.quality, S.scoreStructureQuality, 'quality');
check(sc.pocket, S.scorePocket, 'pocket');
check(sc.pattern, S.scorePhyloPattern, 'pattern');
check(sc.expansion, S.scoreGeneExpansion, 'expansion');
check(sc.coexpression, S.scoreCoexpression, 'coexpression');
check(sc.tissue, S.scoreTissueSpecificity, 'tissue');
check(sc.localization, S.scoreLocalization, 'localization');
check(sc.signal, S.scoreSignalPeptide, 'signal');
check(sc.tm, S.scoreTmDomains, 'tm');
check(sc.pore, S.scorePoreArchitecture, 'pore');

for (const p of fx.pockets) {
  const c = S.residueCounts(p.residues);
  assert.deepEqual([c.basic, c.aromatic, c.acidic, c.polar], p.counts);
  close(S.scoreInsp3(p), p.insp3, 'insp3'); close(S.scoreCadpr(p), p.cadpr, 'cadpr');
  const a = S.pocketSpecificity(p.insp3, p.cadpr, p.insp4, p.nad), b = p.specificity;
  assert.equal(a.primary_ligand, b.primary_ligand); assert.equal(a.is_specific, b.is_specific);
  if (b.specificity_ratio === null) assert.equal(a.specificity_ratio, null); else close(a.specificity_ratio, b.specificity_ratio, 'ratio');
  n++;
}
for (const c of fx.phylo) { assert.deepEqual(S.classifyDistribution(c.presence), { pattern: c.pattern, confidence: c.confidence }); n++; }
for (const d of fx.docking) {
  const a = S.dockingSpecificity(d.insp3_affinity, d.insp4_affinity);
  close(a.delta_affinity, d.delta_affinity, 'delta'); close(a.specificity_ratio, d.specificity_ratio, 'sr');
  assert.equal(a.is_specific, d.is_specific); n++;
}
for (const c of fx.candidates) {
  const a = S.computeScore(c.input), b = c.output;
  for (const k of Object.keys(b.fields)) close(a.fields[k], b.fields[k], `candidate ${k}`);
  close(a.total, b.total, 'total'); assert.equal(a.priority, b.priority); n++;
}
assert.equal(S.priorityClass(60), 'HIGH'); assert.equal(S.priorityClass(59.9), 'MEDIUM'); assert.equal(S.priorityClass(39.9), 'LOW');
assert.equal(S.FIELDS.reduce((s, f) => s + f[1], 0), 100);
console.log(`score_engine: ${n} cases match Python (11 scorer functions, pockets, 256 phylogenetic patterns, docking specificity, whole candidates)`);
