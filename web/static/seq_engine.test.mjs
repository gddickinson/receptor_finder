// node web/static/seq_engine.test.mjs: the JS sequence engine must reproduce src/sequence_analysis.py on every case in
// tests/fixtures/engine_fixture.json (regenerate with scripts/make_engine_fixture.py).
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { hydroProfile, predictTmDomains, scanAll, parseFasta, analyse } from './seq_engine.js';

const fx = JSON.parse(readFileSync(new URL('../../tests/fixtures/engine_fixture.json', import.meta.url)));
let n = 0;
for (const c of fx.sequences) {
  const prof = hydroProfile(c.sequence);
  assert.equal(prof.length, c.profile.length, c.id);
  prof.forEach((v, i) => assert.ok(Math.abs(v - c.profile[i]) <= 1e-12, `${c.id} profile ${i}`));  // bitwise on macOS
  assert.deepEqual(predictTmDomains(c.sequence), c.tm, `${c.id} tm`);
  assert.deepEqual(predictTmDomains(c.sequence, 15, 1.0), c.tm_w15_t1, `${c.id} tm w15`);
  assert.deepEqual(predictTmDomains(c.sequence, 19, 1.6, 15, 0), c.tm_minw15_gap0, `${c.id} tm old rule`);
  assert.deepEqual(scanAll(c.sequence), c.motifs, `${c.id} motifs`);
  n++;
}
const fa = parseFasta('>A first one\nMKL\n  LLV \n>B\n\nAAA\n');
assert.deepEqual(fa, [{ id: 'A', description: 'first one', sequence: 'MKLLLV' }, { id: 'B', description: null, sequence: 'AAA' }]);
const ex = fx.sequences.filter(c => c.id.startsWith('AT')).map(c => ({ id: c.id, description: null, sequence: c.sequence }));
const res = analyse(ex);
assert.deepEqual(res.map(r => r.n_tm), fx.sequences.filter(c => c.id.startsWith('AT')).map(c => c.tm.length));
console.log(`seq_engine: ${n} sequences match Python (profile, TM domains at three settings, 6 motifs); FASTA ok`);
