// Static-site data layer: answers api() calls without the Python server.
// GET paths come from the JSON exported by scripts/build_static.py; POST /analyse and /score run the JS ports
// (seq_engine.js, score_engine.js), which are node-tested against the Python code on tests/fixtures/engine_fixture.json.
import { analyse, parseFasta } from './seq_engine.js';
import { computeScore } from './score_engine.js';

const files = new Map();
function file(rel) {
  if (!files.has(rel)) files.set(rel, fetch('api/' + rel + '.json').then(r => { if (!r.ok) throw new Error(`${r.status} ${rel}`); return r.json(); }));
  return files.get(rel);
}

export async function staticApi(path, body) {
  const clean = path.split('?')[0].replace(/^\/+/, '');
  if (clean === 'analyse') {
    const o = body || {};
    const results = analyse(parseFasta(o.fasta || ''), { window: o.window, threshold: o.threshold,
      minWindows: o.min_windows, mergeGap: o.merge_gap, minLength: o.min_length,
      maxLength: o.max_length, minTm: o.min_tm, maxTm: o.max_tm });
    return { engine: 'javascript', results };
  }
  if (clean === 'score') return computeScore(body || {});
  return file(clean);
}
