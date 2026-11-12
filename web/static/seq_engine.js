// Pure JS port of src/sequence_analysis.py: Kyte-Doolittle TM prediction, regex motif scanning, FASTA parsing and the
// length / TM-count filters. No DOM. Node-tested against the Python code (seq_engine.test.mjs, tests/fixtures).

// Kyte-Doolittle hydropathy (TransmembraneDomainPredictor.HYDROPHOBICITY). Unknown letters count as 0, as in Python.
export const KD = { A: 1.8, R: -4.5, N: -3.5, D: -3.5, C: 2.5, Q: -3.5, E: -3.5, G: -0.4, H: -3.2, I: 4.5,
  L: 3.8, K: -3.9, M: 1.9, F: 2.8, P: -1.6, S: -0.8, T: -0.7, W: -0.9, Y: -1.3, V: 4.2 };

// MotifScanner.MOTIFS (Python regex syntax that JS reads the same way), scanned case-insensitively.
export const MOTIFS = {
  P_loop: 'G[KR].{1,5}G[KR][ST]',
  Walker_B: '[RK].{2,4}[LIVM].{2,4}D',
  EF_hand: 'D[^W][DNS][^ILVFYW][DENSTG][DNQGHRK][^GP][LIVMC][DENQSTAGC].{2}[DE][LIVMFYW]',  // PROSITE PS00018
  Acidic_cluster: '[DE]{3,5}',
  Arg_cluster: '[RK]{2,4}',
  Selectivity_filter: '[TSD][TSD][DE]',
};

// Defaults from the Python signatures: predict_tm_domains(window=19, threshold=1.6, min_windows=1, merge_gap=5), the
// standard Kyte & Doolittle (1982) rule, and SequenceDatabase.filter_database(200, 2000, 2, 10).
export const DEFAULTS = { window: 19, threshold: 1.6, minWindows: 1, mergeGap: 5,
  minLength: 200, maxLength: 2000, minTm: 2, maxTm: 10 };

// numpy's float64 pairwise summation (8 accumulators up to 128 items, then halving), so a window mean lands on exactly
// the same double as np.mean and threshold ties resolve the same way as in Python (checked bitwise in the node test).
function pairwise(a, lo, n) {
  if (n < 8) { let r = -0; for (let i = 0; i < n; i++) r += a[lo + i]; return r; }
  if (n <= 128) {
    const r = a.slice(lo, lo + 8);
    let i = 8;
    for (; i < n - (n % 8); i += 8) for (let j = 0; j < 8; j++) r[j] += a[lo + i + j];
    let res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
    for (; i < n; i++) res += a[lo + i];
    return res;
  }
  let n2 = Math.floor(n / 2); n2 -= n2 % 8;
  return pairwise(a, lo, n2) + pairwise(a, lo + n2, n - n2);
}
export function npMean(a) {
  if (!a.length) return NaN;
  return pairwise(a, 0, a.length) / a.length;
}

// Mean hydropathy of every window (hydro_profile in predict_tm_domains).
export function hydroProfile(sequence, window = DEFAULTS.window) {
  const s = sequence.toUpperCase(), vals = [...s].map(c => KD[c] ?? 0), out = [];
  for (let i = 0; i + window <= vals.length; i++) out.push(npMean(vals.slice(i, i + window)));
  return out;
}

// TransmembraneDomainPredictor.predict_tm_domains: [start, end) residue ranges (0-based, end exclusive).
export function predictTmDomains(sequence, window = DEFAULTS.window, threshold = DEFAULTS.threshold,
  minWindows = DEFAULTS.minWindows, mergeGap = DEFAULTS.mergeGap) {
  if (sequence.length < window) return [];
  const prof = hydroProfile(sequence, window), tm = [];
  let inTm = false, start = 0;
  prof.forEach((h, i) => {
    if (h >= threshold && !inTm) { start = i; inTm = true; }
    else if (h < threshold && inTm) {
      if (i - start >= minWindows) tm.push([start, i - 1 + window]);
      inTm = false;
    }
  });
  if (inTm && prof.length - start >= minWindows) tm.push([start, prof.length - 1 + window]);
  if (tm.length < 2) return tm;
  const merged = [tm[0]];
  for (const cur of tm.slice(1)) {
    const last = merged[merged.length - 1];
    if (cur[0] - last[1] <= mergeGap) merged[merged.length - 1] = [last[0], cur[1]];
    else merged.push(cur);
  }
  return merged;
}

// MotifScanner.scan_sequence: start positions of non-overlapping matches (re.finditer, IGNORECASE).
export function scanSequence(sequence, pattern) {
  const re = new RegExp(pattern, 'gi'), out = [];
  let m;
  while ((m = re.exec(sequence)) !== null) { out.push(m.index); if (m[0] === '') re.lastIndex++; }
  return out;
}

export function scanAll(sequence, motifs = MOTIFS) {
  return Object.fromEntries(Object.entries(motifs).map(([k, p]) => [k, scanSequence(sequence, p)]));
}

// SequenceDatabase.load_from_fasta: lines stripped, header split on the first whitespace run.
export function parseFasta(text) {
  const out = [];
  let cur = null;
  for (const raw of String(text).split(/\r?\n/)) {
    const line = raw.trim();
    if (line.startsWith('>')) {
      if (cur) out.push(cur);
      const h = line.slice(1).trim(), k = h.search(/\s/);
      cur = { id: k < 0 ? h : h.slice(0, k), description: k < 0 ? null : h.slice(k).trim(), sequence: '' };
    } else if (cur) cur.sequence += line;
  }
  if (cur) out.push(cur);
  return out;
}

// Per-sequence analysis plus SequenceDatabase.filter_database (length filter, then TM-count filter).
export function analyse(records, opts = {}) {
  const o = { ...DEFAULTS, ...opts };
  return records.map(r => {
    const tm = predictTmDomains(r.sequence, o.window, o.threshold, o.minWindows, o.mergeGap), len = r.sequence.length;
    const lengthOk = o.minLength <= len && len <= o.maxLength, tmOk = o.minTm <= tm.length && tm.length <= o.maxTm;
    return { id: r.id, description: r.description, length: len, tm, n_tm: tm.length,
      motifs: scanAll(r.sequence), length_ok: lengthOk, tm_ok: tmOk, passes: lengthOk && tmOk };
  });
}
