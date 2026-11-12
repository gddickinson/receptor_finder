// Drawing helpers for the sequence page: hydropathy profile (SVG) and the residue map with TM segments and motif hits.
import { esc } from './core.js';
import { hydroProfile, scanSequence } from './seq_engine.js';

// Hydropathy profile: mean Kyte-Doolittle score of the window starting at each residue, threshold line, TM bands.
export function profileSvg(sequence, tm, window, threshold) {
  const prof = hydroProfile(sequence, window);
  if (!prof.length) return `<p class="mute small">Sequence shorter than the window (${window}): no profile.</p>`;
  const W = 760, H = 220, L = 40, R = 10, T = 10, B = 30, lo = -4.5, hi = 4.5;
  const n = sequence.length, x = i => L + (i / Math.max(1, n)) * (W - L - R), y = v => T + (hi - v) / (hi - lo) * (H - T - B);
  const pts = prof.map((v, i) => `${x(i + window / 2).toFixed(1)},${y(v).toFixed(1)}`).join(' ');
  const grid = [-4, -2, 0, 2, 4].map(v => `<line class="grid-l" x1="${L}" x2="${W - R}" y1="${y(v)}" y2="${y(v)}"/>
    <text x="${L - 6}" y="${y(v) + 4}" text-anchor="end">${v}</text>`).join('');
  const ticks = [0, .25, .5, .75, 1].map(f => Math.round(f * n)).map(i => `<text x="${x(i)}" y="${H - 10}" text-anchor="middle">${i}</text>`).join('');
  const bands = tm.map(([a, b]) => `<rect class="tmband" x="${x(a)}" y="${T}" width="${x(b) - x(a)}" height="${H - T - B}"/>`).join('');
  return `<svg class="plot" viewBox="0 0 ${W} ${H}" role="img" aria-label="Hydropathy profile with ${tm.length} predicted TM segments">
    ${grid}${bands}<line class="thr" x1="${L}" x2="${W - R}" y1="${y(threshold)}" y2="${y(threshold)}"/>
    <polyline class="prof" points="${pts}"/>${ticks}<text x="${W - R}" y="${y(threshold) - 5}" text-anchor="end">threshold ${threshold}</text></svg>
    <p class="small mute">Window mean of the Kyte-Doolittle scale plotted at the window centre; shaded: predicted TM
      segments [start, end) in residue numbers from 0.</p>`;
}

// Residue map: TM residues shaded, residues covered by a match of the chosen motif underlined.
export function residueMap(sequence, tm, pattern) {
  const inTm = new Uint8Array(sequence.length), inMf = new Uint8Array(sequence.length);
  tm.forEach(([a, b]) => inTm.fill(1, a, Math.min(b, sequence.length)));
  if (pattern) {
    const re = new RegExp(pattern, 'gi');
    let m;
    while ((m = re.exec(sequence)) !== null) { inMf.fill(1, m.index, m.index + m[0].length); if (!m[0]) re.lastIndex++; }
  }
  let out = '';
  for (let i = 0; i < sequence.length; i++) {
    const c = [inTm[i] ? 'tm' : '', inMf[i] ? 'mf' : ''].filter(Boolean).join(' ');
    out += (i && i % 10 === 0 ? ' ' : '') + (c ? `<span class="${c}">${esc(sequence[i])}</span>` : esc(sequence[i]));
  }
  return `<div class="seqmap" aria-label="Residue map">${out}</div>`;
}

export const motifCount = (seq, pattern) => scanSequence(seq, pattern).length;
