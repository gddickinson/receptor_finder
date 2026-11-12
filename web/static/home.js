// Home page: hero, facts read from the repo code, and links to the pages.
import { route, app, api, esc, isStatic } from './core.js';

route('home', async () => {
  const [info, au] = await Promise.all([api('/info'), api('/audit')]);
  const f = info.filter_defaults, t = info.tm_defaults;
  app().innerHTML = `<section class="hero"><p class="kicker">Plant Ca²⁺ signalling · convergent evolution</p>
    <h1 class="glow">Receptor Finder</h1>
    <p class="lead">${esc(info.description)}. Plants show InsP₃- and cADPR-induced Ca²⁺ release but lack the animal
      InsP₃ and ryanodine receptors; this toolkit sets out a five-phase search for proteins that could have taken
      their place: sequence filters, structure prediction, pocket scoring, docking and a 100-point ranking.</p>
    <a class="btn" href="#/sequence">Run the sequence filter</a><a class="btn alt" href="#/scoring">Score a candidate</a></section>
    <div class="tiles">
      <div class="tile"><b>${info.n_motifs}</b><span>regex motifs scanned</span></div>
      <div class="tile"><b>${t.window} / ${t.threshold}</b><span>Kyte-Doolittle window / threshold</span></div>
      <div class="tile"><b>${f.min_tm}-${f.max_tm}</b><span>TM domains kept by the default filter</span></div>
      <div class="tile"><b>${info.n_score_fields}</b><span>criteria in the 100-point score</span></div>
      <div class="tile"><b>${info.n_species}</b><span>species in the phylogenetic classifier</span></div>
      <div class="tile"><b>0</b><span>real candidate lists committed</span></div></div>
    <div class="callout"><b>What this site is.</b> The repo is a methods scaffold: it contains no proteome search results
      and no candidate list. The three committed sequences are hypothetical examples, the committed structures are
      placeholders, and the docking and pocket-detection steps return mock values. This site therefore shows the
      method and the parts that are real code (the sequence filters, motif scanning and the scoring rules) running
      ${isStatic() ? 'in your browser' : 'in Python'}, plus an <a href="#/audit">audit</a> of what the rules can and
      cannot detect: for example, the shipped EF-hand pattern finds ${au.ef_hand.detected} of the
      ${au.ef_hand.loops.length} calmodulin EF-hand loops in the example data.</div>
    <div class="grid">
      <a class="card" href="#/method"><b>Method</b><span>The five phases, what each needs, and which parts are implemented or placeholders.</span></a>
      <a class="card" href="#/sequence"><b>Sequence filter</b><span>Paste FASTA: hydropathy profile, predicted TM segments, motif hits and the length / TM-count filter.</span></a>
      <a class="card" href="#/scoring"><b>Scoring</b><span>The 100-point rubric as a calculator, plus the InsP₃ / cADPR pocket rules and the phylogenetic pattern classifier.</span></a>
      <a class="card" href="#/audit"><b>Audit</b><span>Committed files and what they are; benchmarks of the motif, TM and phylogeny rules; ligand structure checks.</span></a>
      <a class="card" href="#/about"><b>About</b><span>How to run the desktop GUI and this web app; credits.</span></a>
    </div>`;
});
