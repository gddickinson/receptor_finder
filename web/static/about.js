// About page: what the project is, how to run the desktop GUI and the web app, credits (from README / CITATION.cff).
import { route, app } from './core.js';

route('about', () => {
  app().innerHTML = `<p class="kicker">About</p><h1>Receptor Finder</h1>
    <p class="lead">A computational toolkit for discovering candidate InsP₃ and cADPR receptors in plant genomes through
      convergent evolution analysis, written as Project #1 of the computational plant biology research proposals in
      <code>docs/</code>.</p>
    <h2>Run the desktop GUI</h2>
    <pre>pip install -r requirements.txt
python launch_gui.py</pre>
    <p>A Tkinter app with tabs for sequences, structures, pockets, docking, scoring and plots (docs/QUICKSTART_GUI.md).</p>
    <h2>Run this web app locally</h2>
    <pre>pip install -r requirements-web.txt
python run_web.py        # http://127.0.0.1:8206</pre>
    <p>The local app runs the sequence filter and the scoring through the repo's Python code. The static site runs
      JavaScript ports of the same functions, checked case by case against the Python in the test suite. Neither calls
      any external database or service; BLAST, HMMER, ColabFold, fpocket and AutoDock Vina are external tools you run
      yourself (the code prints their commands).</p>
    <h2>Limitations (from the README)</h2>
    <ul><li>AlphaFold2 may struggle with multi-domain proteins; truly novel folds may be missed.</li>
      <li>Docking scores are approximate predictions, not measurements.</li>
      <li>Phylogenetic inference needs broad taxon sampling.</li>
      <li>High computational scores do not guarantee function: experimental validation is essential.</li></ul>
    <h2>Credits</h2>
    <p>George D. Dickinson (University of California, Irvine). Methods based on AlphaFold2 (Jumper et al. 2021), fpocket (Le Guilloux et al. 2009) and
      AutoDock Vina (Trott &amp; Olson 2010). Built with AI assistance from Claude (Anthropic). MIT licence.</p>`;
});
