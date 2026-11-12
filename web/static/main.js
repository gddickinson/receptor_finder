// Entry point: registers lazily loaded pages, sets the footer mode text and starts the router.
import { start, lazy, route, app, isStatic, $ } from './core.js';

const PAGES = { 'home.js': ['home'], 'method.js': ['method'], 'sequence.js': ['sequence'], 'scoring.js': ['scoring'],
  'audit.js': ['audit'], 'about.js': ['about'] };
lazy(Object.fromEntries(Object.entries(PAGES).flatMap(([file, names]) => names.map(n => [n, () => import('./' + file)]))));
route('notfound', () => { app().innerHTML = '<h1>Not found</h1><p><a href="#/">Back to the start</a></p>'; });

$('#mode').textContent = isStatic()
  ? 'static site: filters and scoring run in your browser (JS ports tested against the Python)'
  : 'local app: filters and scoring run in Python';
start();
