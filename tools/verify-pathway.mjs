// tools/verify-pathway.mjs — checks for the BSE accreditation learning pathway.
//
// Part 1 (Node): the question bank is well formed (valid answers, unique
// options, explanations), case questions and practice tasks point at real
// pathologies, views and grades, option shuffling keeps the right answer,
// the mock exam spreads across every module, scoring is correct, and numbers
// quoted in explanations agree with what the simulator's report measures.
// Part 2 (browser): opens the pathway, runs a practice task, passes a module
// test answering case questions from the case bar (checking the simulator
// hides the diagnosis and is restored afterwards), sits a short mock exam,
// and checks that progress survives a reload with no console errors.
//
// Usage: node tools/verify-pathway.mjs      (exit 0 = pass, 1 = fail)
import http from 'http';
import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';
import { MODULES, MODULE_BY_ID, PASS_MARK, EXAM, prepare, moduleQuiz, mockExam, score, weakModules, shuffle } from '../js/pathway.js';
import { TTE_VIEWS, EXTRA_VIEWS, TEE_VIEWS } from '../js/views.js';
import { buildReport } from '../js/report.js';

const results = [];
const check = (name, ok, note = '') => results.push({ name, ok: !!ok, note });
const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const html = fs.readFileSync(path.join(ROOT, 'index.html'), 'utf8');
const pathSelect = html.slice(html.indexOf('<select id="pathology"'), html.indexOf('</select>', html.indexOf('<select id="pathology"')));
const PATHS = [...pathSelect.matchAll(/<option value="(\w+)">/g)].map((m) => m[1]);
const VIEWS = new Set([...Object.keys(TTE_VIEWS), ...Object.keys(EXTRA_VIEWS), ...Object.keys(TEE_VIEWS)]);
const GRADES = new Set(['mild', 'moderate', 'severe']);

// seeded PRNG so failures reproduce
const rng = (seed) => () => { seed = (seed * 1664525 + 1013904223) % 4294967296; return seed / 4294967296; };

// ---- part 1: content and logic --------------------------------------------------------
check('ten modules', MODULES.length === 10, MODULES.map((m) => m.id).join(','));
const ids = new Set(), problems = [];
let total = 0, cases = 0;
for (const m of MODULES) {
  if (m.questions.length < 8) problems.push(`${m.id}: only ${m.questions.length} questions`);
  if (!m.objectives.length || !m.learn.length || !m.tasks.length) problems.push(`${m.id}: missing objectives, key points or tasks`);
  for (const q of m.questions) {
    total++;
    if (ids.has(q.id)) problems.push(`duplicate id ${q.id}`);
    ids.add(q.id);
    if (q.options.length !== 4) problems.push(`${q.id}: ${q.options.length} options`);
    if (new Set(q.options).size !== q.options.length) problems.push(`${q.id}: duplicate options`);
    if (!(Number.isInteger(q.answer) && q.answer >= 0 && q.answer < q.options.length)) problems.push(`${q.id}: bad answer index`);
    if (!q.explain || q.explain.length < 20) problems.push(`${q.id}: explanation missing`);
    if (q.case) {
      cases++;
      if (!PATHS.includes(q.case.path)) problems.push(`${q.id}: unknown pathology ${q.case.path}`);
      if (!VIEWS.has(q.case.view)) problems.push(`${q.id}: unknown view ${q.case.view}`);
      if (q.case.grade && !GRADES.has(q.case.grade)) problems.push(`${q.id}: bad grade`);
    }
  }
  for (const t of m.tasks) {
    if (!PATHS.includes(t.sim.path)) problems.push(`${m.id}/${t.id}: unknown pathology ${t.sim.path}`);
    if (!VIEWS.has(t.sim.view)) problems.push(`${m.id}/${t.id}: unknown view ${t.sim.view}`);
    if (t.sim.grade && !GRADES.has(t.sim.grade)) problems.push(`${m.id}/${t.id}: bad grade`);
  }
}
check('question bank is well formed', problems.length === 0, problems.slice(0, 4).join(' | '));
check(`bank size: ${total} questions, ${cases} case questions`, total >= 90 && cases >= 12);
check('pathology list read from index.html', PATHS.length === 11, PATHS.join(','));
check('every pathology appears in a case question', PATHS.every((p) => MODULES.some((m) => m.questions.some((q) => q.case && q.case.path === p))), PATHS.filter((p) => !MODULES.some((m) => m.questions.some((q) => q.case && q.case.path === p))).join(','));

// shuffling keeps the right answer
let shuffleOk = true, positions = new Set();
const r1 = rng(7);
for (const m of MODULES) for (const q of m.questions) for (let k = 0; k < 5; k++) {
  const p = prepare(q, r1);
  if (p.shown[p.correct] !== q.options[q.answer] || p.shown.length !== q.options.length) shuffleOk = false;
  positions.add(p.correct);
}
check('option shuffling keeps the correct answer', shuffleOk);
check('correct answers land in every position', positions.size === 4);
const mq = moduleQuiz('aortic', rng(3));
check('module quiz has every question once', mq.length === MODULE_BY_ID.aortic.questions.length && new Set(mq.map((q) => q.id)).size === mq.length);
const ex = mockExam(EXAM.questions, rng(11));
const perMod = {};
ex.forEach((q) => { perMod[q.module] = (perMod[q.module] || 0) + 1; });
check(`mock exam: ${EXAM.questions} unique questions`, ex.length === EXAM.questions && new Set(ex.map((q) => q.id)).size === ex.length);
check('mock exam covers every module evenly', MODULES.every((m) => perMod[m.id] >= 3 && perMod[m.id] <= 5), JSON.stringify(perMod));
check('mock exams differ between attempts', mockExam(EXAM.questions, rng(12)).map((q) => q.id).join() !== ex.map((q) => q.id).join());

// scoring
const allRight = score(ex, ex.map((q) => q.correct));
check('all correct scores 100% and passes', allRight.pct === 1 && allRight.passed && allRight.correct === EXAM.questions);
const half = score(ex, ex.map((q, i) => (i % 2 ? q.correct : (q.correct + 1) % 4)));
check('half correct fails', half.correct === EXAM.questions / 2 && !half.passed);
const blank = score(ex, ex.map(() => null));
check('unanswered questions score zero', blank.correct === 0);
const atPass = Math.ceil(EXAM.questions * PASS_MARK);
check(`${atPass}/${EXAM.questions} is a pass`, score(ex, ex.map((q, i) => (i < atPass ? q.correct : null))).passed);
check(`${atPass - 1}/${EXAM.questions} is not`, !score(ex, ex.map((q, i) => (i < atPass - 1 ? q.correct : null))).passed);
const aorticWrong = score(ex, ex.map((q) => (q.module === 'aortic' ? null : q.correct)));
check('weak modules: the failed topic is flagged', weakModules(aorticWrong).join() === 'aortic');
check('shuffle does not lose items', shuffle([1, 2, 3, 4, 5], rng(1)).sort().join() === '1,2,3,4,5');

// explanations quote the simulator's own measurements
const P0 = { sex: 'male', age: 45, height: 175, weight: 75 };
const rep = (flags) => buildReport({ ...flags, grade: flags.grade || 'severe' }, P0).m;
const quoted = (id) => MODULES.flatMap((m) => m.questions).find((q) => q.id === id).explain;
const findQ = (mid, pathKey) => MODULE_BY_ID[mid].questions.find((q) => q.case && q.case.path === pathKey).id;
const near = (txt, re, v, tol) => { const m = txt.match(re); return m && Math.abs(+m[1] - v) <= tol; };
check('DCM case quotes the model LVEF', near(quoted(findQ('lvfunc', 'dcm')), /LVEF ≈ (\d+)%/, rep({ dilated: true }).lvef, 1));
check('RWMA case quotes the model LVEF', near(quoted(findQ('lvfunc', 'rwma')), /\(~(\d+)%\)/, rep({ rwma: 'septal' }).lvef, 1));
check('MR case quotes the model regurgitant fraction', near(quoted(findQ('mitral', 'mr')), /~(\d+)%/, rep({ mr: true }).mrRf, 1));
const ph = rep({ rvpo: true, tr: true });
check('PH case quotes the model TR Vmax and PASP', near(quoted(findQ('rightheart', 'phtn')), /TR Vmax ([\d.]+) m\/s/, ph.trVmax, 0.05) && near(quoted(findQ('rightheart', 'phtn')), /PASP ~(\d+)/, ph.pasp, 1));
check('Effusion case quotes the model depth', near(quoted(findQ('pericongenital', 'effusion')), /measures ([\d.]+) cm/, rep({ effusion: true }).effusion, 0.05));
const as = rep({ aorticStenosis: true, lvh: true });
const asQ = MODULE_BY_ID.aortic.questions.find((q) => /^Peak velocity/.test(q.q)).q;
check('AS grading question matches the model', near(asQ, /Peak velocity ([\d.]+) m\/s/, as.avVmax, 0.05)
  && near(asQ, /mean gradient (\d+) mmHg/, as.avMeanGrad, 1) && near(asQ, /AVA ([\d.]+) cm²/, as.ava, 0.01));

// ---- part 2: the dialog in the browser ---------------------------------------------
let chromium;
try { ({ chromium } = await import('playwright')); } catch (e) { chromium = null; }
if (!chromium) {
  check('browser part ran (playwright importable)', false, 'playwright not found: npm install playwright');
} else {
  const MIME = { '.html': 'text/html', '.js': 'text/javascript', '.css': 'text/css', '.svg': 'image/svg+xml', '.png': 'image/png', '.webmanifest': 'application/manifest+json' };
  const PORT = +(process.env.PORT || 8098);
  const server = http.createServer((req, res) => {
    let p = decodeURIComponent(req.url.split('?')[0]);
    if (p === '/') p = '/index.html';
    fs.readFile(path.join(ROOT, p), (e, d) => {
      if (e) { res.writeHead(404); res.end('not found'); return; }
      res.writeHead(200, { 'Content-Type': MIME[path.extname(p)] || 'application/octet-stream' });
      res.end(d);
    });
  });
  await new Promise((ok) => server.listen(PORT, ok));
  const browser = await chromium.launch({
    executablePath: process.env.CHROME_PATH || [
      '/opt/pw-browsers/chromium', '/opt/pw-browsers/chromium-1194/chrome-linux/chrome',
    ].find((p) => fs.existsSync(p)),
    args: ['--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--no-sandbox'],
  });
  const ctx = await browser.newContext({ viewport: { width: 1440, height: 900 } });
  const page = await ctx.newPage();
  const errors = [];
  page.on('console', (msg) => { if (msg.type() === 'error' && !/favicon/.test(msg.text())) errors.push(msg.text()); });
  page.on('pageerror', (e) => errors.push(e.message));
  let dialogs = 0;
  page.on('dialog', (d) => { dialogs++; d.accept(); });
  await page.goto(`http://localhost:${PORT}/`, { waitUntil: 'networkidle' });
  await page.waitForTimeout(1200);
  const got = await page.$('#helpGot'); if (got) await got.click().catch(() => {});
  await page.waitForTimeout(300);
  const shot = async (name) => { if (process.env.SHOT) await page.screenshot({ path: `${process.env.SHOT}-${name}.png` }); };

  // the bank, so the test can pick the right answer from what is on screen
  const BANK = MODULES.flatMap((m) => m.questions.map((q) => ({ q: q.q, options: q.options, right: q.options[q.answer] })));
  const rightAnswer = async (sel) => page.evaluate(([bank, sel]) => {
    const q = document.querySelector(sel.q).textContent.trim();
    const shown = [...document.querySelectorAll(sel.opts)].map((b) => b.textContent.slice(1).trim());
    const hit = bank.find((x) => x.q === q && x.options.length === shown.length && x.options.every((o) => shown.includes(o)));
    return hit ? shown.indexOf(hit.right) : -1;
  }, [BANK, sel]);
  const MAIN = { q: '#pwMain .pw-q', opts: '#pwMain .pw-opt' }, BAR = { q: '#pwCbQ', opts: '#pwCbOpts .pw-opt' };

  await page.click('#pathwayBtn');
  const ov = await page.evaluate(() => ({
    open: !document.getElementById('pathwayOverlay').classList.contains('hidden'),
    nav: document.querySelectorAll('#pwNav .pw-nav-item').length,
    disclaimer: /not produced by or affiliated with the British Society of Echocardiography/.test(document.getElementById('pwMain').textContent),
  }));
  check('pathway opens from the header', ov.open);
  check('nav lists overview, 10 modules and the mock exam', ov.nav === 12, String(ov.nav));
  check('overview shows the independence disclaimer', ov.disclaimer);
  await shot('overview');

  // learn + practise
  await page.click('#pwNav [data-nav="aortic"]');
  check('module Learn tab shows key points', (await page.$$('#pwMain .pw-key li')).length >= 5);
  await page.click('#pwMain [data-mtab="practise"]');
  await page.click('#pwMain [data-open]');
  await page.waitForTimeout(400);
  const task = await page.evaluate(() => ({
    closed: document.getElementById('pathwayOverlay').classList.contains('hidden'),
    path: document.getElementById('pathology').value,
    grade: document.querySelector('#severityCtl .seg-btn.active')?.dataset.grade,
    view: document.getElementById('viewTag').textContent,
  }));
  check('practice task loads its case and closes the pathway', task.closed && task.path === 'as' && task.grade === 'moderate' && task.view === 'A5C', JSON.stringify(task));

  // module test with case questions, answered from the case bar
  await page.click('#pathwayBtn');
  await page.click('#pwNav [data-nav="views"]');
  await page.click('#pwMain [data-mtab="test"]');
  await page.click('#pwStartModule');
  const n = MODULE_BY_ID.views.questions.length;
  let hiddenOk = true, barOk = true, viewHidden = true, caseSeen = 0;
  for (let i = 0; i < n; i++) {
    const isCase = !!(await page.$('#pwWatch'));
    if (isCase) {
      caseSeen++;
      const st = await page.evaluate(() => ({
        sel: document.getElementById('pathology').selectedIndex, dis: document.getElementById('pathology').disabled,
        rep: document.getElementById('reportBtn').disabled, chal: document.getElementById('challengeBtn').disabled,
        tag: document.getElementById('viewTag').style.visibility, act: document.querySelectorAll('[data-view].active').length,
        labels: document.getElementById('labelsOn').checked,
      }));
      if (!(st.sel === -1 && st.dis && st.rep && st.chal && !st.labels)) hiddenOk = false;
      if (!(st.tag === 'hidden' && st.act === 0)) viewHidden = false;
      await page.click('#pwWatch');
      const vis = await page.evaluate(() => !document.getElementById('pwCaseBar').classList.contains('hidden') && document.getElementById('pathwayOverlay').classList.contains('hidden'));
      if (!vis) barOk = false;
      if (caseSeen === 1) await shot('casebar');
      const k = await rightAnswer(BAR);
      await page.click(`#pwCbOpts [data-cbopt="${k}"]`);
    } else {
      const k = await rightAnswer(MAIN);
      await page.click(`#pwMain [data-opt="${k}"]`);
    }
    if (i === 0) await shot('feedback');
    await page.click('#pwNext');
  }
  check('case questions were included', caseSeen === MODULE_BY_ID.views.questions.filter((q) => q.case).length, String(caseSeen));
  check('a case question hides the diagnosis and labels and locks the report', hiddenOk);
  check('"name this view" cases hide the view name', viewHidden);
  check('Watch the case swaps the dialog for the case bar', barOk);
  const res = await page.evaluate(() => ({ score: document.querySelector('.pw-score')?.textContent || '', status: document.querySelector('#pwNav [data-nav="views"] .pw-st')?.textContent }));
  check('all-correct module test passes', new RegExp(`${n}/${n}`).test(res.score) && /pass/.test(res.score), res.score);
  check('module marked passed in the nav', res.status === '✓');
  const restored = await page.evaluate(() => ({ labels: document.getElementById('labelsOn').checked, path: document.getElementById('pathology').value, dis: document.getElementById('pathology').disabled, rep: document.getElementById('reportBtn').disabled, tag: document.getElementById('viewTag').style.visibility }));
  check('simulator restored after the test', restored.labels && restored.path === 'as' && !restored.dis && !restored.rep && restored.tag !== 'hidden', JSON.stringify(restored));
  await shot('result');

  // leaving a test part-way asks first
  await page.click('#pwNav [data-nav="physics"]');
  await page.click('#pwMain [data-mtab="test"]');
  await page.click('#pwStartModule');
  const d0 = dialogs;
  await page.click('#pwNav [data-nav="overview"]');
  check('leaving a test part-way asks for confirmation', dialogs === d0 + 1 && !!(await page.$('#pwMain h3')));

  // mock exam: timer, case-question Escape, finish with blanks
  await page.click('#pwNav [data-nav="exam"]');
  await page.click('#pwStartExam');
  const timer = await page.$eval('.pw-timer', (el) => el.textContent);
  check('mock exam shows a running timer', /⏱ (59|60):\d\d/.test(timer), timer);
  check(`mock exam has ${EXAM.questions} questions`, (await page.$$('#pwMain .pw-dot')).length === EXAM.questions);
  // find a case question and press Escape on it
  let escOk = null;
  for (let i = 0; i < EXAM.questions && escOk === null; i++) {
    if (await page.$('#pwWatch')) {
      await page.keyboard.press('Escape');
      escOk = await page.evaluate(() => !document.getElementById('pwCaseBar').classList.contains('hidden'));
      await page.click('#pwCbBack');
    } else {
      const k = await rightAnswer(MAIN);
      await page.click(`#pwMain [data-opt="${k}"]`);
      if (await page.$('#pwNext')) await page.click('#pwNext');
    }
  }
  check('Escape on a case question minimises to the case bar', escOk === true);
  await shot('exam');
  await page.click('#pwFinish');
  await page.waitForTimeout(300);
  const exRes = await page.evaluate(() => ({ score: document.querySelector('.pw-score')?.textContent || '', rows: document.querySelectorAll('.pw-bymod tbody tr').length, review: document.querySelectorAll('.pw-rev').length }));
  check('finishing early scores the exam with a per-module breakdown', /\/40/.test(exRes.score) && exRes.rows === 10 && exRes.review === EXAM.questions, JSON.stringify(exRes));

  // persistence
  await page.reload({ waitUntil: 'networkidle' });
  await page.waitForTimeout(800);
  await page.click('#pathwayBtn');
  const persisted = await page.evaluate(() => ({
    views: document.querySelector('#pwNav [data-nav="views"] .pw-st')?.textContent,
    hist: document.querySelectorAll('.pw-history li').length,
    progress: document.getElementById('pwProgressText').textContent,
  }));
  check('progress and exam history survive a reload', persisted.views === '✓' && persisted.hist === 1 && /1 of 10 modules passed/.test(persisted.progress), JSON.stringify(persisted));
  await page.keyboard.press('Escape');
  check('Escape closes the pathway', await page.evaluate(() => document.getElementById('pathwayOverlay').classList.contains('hidden')));
  check('no console errors', errors.length === 0, errors.slice(0, 3).join(' | '));
  await browser.close(); server.close();
}

for (const x of results) console.log(`  ${x.ok ? 'PASS' : 'FAIL'} ${x.name}${x.note ? '  [' + x.note + ']' : ''}`);
const failed = results.filter((x) => !x.ok).length;
console.log(`\n  ${results.length - failed} passed, ${failed} failed.`);
process.exit(failed ? 1 : 0);
