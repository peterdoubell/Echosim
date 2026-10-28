// tools/verify-report.mjs — checks for the structured echo report.
//
// Part 1 (Node, no browser): builds the report for every pathology and checks
// that each case's conclusion names its lesion, that the graded findings agree
// with the BSE 2020 bands in js/references.js, that indexing follows the patient
// (BSA, height, sex), and that every template fills without leaving a {KEY}.
// Part 2 (browser): opens the report dialog from the Measurements panel, checks
// the tables, findings and conclusion render, that the patient form re-grades
// the case, that a template inserts into the write-up, that challenge mode locks
// the report, and that nothing logs a console error.
//
// Usage: node tools/verify-report.mjs      (exit 0 = pass, 1 = fail)
import http from 'http';
import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';
import { buildReport, fillTemplate, TEMPLATES, reportText } from '../js/report.js';
import { gradeGraded, gradeLvef, gradeLavi, bsa } from '../js/references.js';

const results = [];
const check = (name, ok, note = '') => results.push({ name, ok: !!ok, note });

// ---- part 1: report content ------------------------------------------------------
const CASES = {
  normal: [{}, /structurally normal/i],
  as: [{ aorticStenosis: true, lvh: true, grade: 'severe' }, /severe aortic stenosis/i],
  mr: [{ mr: true, grade: 'severe' }, /severe mitral regurgitation/i],
  ms: [{ mitralStenosis: true, grade: 'severe' }, /severe mitral stenosis/i],
  dcm: [{ dilated: true }, /severely impaired lv systolic function/i],
  rwma: [{ rwma: 'septal' }, /regional wall motion abnormality/i],
  effusion: [{ effusion: true }, /pericardial effusion/i],
  vsd: [{ vsd: true }, /vsd/i],
  asd: [{ asd: true }, /secundum asd/i],
  tr: [{ tr: true, grade: 'severe' }, /tricuspid regurgitation/i],
  phtn: [{ rvpo: true, tr: true }, /high echocardiographic probability of pulmonary hypertension/i],
};
const P0 = { sex: 'male', age: 45, height: 175, weight: 75 };
for (const [k, [path, re]] of Object.entries(CASES)) {
  const rep = buildReport(path, P0);
  const concl = rep.conclusion.join(' ');
  check(`${k}: conclusion names the lesion`, re.test(concl), concl.slice(0, 120));
  // the conclusion of an abnormal case is not the normal-study sentence
  if (k !== 'normal') check(`${k}: not reported as a normal study`, !/structurally normal/i.test(concl));
  check(`${k}: every finding section present`, rep.findings.length >= 8, `${rep.findings.length} sections`);
  const text = reportText(rep, 'test write-up');
  check(`${k}: plain-text report renders`, /CONCLUSION/.test(text) && /LEARNER REPORT/.test(text) && !/undefined|NaN/.test(text));
}
// graded findings use the BSE bands of the patient's sex
{
  const rep = buildReport({}, P0), m = rep.m;
  const lv = rep.rows.find((g) => g.title === 'Left ventricle').rows;
  const row = (l) => lv.find((r) => r.label === l);
  check('LVIDd grade matches the BSE male band', row('LVIDd').level === gradeGraded('lvidd', m.lvidd, 'male').level);
  check('LVEF grade matches the BSE LVEF band', row('LVEF').word === gradeLvef(m.lvef).word, `${m.lvef.toFixed(0)}% -> ${row('LVEF').word}`);
  const B = bsa(175, 75);
  check('BSA (Mosteller) for 175 cm / 75 kg', Math.abs(rep.patient.bsa - 1.909) < 0.01, rep.patient.bsa.toFixed(3));
  const la = rep.rows.find((g) => g.title === 'Atria').rows[0];
  check('LA volume indexed to BSA', la.idx === `${(m.lavol / B).toFixed(0)} mL/m²` && la.word === gradeLavi(m.lavol / B).word, la.idx);
  // the same heart reported as a smaller woman: indexed values rise, female bands apply
  const f = buildReport({}, { sex: 'female', age: 45, height: 160, weight: 55 });
  const lvF = f.rows.find((g) => g.title === 'Left ventricle').rows.find((r) => r.label === 'LVIDd');
  check('female patient uses the female LVIDd band', lvF.ref === '35–51', lvF.ref);
  const aoM = rep.rows.find((g) => g.title === 'Aorta').rows[0], aoF = f.rows.find((g) => g.title === 'Aorta').rows[0];
  check('aortic dimensions indexed to height', parseFloat(aoF.idx) > parseFloat(aoM.idx), `${aoM.idx} (175 cm) vs ${aoF.idx} (160 cm)`);
  // DCM: LVEF quoted as a number with BSE wording, never "mild/moderate dysfunction"
  const d = buildReport({ dilated: true }, P0);
  const lvTxt = d.findings.find((x) => x.section === 'Left ventricle').text;
  check('DCM LVEF quoted as a number with BSE wording', /LVEF \d+% \(severely impaired\)/.test(lvTxt) && !/(mild|moderate)(ly)? (lv )?dysfunction/i.test(lvTxt), lvTxt.slice(0, 90));
  // AS: Vmax, mean gradient and AVA agree on severity (concordant grading)
  const a = buildReport({ aorticStenosis: true, lvh: true, grade: 'severe' }, P0).m;
  check('severe AS: Vmax >= 4, mean >= 40, AVA < 1.0 (concordant)', a.avVmax >= 4 && a.avMeanGrad >= 40 && a.ava < 1.0,
    `Vmax ${a.avVmax.toFixed(1)} m/s, mean ${a.avMeanGrad.toFixed(0)} mmHg, AVA ${a.ava.toFixed(2)} cm²`);
}
// templates fill every placeholder
{
  const rep = buildReport({ mitralStenosis: true, grade: 'severe' }, P0);
  for (const k of Object.keys(TEMPLATES)) {
    const t = fillTemplate(k, rep);
    check(`template "${k}" fills every {KEY}`, !/\{[A-Z0-9]+\}/.test(t), (t.match(/\{[A-Z0-9]+\}/g) || []).join(' '));
  }
}

// ---- part 2: the dialog in the browser ---------------------------------------------
let chromium;
try { ({ chromium } = await import('playwright')); } catch (e) { chromium = null; }
if (!chromium) {
  check('browser part ran (playwright importable)', false, 'playwright not found: npm install playwright');
} else {
  const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
  const MIME = { '.html': 'text/html', '.js': 'text/javascript', '.css': 'text/css', '.svg': 'image/svg+xml', '.png': 'image/png', '.webmanifest': 'application/manifest+json' };
  const PORT = +(process.env.PORT || 8096);
  const server = http.createServer((req, res) => {
    let p = decodeURIComponent(req.url.split('?')[0]);
    if (p === '/') p = '/index.html';
    fs.readFile(path.join(ROOT, p), (e, d) => {
      if (e) { res.writeHead(404); res.end('not found'); return; }
      res.writeHead(200, { 'Content-Type': MIME[path.extname(p)] || 'application/octet-stream' });
      res.end(d);
    });
  });
  await new Promise((r) => server.listen(PORT, r));
  const browser = await chromium.launch({
    executablePath: process.env.CHROME_PATH || [
      '/opt/pw-browsers/chromium', '/opt/pw-browsers/chromium-1194/chrome-linux/chrome',
    ].find((p) => fs.existsSync(p)),
    args: ['--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--no-sandbox'],
  });
  const page = await browser.newPage({ viewport: { width: 1440, height: 880 } });
  const errors = [];
  page.on('console', (m) => { if (m.type() === 'error' && !/favicon/.test(m.text())) errors.push(m.text()); });
  page.on('pageerror', (e) => errors.push(e.message));
  await page.goto(`http://localhost:${PORT}/`, { waitUntil: 'networkidle' });
  await page.waitForTimeout(1200);
  const got = await page.$('#helpGot'); if (got) await got.click().catch(() => {});
  await page.waitForTimeout(300);

  const openFor = async (key) => {
    await page.selectOption('#pathology', key);
    await page.waitForTimeout(400);
    await page.click('#reportBtn');
    await page.waitForFunction(() => document.querySelectorAll('#reportTables table').length > 0, null, { timeout: 15000 });
  };
  await openFor('as');
  const ui = await page.evaluate(() => ({
    tables: document.querySelectorAll('#reportTables table').length,
    rows: document.querySelectorAll('#reportTables tbody tr').length,
    findings: document.querySelectorAll('#reportFindings dt').length,
    concl: document.getElementById('reportConclusion').textContent,
    focusInside: document.getElementById('reportOverlay').contains(document.activeElement),
  }));
  check('dialog: measurement tables render', ui.tables >= 5 && ui.rows >= 20, `${ui.tables} tables, ${ui.rows} rows`);
  check('dialog: findings render', ui.findings >= 8, `${ui.findings} sections`);
  check('dialog: AS conclusion', /severe aortic stenosis/i.test(ui.concl), ui.concl.slice(0, 80));
  check('dialog: focus moves into the dialog', ui.focusInside);
  // patient form re-grades: female, 160 cm / 55 kg
  await page.selectOption('#rpSex', 'female');
  await page.fill('#rpHeight', '160'); await page.dispatchEvent('#rpHeight', 'change');
  await page.fill('#rpWeight', '55'); await page.dispatchEvent('#rpWeight', 'change');
  await page.waitForTimeout(300);
  const fem = await page.evaluate(() => ({ bsa: document.getElementById('rpBsa').textContent, src: document.getElementById('reportSrc').textContent }));
  check('dialog: BSA follows the patient form', fem.bsa === '1.56', `BSA ${fem.bsa}`);
  check('dialog: report re-graded for a female patient', /female/.test(fem.src));
  // template inserted into the write-up
  await page.click('[data-rtab="write"]');
  await page.click('#rpInsert');
  const wu = await page.$eval('#rpWriteUp', (el) => el.value);
  check('dialog: AS template inserted with values', /Vmax \d\.\d m\/s/.test(wu) && !/\{[A-Z0-9]+\}/.test(wu), wu.split('\n')[0].slice(0, 90));
  await page.click('[data-rtab="guide"]');
  const guideVisible = await page.$eval('[data-rpanel="guide"]', (el) => !el.classList.contains('hidden') && el.textContent.includes('Harkness'));
  check('dialog: reporting guide tab shows the cited guide', guideVisible);
  await page.click('[data-rtab="report"]');
  if (process.env.SHOT) await page.screenshot({ path: process.env.SHOT });
  // Escape closes and returns focus to the report button
  await page.keyboard.press('Escape');
  await page.waitForTimeout(200);
  const closed = await page.evaluate(() => document.getElementById('reportOverlay').classList.contains('hidden') && document.activeElement.id === 'reportBtn');
  check('dialog: Escape closes it and restores focus', closed);
  // the panel's normal range follows the patient's sex
  const ref = await page.$eval('#refLvidd', (el) => el.textContent);
  check('panel: LVIDd normal range follows the patient (female)', /3\.5–5\.1/.test(ref), ref);
  // restore a male default so later runs start clean
  await page.evaluate(() => localStorage.removeItem('echosim_patient'));
  // challenge mode locks the report (it would name the lesion)
  const exp = await page.$('.ctl-group h3:text-is("Training")');
  await page.evaluate(() => { for (const h of document.querySelectorAll('.ctl-group h3')) if (/^Training/.test(h.firstChild.textContent.trim()) && h.parentElement.classList.contains('collapsed')) h.click(); });
  await page.click('#challengeBtn');
  await page.waitForTimeout(300);
  const locked = await page.$eval('#reportBtn', (el) => el.disabled);
  check('challenge mode disables the report', locked);
  await page.click('#challengeBtn');
  await page.waitForTimeout(300);
  check('report re-enabled after the challenge', !(await page.$eval('#reportBtn', (el) => el.disabled)));
  check('no console errors', errors.length === 0, errors.slice(0, 3).join(' | '));
  await browser.close(); server.close();
}

for (const r of results) console.log(`  ${r.ok ? 'PASS' : 'FAIL'} ${r.name}${r.note ? '  [' + r.note + ']' : ''}`);
const failed = results.filter((r) => !r.ok).length;
console.log(`\n  ${results.length - failed} passed, ${failed} failed.`);
process.exit(failed ? 1 : 0);
