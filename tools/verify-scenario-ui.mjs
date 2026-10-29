// Scenario/control integration checks. Run with `node tools/verify-scenario-ui.mjs`.
// CHROME_PATH selects an installed browser; otherwise Playwright's browser is used.
import { chromium } from 'playwright';
import http from 'node:http';
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { hemoSummary } from '../js/cardiac-model.js';

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const OUT = path.join(ROOT, '.review/scenario-ui');
fs.mkdirSync(OUT, { recursive: true });
const rows = [];
const check = (name, ok, detail = '') => {
  rows.push({ name, ok: !!ok, detail });
  console.log(`  ${ok ? 'PASS' : 'FAIL'} ${name}${detail ? ` — ${detail}` : ''}`);
};
const mime = { '.html': 'text/html', '.js': 'text/javascript', '.css': 'text/css', '.svg': 'image/svg+xml', '.png': 'image/png', '.webmanifest': 'application/manifest+json' };
const server = http.createServer((req, res) => {
  const name = decodeURIComponent(req.url.split('?')[0]);
  const file = path.join(ROOT, name === '/' ? 'index.html' : name);
  fs.readFile(file, (err, data) => {
    res.writeHead(err ? 404 : 200, { 'Content-Type': mime[path.extname(file)] || 'application/octet-stream' });
    res.end(err ? 'not found' : data);
  });
});
await new Promise((resolve) => server.listen(0, '127.0.0.1', resolve));
let browser;
try {
  browser = await chromium.launch({
    executablePath: process.env.CHROME_PATH || [
      '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',
      '/opt/pw-browsers/chromium',
      '/opt/pw-browsers/chromium-1194/chrome-linux/chrome',
    ].find((p) => fs.existsSync(p)),
    args: ['--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--no-sandbox'],
  });
  const page = await browser.newPage({ viewport: { width: 1440, height: 900 } });
  const errors = [];
  page.on('pageerror', (err) => errors.push(err.message));
  page.on('console', (msg) => { if (msg.type() === 'error' && !/favicon/.test(msg.text())) errors.push(msg.text()); });
  page.on('dialog', (dialog) => dialog.accept());
  await page.addInitScript(() => {
    localStorage.setItem('echosim_seen', '1');
    localStorage.setItem('echosim_panels_v1', JSON.stringify({ Training: false }));
  });
  await page.goto(`http://127.0.0.1:${server.address().port}/`, { waitUntil: 'networkidle' });
  await page.waitForFunction(() => window.echosim?.echo.metrics.ef != null);

  async function phase(value) {
    await page.$eval('#scrub', (el, value) => { el.value = String(value); el.dispatchEvent(new Event('input', { bubbles: true })); }, value);
    await page.evaluate(() => new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve))));
  }
  async function settledAnatomy() {
    await page.locator('#rebuildOverlay').waitFor({ state: 'hidden', timeout: 60000 });
    // The bake binds frames before hiding its overlay; let the live renderer
    // paint those frames before capturing proof of the completed result.
    await page.evaluate(() => new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve))));
  }
  const snapshot = () => page.evaluate(() => ({
    path: { ...window.echosim.state.path },
    ef: Number.parseFloat(document.getElementById('mEf').textContent),
    hemoEf: Number.parseFloat(document.getElementById('hEf').textContent),
    gls: document.getElementById('sGls').textContent,
    pattern: document.getElementById('sPattern').textContent,
    bull: document.getElementById('bullseye').toDataURL(),
    pv: document.getElementById('pvLoop').toDataURL(),
    ecg: document.getElementById('twelveLead').toDataURL(),
  }));

  // Every cardiac phase must occupy the same header slot at desktop/tablet/mobile widths.
  for (const width of [1440, 1200, 1100, 768, 390, 320]) {
    await page.setViewportSize({ width, height: 900 });
    const boxes = [];
    for (const p of [0, 200, 460, 580, 740, 900]) {
      await phase(p);
      boxes.push(await page.evaluate(() => {
        const box = (id) => { const r = document.getElementById(id).getBoundingClientRect(); return [r.x, r.y, r.width, r.height]; };
        const name = document.getElementById('phaseName');
        // A fixed-width span can still hide text motion caused by right or
        // centre alignment. Measure the rendered glyphs, not just their slot.
        const range = document.createRange();
        range.selectNodeContents(name);
        const text = range.getBoundingClientRect();
        return { ecg: box('ecg'), button: box('twelveBtn'), name: box('phaseName'), textX: text.x, label: name.textContent };
      }));
    }
    const stable = boxes.every((b) => ['ecg', 'button', 'name'].every((key) => b[key].every((v, i) => Math.abs(v - boxes[0][key][i]) < 1)));
    const textStartStable = boxes.every((b) => Math.abs(b.textX - boxes[0].textX) < 1 && Math.abs(b.textX - b.name[0]) < 1);
    const fits = await page.evaluate(() => document.querySelector('.topbar').scrollWidth <= document.documentElement.clientWidth);
    check(`${width}px header does not move across all six phase labels`, stable, boxes.map((b) => `${b.label}: ${b.button[0].toFixed(1)}`).join(', '));
    check(`${width}px phase text keeps a fixed left-aligned starting position`, textStartStable, boxes.map((b) => `${b.label}: ${b.textX.toFixed(1)}`).join(', '));
    check(`${width}px header fits the viewport`, fits);
  }
  await page.screenshot({ path: path.join(OUT, 'mobile-header.png'), fullPage: true });
  await page.setViewportSize({ width: 1440, height: 900 });

  // Preset changes must change the solved model AND each cached display.
  await page.selectOption('#pathology', 'dcm');
  check('DCM exposes EF selector and hides region selector', await page.locator('#dcmSeverity').isVisible() && !(await page.locator('#rwmaRegion').isVisible()));
  const dcms = [];
  for (const severity of ['severe', 'impaired', 'borderline', 'severe']) {
    await page.selectOption('#dcmSeverity', severity);
    await phase(250);
    const s = await snapshot(), expected = Math.round(hemoSummary({ dilated: true, dcmSeverity: severity }).EF);
    check(`DCM ${severity}: model, measured EF and haemodynamics agree`, s.path.dcmSeverity === severity && Math.abs(s.ef - expected) <= 1 && s.hemoEf === expected, `echo ${s.ef}%, haemodynamics ${s.hemoEf}%, model ${expected}%`);
    dcms.push(s);
  }
  check('DCM severity changes the pressure-volume loop and strain bullseye', new Set(dcms.slice(0, 3).map((s) => s.pv)).size === 3 && new Set(dcms.slice(0, 3).map((s) => s.bull)).size === 3);
  check('returning to severe DCM restores its EF and regional strain', dcms[0].ef === dcms[3].ef && dcms[0].bull === dcms[3].bull);
  await settledAnatomy();
  await page.screenshot({ path: path.join(OUT, 'dcm-severity.png') });

  await page.selectOption('#pathology', 'rwma');
  check('RWMA exposes region selector and hides EF selector', await page.locator('#rwmaRegion').isVisible() && !(await page.locator('#dcmSeverity').isVisible()));
  await page.click('#twelveBtn');
  const regions = [];
  for (const region of ['septal', 'anteroseptal', 'anterior', 'lateral', 'inferior', 'inferolateral']) {
    // The open ECG dialog covers the controls: close, change, and reopen it.
    await page.click('#twelveClose');
    await page.selectOption('#rwmaRegion', region);
    await phase(250);
    await page.click('#twelveBtn');
    await page.evaluate(() => new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve))));
    const s = await snapshot();
    check(`RWMA ${region}: selected region reaches the model`, s.path.rwma === region);
    regions.push(s);
  }
  check('all six RWMA presets produce distinct strain bullseyes', new Set(regions.map((s) => s.bull)).size === 6);
  check('all six RWMA presets refresh the matching 12-lead ECG while paused', new Set(regions.map((s) => s.ecg)).size === 6);
  await page.screenshot({ path: path.join(OUT, 'regional-ecg.png') });
  await page.click('#twelveClose');
  await settledAnatomy();
  await page.screenshot({ path: path.join(OUT, 'regional-control.png') });

  // Unknown-case workflows must hide added diagnostic controls and restore the prior preset.
  for (const [key, selector, value, field] of [
    ['rwma', '#rwmaRegion', 'inferior', 'rwma'],
    ['dcm', '#dcmSeverity', 'impaired', 'dcmSeverity'],
  ]) {
    await page.selectOption('#pathology', key);
    await page.selectOption(selector, value);
    await page.click('#pathwayBtn');
    await page.click('#pwNav [data-nav="lvfunc"]');
    await page.click('#pwMain [data-mtab="test"]');
    await page.click('#pwStartModule');
    let reachedCase = false;
    for (let i = 0; i < 12; i++) {
      if (await page.locator('#pwWatch').count()) { reachedCase = true; break; }
      if (!(await page.locator('#pwMain .pw-opt').count())) break;
      await page.locator('#pwMain .pw-opt').first().click();
      await page.click('#pwNext');
    }
    check(`${key}: pathway test reaches a case question`, reachedCase);
    if (reachedCase) {
      const controlsHidden = await page.evaluate(() => ['rwmaRegion', 'dcmSeverity'].every((id) => !document.getElementById(id).checkVisibility()));
      check(`${key}: pathway case hides both diagnostic selectors`, controlsHidden);
    }
    await page.click('#pwNav [data-nav="overview"]');
    await page.click('#pathwayClose');
    const restored = await page.evaluate((field) => window.echosim.state.path[field], field);
    check(`${key}: leaving pathway test restores selected subtype`, restored === value && await page.locator(selector).inputValue() === value, `${field}=${restored}`);
  }
  await page.click('#challengeBtn');
  check('challenge hides region and EF selectors', !(await page.locator('#rwmaRegion').isVisible()) && !(await page.locator('#dcmSeverity').isVisible()));
  await page.click('#challengeBtn');
  await page.selectOption('#pathology', 'normal');
  check('normal heart hides both scenario selectors', !(await page.locator('#rwmaRegion').isVisible()) && !(await page.locator('#dcmSeverity').isVisible()));
  check('no browser runtime errors', errors.length === 0, errors.slice(0, 4).join(' | '));
} catch (error) {
  check('test run completes', false, error.stack || error.message);
} finally {
  if (browser) await browser.close();
  await new Promise((resolve) => server.close(resolve));
}
fs.writeFileSync(path.join(OUT, 'results.json'), JSON.stringify(rows, null, 2));
const failures = rows.filter((r) => !r.ok).length;
console.log(`\n${rows.length - failures} passed, ${failures} failed. Screenshots: ${OUT}`);
process.exitCode = failures ? 1 : 0;
