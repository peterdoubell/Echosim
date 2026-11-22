// tools/verify-measure.mjs — end-to-end check of the learner measurement tools.
//
// Loads the app, freezes at end-diastole/end-systole, drives the caliper and the
// Simpson biplane-disc EF tracer through real mouse clicks, and asserts that a
// distance, an EF and an accuracy-vs-truth figure all come out with zero console
// errors. The synthetic traces are approximate, so this validates the PIPELINE
// (tool state machine -> disc integration -> readout -> accuracy comparison), not
// clinical accuracy; the disc maths itself is exercised separately.
//
// This lives in tools/ rather than a scratch directory deliberately: an earlier
// copy sat outside the committed suite and silently stopped working when the
// console panels became collapsible — nothing ran it, so nothing caught it.
//
// Usage: node tools/verify-measure.mjs      (exit 0 = pass, 1 = fail)
import { chromium } from 'playwright';
import http from 'http';
import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const MIME = { '.html': 'text/html', '.js': 'text/javascript', '.css': 'text/css' };
const PORT = +(process.env.PORT || 8094);
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
await page.click('[data-view="A4C"]'); await page.waitForTimeout(800);

const overlayBox = () => page.$eval('#echoOverlay', (el) => {
  const r = el.getBoundingClientRect(); return { x: r.x, y: r.y, w: r.width, h: r.height };
});
async function setScrub(v) {
  await page.$eval('#scrub', (el, x) => {
    el.value = String(x); el.dispatchEvent(new Event('input', { bubbles: true }));
  }, v);
  await page.waitForTimeout(300);
}
// Trace hinge -> apex -> hinge: an open contour with both endpoints at the mitral
// plane and the apex at the far end, which is how the disc integrator expects a
// clinical trace to arrive.
async function clickTrace(cxFrac, cyFrac, wFrac, hFrac, n) {
  const box = await overlayBox();
  const cx = box.x + box.w * cxFrac, cy = box.y + box.h * cyFrac;
  const rx = box.w * wFrac, ry = box.h * hFrac, D = Math.PI / 180;
  const pt = (deg) => [cx + Math.cos(deg * D) * rx, cy + Math.sin(deg * D) * ry];
  for (let i = 0; i < n; i++) {
    const [x, y] = pt(115 + 310 * i / (n - 1));
    await page.mouse.click(x, y); await page.waitForTimeout(25);
  }
  const [lx, ly] = pt(425);
  await page.mouse.dblclick(lx, ly); await page.waitForTimeout(150);
}
// Console sections collapse, and a collapsed section's buttons are not clickable.
// Expand by activating the header, exactly as a user would.
async function expandPanel(title) {
  await page.evaluate((t) => {
    for (const h of document.querySelectorAll('.ctl-group h3')) {
      const label = h.firstChild && h.firstChild.textContent ? h.firstChild.textContent.trim() : '';
      if (label === t && h.parentElement.classList.contains('collapsed')) h.click();
    }
  }, title);
  await page.waitForTimeout(200);
}

await expandPanel('Caliper tools');

// caliper: two clicks -> a distance in cm
await setScrub(520);
await page.click('[data-tool="caliper"]'); await page.waitForTimeout(120);
{
  const box = await overlayBox();
  await page.mouse.click(box.x + box.w * 0.36, box.y + box.h * 0.34); await page.waitForTimeout(60);
  await page.mouse.click(box.x + box.w * 0.50, box.y + box.h * 0.34); await page.waitForTimeout(120);
}
// Simpson: trace end-diastole, then end-systole
await setScrub(520); await page.click('[data-tool="simpsonED"]'); await page.waitForTimeout(120);
await clickTrace(0.42, 0.45, 0.13, 0.26, 14);
await setScrub(200); await page.click('[data-tool="simpsonES"]'); await page.waitForTimeout(120);
await clickTrace(0.42, 0.45, 0.10, 0.20, 14);

const res = await page.evaluate(() => ({
  dist: document.getElementById('measDist').textContent,
  ef: document.getElementById('measSimpEf').textContent,
  tag: document.getElementById('measSimpTag').textContent,
  err: document.getElementById('measErr').textContent,
  truth: document.getElementById('mEf').textContent,
}));

// ---- trace panel: M-mode (TAPSE by a two-point caliper), tissue Doppler, IVC sniff ----
// The M-mode caliper is driven by real clicks on the trace at the annulus echo's
// deepest and shallowest points; it must return the excursion the trace shows,
// and that must be the model's TAPSE. Tissue Doppler must read s', e', a' and
// E/e'. In the subcostal IVC view the sniff trace must show the collapse.
const shot = async (name) => { if (process.env.SHOT_DIR) await page.locator('.echo-panel, #spectral').first().screenshot({ path: path.join(process.env.SHOT_DIR, name) }).catch(() => {}); };
const pickTrace = async (mode, line) => {
  await page.selectOption('#traceMode', mode); await page.waitForTimeout(150);
  if (line) { await page.selectOption('#traceLine', line); await page.waitForTimeout(150); }
};
await page.click('[data-view="A4C"]'); await page.waitForTimeout(600);
await pickTrace('mmode', 'tapse');
await page.waitForFunction(() => { const m = window.echosim.echo.mmodeData; return !!(m && m.key === 'tapse' && m.read); }, null, { timeout: 60000 });
const mm = await page.evaluate(() => { const m = window.echosim.echo.mmodeData; return { d0: m.d0, d1: m.d1, plotH: m.plotH, W: m.W, min: m.read.min, max: m.read.max, exc: m.excursion }; });
{
  const r = await page.$eval('#spectral', (el) => { const b = el.getBoundingClientRect(); return { x: b.x, y: b.y, w: b.width, h: b.height, H: el.height }; });
  const yOf = (d) => r.y + (d - mm.d0) / (mm.d1 - mm.d0) * mm.plotH * r.h / r.H;
  await page.mouse.click(r.x + r.w * 0.3, yOf(mm.min)); await page.waitForTimeout(80);
  await page.mouse.click(r.x + r.w * 0.6, yOf(mm.max)); await page.waitForTimeout(300);
}
const mmCal = await page.evaluate(() => (window.echosim.echo.mmodeData.caliper || {}).dDepth);
const tapseModel = await page.evaluate(async () => { const { measureModel } = await import('./js/report.js'); return measureModel(window.echosim.state.path).tapse / 10; });
await shot('mmode_tapse.png');
await pickTrace('tdi', 'mvLat');
await page.waitForFunction(() => { const t = window.echosim.echo.tdiData; return !!(t && t.site === 'mvLat'); }, null, { timeout: 60000 });
const tdi = await page.evaluate(() => { const t = window.echosim.echo.tdiData; return { s: t.s, e: t.e, a: t.a, EoverE: t.EoverE }; });
await shot('tdi_mvlat.png');
// the PLAX lines carry an automatic reading (LV dimensions; aortic root and LA)
await page.click('[data-view="PLAX"]'); await page.waitForTimeout(600);
await pickTrace('mmode', 'lv');
await page.waitForFunction(() => { const m = window.echosim.echo.mmodeData; return !!(m && m.key === 'lv' && m.read); }, null, { timeout: 60000 });
const plaxLv = await page.evaluate(async () => { const { measureModel } = await import('./js/report.js'); const r = window.echosim.echo.mmodeData.read, m = measureModel(window.echosim.state.path); return { lvidd: r.lvidd, lvids: r.lvids, ivsd: r.ivsd, pwd: r.pwd, model: m.lvidd / 10 }; });
await shot('mmode_plax_lv.png');
await pickTrace('mmode', 'ao');
await page.waitForFunction(() => { const m = window.echosim.echo.mmodeData; return !!(m && m.key === 'ao' && m.read); }, null, { timeout: 60000 });
const plaxAo = await page.evaluate(() => { const r = window.echosim.echo.mmodeData.read; return { ao: r.ao, la: r.la_max }; });
await shot('mmode_plax_ao.png');
await page.click('[data-view="SC_IVC"]'); await page.waitForTimeout(600);
await pickTrace('mmode', 'ivc');
await page.waitForFunction(() => { const m = window.echosim.echo.mmodeData; return !!(m && m.key === 'ivc' && m.ivc); }, null, { timeout: 60000 });
const ivcM = await page.evaluate(() => window.echosim.echo.mmodeData.ivc);
await shot('mmode_ivc.png');
// the Sniff button narrows the live IVC and it recovers
const ivcR = () => page.evaluate(async () => { const { geometryAt } = await import('./js/cardiac-model.js'); const v = geometryAt(0, window.echosim.state.path).A.ivc; return v.r / v.rMax; });
await page.evaluate(() => { const b = document.getElementById('playBtn'); if (/Play/.test(b.textContent)) b.click(); });
await page.click('#sniffBtn');
let sniffMin = 1;
for (let i = 0; i < 30; i++) { sniffMin = Math.min(sniffMin, await ivcR()); await page.waitForTimeout(60); }
await shot('sniff_after.png');
await pickTrace('doppler');
const mres = { mmodeTapse: mm.exc, caliper: mmCal, tapseModel, tdi, ivc: ivcM, sniffMin, plaxLv, plaxAo };

const efNum = parseFloat(res.ef), distNum = parseFloat(res.dist);
const checks = {
  'no console errors': errors.length === 0,
  'caliper returned a distance': Number.isFinite(distNum) && distNum > 0,
  'Simpson returned a plausible EF': Number.isFinite(efNum) && efNum > 0 && efNum < 100,
  'accuracy compared against truth': res.err.includes('truth'),
  'M-mode caliper returns the traced TAPSE (within 1 mm)': Number.isFinite(mmCal) && Math.abs(mmCal - mm.exc) <= 0.1,
  'M-mode TAPSE is the model TAPSE (within 2 mm)': Math.abs(mm.exc - tapseModel) <= 0.2,
  'PLAX M-mode reads LVIDd automatically (within 3 mm of the report)': Number.isFinite(plaxLv.lvidd) && Math.abs(plaxLv.lvidd - plaxLv.model) <= 0.3 && plaxLv.lvids < plaxLv.lvidd && plaxLv.ivsd > 0.5 && plaxLv.pwd > 0.5,
  'PLAX M-mode reads the aortic root and LA automatically': plaxAo.ao > 2 && plaxAo.ao < 4 && plaxAo.la > 2.5 && plaxAo.la < 4.2,
  "tissue Doppler reads s', e', a' and E/e'": tdi.s > 3 && tdi.e > 3 && tdi.a > 1 && tdi.EoverE > 3 && tdi.EoverE < 15,
  'IVC M-mode shows the sniff collapse (normal: > 50 %)': ivcM.collapse > 0.5 && ivcM.dMax > 1.2 && ivcM.dMax <= 2.1,
  'Sniff button narrows the live IVC by half': sniffMin < 0.6,
};
const pass = Object.values(checks).every(Boolean);
for (const [name, ok] of Object.entries(checks)) console.log(`  ${ok ? 'PASS' : 'FAIL'} ${name}`);
console.log(JSON.stringify({ pass, result: res, trace: mres, errors }, null, 2));

await browser.close();
server.close();
process.exit(pass ? 0 : 1);
