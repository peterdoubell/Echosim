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
const PORT = 8094;
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
  executablePath: process.env.CHROME_PATH || '/opt/pw-browsers/chromium-1194/chrome-linux/chrome',
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

const efNum = parseFloat(res.ef), distNum = parseFloat(res.dist);
const checks = {
  'no console errors': errors.length === 0,
  'caliper returned a distance': Number.isFinite(distNum) && distNum > 0,
  'Simpson returned a plausible EF': Number.isFinite(efNum) && efNum > 0 && efNum < 100,
  'accuracy compared against truth': res.err.includes('truth'),
};
const pass = Object.values(checks).every(Boolean);
for (const [name, ok] of Object.entries(checks)) console.log(`  ${ok ? 'PASS' : 'FAIL'} ${name}`);
console.log(JSON.stringify({ pass, result: res, errors }, null, 2));

await browser.close();
server.close();
process.exit(pass ? 0 : 1);
