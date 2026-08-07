// tools/shoot.mjs — deterministic screenshot capture of EchoSim.
// Usage: node tools/shoot.mjs <outDir>
// Serves the repo, drives the UI to fixed views/pathologies/phases and captures
// full-page + echo-cropped PNGs plus any runtime console/page errors.
import { chromium } from 'playwright';
import http from 'http';
import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const OUT = path.resolve(process.argv[2] || '.review/shots');
fs.mkdirSync(OUT, { recursive: true });

const MIME = { '.html': 'text/html', '.js': 'text/javascript', '.css': 'text/css', '.json': 'application/json' };
const server = http.createServer((req, res) => {
  let p = decodeURIComponent(req.url.split('?')[0]);
  if (p === '/') p = '/index.html';
  fs.readFile(path.join(ROOT, p), (e, d) => {
    if (e) { res.writeHead(404); res.end('nf'); return; }
    res.writeHead(200, { 'Content-Type': MIME[path.extname(p)] || 'application/octet-stream' });
    res.end(d);
  });
});
const PORT = 8091;
await new Promise((r) => server.listen(PORT, r));

const browser = await chromium.launch({
  executablePath: '/opt/pw-browsers/chromium-1194/chrome-linux/chrome',
  args: ['--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--no-sandbox'],
});
const page = await browser.newPage({ viewport: { width: 1440, height: 880 } });
const errors = [];
page.on('console', (m) => { if (m.type() === 'error' && !/favicon/.test(m.text())) errors.push('console: ' + m.text()); });
page.on('pageerror', (e) => errors.push('pageerror: ' + e.message));

await page.goto(`http://localhost:${PORT}/`, { waitUntil: 'networkidle' });
await page.waitForTimeout(1200);
// dismiss first-run help overlay if present
const got = await page.$('#helpGot');
if (got) { try { await got.click(); } catch (e) {} }
await page.waitForTimeout(300);

async function setPhase(frac) {
  await page.$eval('#scrub', (el, v) => {
    el.value = String(Math.round(v * 1000));
    el.dispatchEvent(new Event('input', { bubbles: true }));
  }, frac);
  await page.waitForTimeout(500);
}
async function setView(v) {
  await page.click(`[data-view="${v}"]`);
  await page.waitForTimeout(900);
}
async function setPath(v) {
  await page.selectOption('#pathology', v);
  await page.waitForTimeout(900);
}
async function shot(name) {
  await page.screenshot({ path: path.join(OUT, name + '.png') });
  // echo-only crop for close inspection
  const echo = await page.$('.viewEcho');
  if (echo) await echo.screenshot({ path: path.join(OUT, name + '_echo.png') });
}

// scenarios: [label, pathology, view, phaseFraction]
const scenes = [
  ['normal_plax_dia', 'normal', 'PLAX', 0.52],
  ['normal_psax_sys', 'normal', 'PSAX', 0.20],
  ['normal_a4c_dia', 'normal', 'A4C', 0.52],
  ['mr_a4c_sys', 'mr', 'A4C', 0.18],
  ['as_plax_sys', 'as', 'PLAX', 0.18],
  ['dcm_a4c_sys', 'dcm', 'A4C', 0.20],
  ['rwma_a4c_sys', 'rwma', 'A4C', 0.20],
  ['effusion_plax', 'effusion', 'PLAX', 0.30],
  ['vsd_a4c_sys', 'vsd', 'A4C', 0.18],
  ['asd_subcostal', 'asd', 'SUBCOSTAL', 0.30],
];

for (const [label, pth, view, ph] of scenes) {
  await setPath(pth);   // changing pathology auto-switches view
  await setView(view);  // enforce the intended view
  await setPhase(ph);
  await shot(label);
}

fs.writeFileSync(path.join(OUT, 'errors.json'), JSON.stringify(errors, null, 2));
console.log(JSON.stringify({ out: OUT, images: scenes.length, errors }, null, 2));
await browser.close();
server.close();
