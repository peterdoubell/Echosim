// tools/shoot-hemo.mjs — capture the Haemodynamics (PV loop) + Image-quality panels
// and dump their live DOM values, for a set of pathologies. Verifies the UI wiring.
import { chromium } from 'playwright';
import http from 'http';
import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const OUT = path.resolve(process.argv[2] || '.review/hemo');
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
const PORT = 8092;
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
const got = await page.$('#helpGot');
if (got) { try { await got.click(); } catch (e) {} }
await page.waitForTimeout(300);

const readVals = () => page.evaluate(() => {
  const t = (id) => (document.getElementById(id) || {}).textContent || '—';
  return {
    EDV: t('hEdv'), ESV: t('hEsv'), SV: t('hSv'), EF: t('hEf'),
    Plv: t('hPlv'), Pao: t('hPao'), grad: t('hGrad'), RF: t('hRf'),
    gCNR: t('qGcnr'), CNR: t('qCnr'), SNR: t('qSnr'),
  };
});

const cases = ['normal', 'as', 'mr', 'dcm', 'rwma', 'ms'];
const table = {};
for (const c of cases) {
  await page.selectOption('#pathology', c);
  await page.waitForTimeout(1000);
  table[c] = await readVals();
}

// screenshot the console scrolled to the Haemodynamics panel (normal case)
await page.selectOption('#pathology', 'as');
await page.waitForTimeout(1000);
const pv = await page.$('#pvLoop');
if (pv) {
  const grp = await page.evaluateHandle(() => document.getElementById('pvLoop').closest('.ctl-group'));
  await grp.scrollIntoViewIfNeeded();
  await page.waitForTimeout(400);
  // capture the two new panels together
  const box = await page.evaluate(() => {
    const g = document.getElementById('pvLoop').closest('.ctl-group');
    const iq = document.getElementById('qGcnr').closest('.ctl-group');
    const r1 = g.getBoundingClientRect(), r2 = iq.getBoundingClientRect();
    return { x: Math.min(r1.x, r2.x), y: r1.y, w: Math.max(r1.width, r2.width), h: (r2.y + r2.height) - r1.y };
  });
  await page.screenshot({ path: path.join(OUT, 'hemo_panels_as.png'),
    clip: { x: box.x, y: box.y, width: box.w, height: box.h } });
}

fs.writeFileSync(path.join(OUT, 'values.json'), JSON.stringify(table, null, 2));
console.log(JSON.stringify({ out: OUT, errors, table }, null, 2));
await browser.close();
server.close();
