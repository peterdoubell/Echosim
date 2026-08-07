// tools/shoot-ecg.mjs — screenshot the ECG strip + 12-lead modal (normal & AF).
import { chromium } from 'playwright';
import http from 'http'; import fs from 'fs'; import path from 'path';
import { fileURLToPath } from 'url';
const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const OUT = path.join(ROOT, '.review/ecg'); fs.mkdirSync(OUT, { recursive: true });
const MIME = { '.html': 'text/html', '.js': 'text/javascript', '.css': 'text/css', '.json': 'application/json' };
const srv = http.createServer((rq, rs) => { let p = decodeURIComponent(rq.url.split('?')[0]); if (p === '/') p = '/index.html'; fs.readFile(path.join(ROOT, p), (e, d) => { if (e) { rs.writeHead(404); rs.end(); return; } rs.writeHead(200, { 'Content-Type': MIME[path.extname(p)] || 'application/octet-stream' }); rs.end(d); }); });
await new Promise((r) => srv.listen(8123, r));
const b = await chromium.launch({ executablePath: '/opt/pw-browsers/chromium-1194/chrome-linux/chrome', args: ['--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--no-sandbox'] });
const pg = await b.newPage({ viewport: { width: 1440, height: 880 } });
const errs = []; pg.on('pageerror', (e) => errs.push(e.message)); pg.on('console', (m) => { if (m.type() === 'error' && !/favicon/.test(m.text())) errs.push('console: ' + m.text()); });
await pg.goto('http://localhost:8123/', { waitUntil: 'networkidle' });
const got = await pg.$('#helpGot'); if (got) await got.click().catch(() => {});
await pg.waitForTimeout(3500);                    // let the ECG buffer fill
await pg.locator('.cycle-readout').screenshot({ path: path.join(OUT, 'strip_sinus.png') });
await pg.click('#twelveBtn'); await pg.waitForTimeout(3200);
await pg.locator('.twelve-card').screenshot({ path: path.join(OUT, 'twelve_sinus.png') });
await pg.click('#twelveClose');
// AF
await pg.selectOption('#rhythm', 'afib'); await pg.waitForTimeout(3800);
await pg.click('#twelveBtn'); await pg.waitForTimeout(3200);
await pg.locator('.twelve-card').screenshot({ path: path.join(OUT, 'twelve_afib.png') });
await pg.click('#twelveClose');
// STEMI
await pg.selectOption('#rhythm', 'stemi'); await pg.waitForTimeout(3800);
await pg.click('#twelveBtn'); await pg.waitForTimeout(3000);
await pg.locator('.twelve-card').screenshot({ path: path.join(OUT, 'twelve_stemi.png') });
console.log(JSON.stringify({ errors: errs }, null, 2));
await b.close(); srv.close();
