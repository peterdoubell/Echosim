// tools/review-pack.mjs — build the image pack for an expert-panel review.
//
// Produces, under <outDir>:
//   tissue/<scene>/<VIEW>.png  tissue-label maps of each imaging plane (classify()
//                              over the sector; legend in README.md), in the
//                              echo display's orientation: depth down, probe.lat
//                              to the right, 0.66 rad half-angle sector
//   echo/<scene>/<VIEW>.png    the rendered Echo Display panel (B-mode, colour,
//                              spectral trace) for each view at a fixed phase
//   3d/<VIEW>.png              the 3D anatomy panel with the imaging plane
//   audit_*.txt, hemodynamics.txt
// Scenes match the pack used for reviews 8 and 9 (normal ED/ES and TEE, then
// each pathology at the phase where its lesion shows).
//
// Usage: node tools/review-pack.mjs <outDir>
//        node tools/review-pack.mjs <outDir> --scene <pathology> <VIEW,VIEW,...> <phase>
//          renders just that scene (echo + tissue) into <outDir>/echo/<pathology>_<phase>/ and
//          <outDir>/tissue/<pathology>_<phase>/, for checking a change quickly. Env PORT picks the server port.
import fs from 'fs';
import path from 'path';
import http from 'http';
import zlib from 'zlib';
import { execFileSync } from 'child_process';
import { fileURLToPath } from 'url';
import { geometryAt, classify, TISSUE } from '../js/cardiac-model.js';
import { TTE_VIEWS, EXTRA_VIEWS, TEE_VIEWS } from '../js/views.js';
import { hemoSummary } from '../js/hemodynamics.js';

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const OUT = path.resolve(process.argv[2] || 'review_pack');
const ONLY = process.env.ONLY || '';           // 'tissue', 'echo', '3d' or 'audit' to build one part
fs.mkdirSync(OUT, { recursive: true });

const PATHS = {
  normal: {}, dcm: { dilated: true }, as: { aorticStenosis: true, lvh: true }, mr: { mr: true },
  ms: { mitralStenosis: true }, phtn: { rvpo: true, tr: true }, effusion: { effusion: true },
  asd: { asd: true }, vsd: { vsd: true }, rwma: { rwma: 'septal' }, tr: { tr: true },
};
const TTE_ALL = 'PLAX,PSAX,PSAX_AV,PSAX_MV,A4C,A5C,A2C,A3C,SUBCOSTAL,RVIT,SC_IVC,SSN';
const TEE_ALL = 'ME4C,ME2C,MELAX,MEAVSAX,MEBICAVAL,MERVIO,MELAA,TGSAX,DESCAO';
// [scene dir, pathology key, views, phase]
let ECHO_SCENES = [
  ['normal_ED', 'normal', TTE_ALL, 0.99], ['normal_ES', 'normal', 'PLAX,PSAX,PSAX_AV,PSAX_MV,A4C,A2C,SUBCOSTAL', 0.42],
  ['normal_TEE', 'normal', TEE_ALL, 0.99],
  ['dcm', 'dcm', 'A4C,PSAX,PLAX', 0.99], ['as', 'as', 'PLAX,PSAX_AV', 0.2], ['mr', 'mr', 'A4C,PLAX', 0.25],
  ['ms', 'ms', 'PSAX_MV,A4C,PLAX', 0.62], ['phtn', 'phtn', 'PSAX,A4C', 0.3], ['eff', 'effusion', 'PLAX,A4C,SUBCOSTAL', 0.99],
  ['asd', 'asd', 'SUBCOSTAL,A4C,MEBICAVAL', 0.5], ['vsd', 'vsd', 'A4C,PSAX', 0.25], ['rwma', 'rwma', 'PSAX,A4C', 0.42],
  ['tr', 'tr', 'A4C,RVIT', 0.25],
];
let TISSUE_SCENES = [
  ['normal_ED', 'normal', `${TTE_ALL},${TEE_ALL}`, 0.99], ['normal_ES', 'normal', `${TTE_ALL},${TEE_ALL}`, 0.42],
  ...['dcm', 'as', 'mr', 'ms', 'phtn', 'effusion', 'asd', 'vsd', 'rwma', 'tr'].flatMap((p) => [
    [`${p}_ED`, p, 'PLAX,PSAX,PSAX_AV,PSAX_MV,A4C,A2C,SUBCOSTAL,ME4C,MEBICAVAL', 0.99],
    [`${p}_ES`, p, 'PLAX,PSAX,A4C,SUBCOSTAL', 0.42],
  ]),
];

// quick single-scene mode
let QUICK = false;
const si = process.argv.indexOf('--scene');
if (si > 0) {
  const [key, views, ph] = process.argv.slice(si + 1, si + 4);
  if (!PATHS[key] || !views || !(+ph >= 0)) { console.error('usage: --scene <pathology> <VIEW,...> <phase>'); process.exit(2); }
  const dir = `${key}_${ph}`;
  ECHO_SCENES = [[dir, key, views, +ph]];
  TISSUE_SCENES = [[dir, key, views, +ph]];
  QUICK = true;
}

// ---- tissue maps -------------------------------------------------------------
const COLOURS = {
  [TISSUE.OUTSIDE]: [0, 0, 0], [TISSUE.MYO]: [224, 128, 128], [TISSUE.LV]: [127, 0, 0], [TISSUE.RV]: [30, 60, 170],
  [TISSUE.LA]: [255, 150, 150], [TISSUE.RA]: [120, 170, 255], [TISSUE.AORTA]: [235, 110, 0], [TISSUE.VALVE]: [255, 240, 0],
  [TISSUE.PERICARDIUM]: [0, 230, 230], [TISSUE.LUNG]: [70, 70, 70], [TISSUE.PERI_LINE]: [255, 255, 255],
  [TISSUE.LIVER]: [139, 90, 43], [TISSUE.VEIN]: [110, 0, 230], [TISSUE.FAT]: [230, 200, 150], [TISSUE.VWALL]: [170, 170, 170],
};
const crcTable = new Int32Array(256).map((_, n) => { let c = n; for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1; return c; });
const crc32 = (buf) => { let c = -1; for (const b of buf) c = crcTable[(c ^ b) & 255] ^ (c >>> 8); return (c ^ -1) >>> 0; };
function png(w, h, rgb) {
  const raw = Buffer.alloc((w * 3 + 1) * h);
  for (let y = 0; y < h; y++) { raw[y * (w * 3 + 1)] = 0; rgb.copy(raw, y * (w * 3 + 1) + 1, y * w * 3, (y + 1) * w * 3); }
  const chunk = (type, data) => {
    const len = Buffer.alloc(4); len.writeUInt32BE(data.length);
    const td = Buffer.concat([Buffer.from(type), data]);
    const crc = Buffer.alloc(4); crc.writeUInt32BE(crc32(td));
    return Buffer.concat([len, td, crc]);
  };
  const ihdr = Buffer.alloc(13); ihdr.writeUInt32BE(w, 0); ihdr.writeUInt32BE(h, 4); ihdr[8] = 8; ihdr[9] = 2;
  return Buffer.concat([Buffer.from([137, 80, 78, 71, 13, 10, 26, 10]), chunk('IHDR', ihdr), chunk('IDAT', zlib.deflateSync(raw)), chunk('IEND', Buffer.alloc(0))]);
}
const viewOf = (name) => TTE_VIEWS[name] || EXTRA_VIEWS[name] || TEE_VIEWS[name];
function tissueMap(name, phase, flags, size = 420) {
  const v = viewOf(name), G = geometryAt(phase, flags);
  const p = v.probe();
  let pos = p.pos;
  if (v.track) { const o = v.track(G.A); pos = pos.map((x, i) => x + o[i]); }
  const depth = v.depth || 16, half = 0.66, W = size, H = size;
  const sx = (2 * depth * Math.sin(half)) / W, sy = depth / H;
  const rgb = Buffer.alloc(W * H * 3);
  for (let j = 0; j < H; j++) for (let i = 0; i < W; i++) {
    const lat = (i - W / 2) * sx, along = j * sy;
    const o = (j * W + i) * 3;
    if (Math.atan2(Math.abs(lat), along) > half || Math.hypot(lat, along) > depth) { rgb[o] = rgb[o + 1] = rgb[o + 2] = 20; continue; }
    const x = pos[0] + p.dir[0] * along + p.lat[0] * lat, y = pos[1] + p.dir[1] * along + p.lat[1] * lat, z = pos[2] + p.dir[2] * along + p.lat[2] * lat;
    const c = COLOURS[classify(x, y, z, G, flags).tissue] || [255, 0, 255];
    rgb[o] = c[0]; rgb[o + 1] = c[1]; rgb[o + 2] = c[2];
  }
  return png(W, H, rgb);
}
if (!ONLY || ONLY === 'tissue') {
  for (const [dir, key, views, phase] of TISSUE_SCENES) {
    const d = path.join(OUT, 'tissue', dir); fs.mkdirSync(d, { recursive: true });
    for (const v of views.split(',')) fs.writeFileSync(path.join(d, `${v}.png`), tissueMap(v, phase, { ...PATHS[key], grade: 'severe' }));
  }
  console.log('tissue maps done');
}

// ---- echo and 3D screenshots ---------------------------------------------------------
if (!ONLY || ONLY === 'echo' || ONLY === '3d') {
  const { chromium } = await import('playwright');
  const MIME = { '.html': 'text/html', '.js': 'text/javascript', '.css': 'text/css', '.svg': 'image/svg+xml', '.png': 'image/png', '.webmanifest': 'application/manifest+json' };
  const PORT = +(process.env.PORT || 8140);
  const server = http.createServer((req, res) => {
    let p = decodeURIComponent(req.url.split('?')[0]); if (p === '/') p = '/index.html';
    fs.readFile(path.join(ROOT, p), (e, d) => { if (e) { res.writeHead(404); res.end(); return; } res.writeHead(200, { 'Content-Type': MIME[path.extname(p)] || 'application/octet-stream' }); res.end(d); });
  });
  await new Promise((r) => server.listen(PORT, r));
  const browser = await chromium.launch({
    executablePath: process.env.CHROME_PATH || ['/opt/pw-browsers/chromium', '/opt/pw-browsers/chromium-1194/chrome-linux/chrome'].find((p) => fs.existsSync(p)),
    args: ['--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--no-sandbox'],
  });
  const page = await browser.newPage({ viewport: { width: 1440, height: 900 } });
  const errors = [];
  page.on('console', (m) => { if (m.type() === 'error' && !/favicon/.test(m.text())) errors.push(m.text()); });
  page.on('pageerror', (e) => errors.push(e.message));
  await page.goto(`http://localhost:${PORT}/`, { waitUntil: 'networkidle' });
  await page.waitForTimeout(1500);
  const got = await page.$('#helpGot'); if (got) await got.click().catch(() => {});
  await page.evaluate(() => { const c = document.getElementById('labelsOn'); if (c.checked) c.click(); });
  const setScene = async (key, view, phase) => {
    await page.selectOption('#pathology', key);
    await page.waitForTimeout(500);
    await page.evaluate((v) => window.echosim.setView(v), view);
    await page.$eval('#scrub', (el, v) => { el.value = String(Math.round(v * 1000)); el.dispatchEvent(new Event('input', { bubbles: true })); }, phase);
    await page.waitForTimeout(1300);   // several frames: persistence settles
  };
  if (!ONLY || ONLY === 'echo') {
    for (const [dir, key, views, phase] of ECHO_SCENES) {
      const d = path.join(OUT, 'echo', dir); fs.mkdirSync(d, { recursive: true });
      for (const v of views.split(',')) { await setScene(key, v, phase); await (await page.$('.viewEcho')).screenshot({ path: path.join(d, `${v}.png`) }); }
      console.log('echo', dir);
    }
  }
  if (!QUICK && (!ONLY || ONLY === '3d')) {
    const d = path.join(OUT, '3d'); fs.mkdirSync(d, { recursive: true });
    for (const v of ['PLAX', 'A4C', 'PSAX', 'SUBCOSTAL']) { await setScene('normal', v, 0.99); await (await page.$('.panel')).screenshot({ path: path.join(d, `${v}.png`) }); }
    console.log('3d done');
  }
  fs.writeFileSync(path.join(OUT, 'errors.json'), JSON.stringify(errors, null, 2));
  await browser.close(); server.close();
}

// ---- audits ------------------------------------------------------------------------
if (!QUICK && (!ONLY || ONLY === 'audit')) {
  const run = (tool) => { try { return execFileSync('node', [path.join(ROOT, 'tools', tool)], { cwd: ROOT, encoding: 'utf8', maxBuffer: 1 << 26 }); } catch (e) { return (e.stdout || '') + (e.stderr || ''); } };
  fs.writeFileSync(path.join(OUT, 'audit_anatomy.txt'), run('verify-anatomy.mjs'));
  fs.writeFileSync(path.join(OUT, 'audit_atrial_phase.txt'), run('verify-atrial-phase.mjs'));
  fs.writeFileSync(path.join(OUT, 'audit_timing.txt'), run('verify-timing.mjs'));
  fs.writeFileSync(path.join(OUT, 'hemodynamics.txt'), Object.entries(PATHS).map(([k, p]) => `${k} ${JSON.stringify(p)} ${JSON.stringify(hemoSummary(p))}`).join('\n') + '\n');
  console.log('audits done');
}
fs.writeFileSync(path.join(OUT, 'tissue_legend.json'), JSON.stringify(Object.fromEntries(Object.entries(TISSUE).map(([k, v]) => [k, COLOURS[v]])), null, 1));
console.log('pack:', OUT);
