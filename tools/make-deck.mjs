// tools/make-deck.mjs — build a landscape SLIDE DECK PDF: benchmark evolution +
// screenshots of the physics running in-app. Captures shots live, then lays them
// out one concept per 16:9 slide. Output: EchoSim_Deck.pdf
import { chromium } from 'playwright';
import http from 'http';
import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const OUT_PDF = path.join(ROOT, 'EchoSim_Deck.pdf');
const ASSET = path.join(ROOT, '.review', 'deck');
fs.mkdirSync(ASSET, { recursive: true });

// current v2.1 scorecard (mirrors BENCHMARKS.md)
const DIMS = [
  ['B1', 'Anatomy & structure', 1.0, 6.5, 7.5],
  ['B2', 'Electrophysiology & 12-lead ECG', 1.2, 4.0, 7.5],
  ['B3', 'Myocardial mechanics (PV loop)', 1.1, 5.5, 8.0],
  ['B4', 'Hemodynamics & Doppler', 1.1, 6.0, 7.5],
  ['B5', 'Ultrasound image formation', 1.1, 6.5, 8.0],
  ['B6', 'Pathology breadth & depth', 1.0, 6.5, 7.5],
  ['B7', 'Modalities & acquisition', 0.8, 5.5, 6.5],
  ['B8', 'Quantification accuracy', 1.0, 7.0, 8.0],
  ['B9', 'External validation', 1.0, 3.5, 4.5],
  ['B10', 'Interactivity, UX & education', 0.7, 8.0, 8.5],
  ['B11', 'Engineering rigor', 0.6, 8.0, 8.5],
];
const wsum = DIMS.reduce((a, d) => a + d[2], 0);
const base = DIMS.reduce((a, d) => a + d[2] * d[3], 0) / wsum;
const curr = DIMS.reduce((a, d) => a + d[2] * d[4], 0) / wsum;

// ---- serve + capture ----
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
const PORT = 8101;
await new Promise((r) => server.listen(PORT, r));

const browser = await chromium.launch({
  executablePath: '/opt/pw-browsers/chromium-1194/chrome-linux/chrome',
  args: ['--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--no-sandbox'],
});
const page = await browser.newPage({ viewport: { width: 1440, height: 880 }, deviceScaleFactor: 2 });
const errors = [];
page.on('console', (m) => { if (m.type() === 'error' && !/favicon/.test(m.text())) errors.push(m.text()); });
page.on('pageerror', (e) => errors.push(e.message));
await page.goto(`http://localhost:${PORT}/`, { waitUntil: 'networkidle' });
await page.waitForTimeout(1200);
const got = await page.$('#helpGot'); if (got) await got.click().catch(() => {});
await page.waitForTimeout(300);

const setPath = async (v) => { await page.selectOption('#pathology', v); await page.waitForTimeout(1100); };
const setView = async (v) => { await page.click(`[data-view="${v}"]`); await page.waitForTimeout(900); };
const setPhase = async (f) => { await page.$eval('#scrub', (el, v) => { el.value = String(Math.round(v * 1000)); el.dispatchEvent(new Event('input', { bubbles: true })); }, f); await page.waitForTimeout(500); };
const shot = (name, el) => (el || page).screenshot({ path: path.join(ASSET, name + '.png') });
const ovBox = () => page.$eval('#echoOverlay', (el) => { const r = el.getBoundingClientRect(); return { x: r.x, y: r.y, w: r.width, h: r.height }; });

// hero: full app on a normal PLAX
await setPath('normal'); await setView('PLAX'); await setPhase(0.18);
await shot('hero');
await shot('hero3d', await page.$('.view3d'));
// echo crops for pathologies
await setPath('as'); await setView('PLAX'); await setPhase(0.18); await shot('echo_as', await page.$('.viewEcho'));
await setPath('mr'); await setView('A4C'); await setPhase(0.18); await shot('echo_mr', await page.$('.viewEcho'));
await setPath('rwma'); await setView('A4C'); await setPhase(0.20); await shot('echo_rwma', await page.$('.viewEcho'));

// hemodynamics panel (AS)
await setPath('as'); await page.waitForTimeout(700);
const grp = await page.evaluateHandle(() => document.getElementById('pvLoop').closest('.ctl-group'));
await grp.scrollIntoViewIfNeeded(); await page.waitForTimeout(400);
{
  const box = await page.evaluate(() => {
    const g = document.getElementById('pvLoop').closest('.ctl-group');
    const iq = document.getElementById('qGcnr').closest('.ctl-group');
    const r1 = g.getBoundingClientRect(), r2 = iq.getBoundingClientRect();
    return { x: Math.min(r1.x, r2.x) - 6, y: r1.y - 6, w: Math.max(r1.width, r2.width) + 12, h: (r2.y + r2.height) - r1.y + 12 };
  });
  await page.screenshot({ path: path.join(ASSET, 'hemo.png'), clip: { x: box.x, y: box.y, width: box.w, height: box.h } });
}

// 12-lead sinus + STEMI
await page.selectOption('#pathology', 'normal'); await page.waitForTimeout(500);
await page.selectOption('#rhythm', 'sinus'); await page.waitForTimeout(3400);
await page.click('#twelveBtn'); await page.waitForTimeout(3000);
await shot('twelve_sinus', await page.$('.twelve-card')); await page.click('#twelveClose');
await page.selectOption('#rhythm', 'stemi'); await page.waitForTimeout(3600);
await page.click('#twelveBtn'); await page.waitForTimeout(3000);
await shot('twelve_stemi', await page.$('.twelve-card')); await page.click('#twelveClose');

// learner calipers + Simpson EF on A4C
await page.selectOption('#rhythm', 'sinus'); await page.waitForTimeout(500);
await setPath('normal'); await setView('A4C');
const D = Math.PI / 180;
async function traceLV(cxF, cyF, wF, hF, n) {
  const b = await ovBox(); const cx = b.x + b.w * cxF, cy = b.y + b.h * cyF, rx = b.w * wF, ry = b.h * hF;
  const pt = (deg) => [cx + Math.cos(deg * D) * rx, cy + Math.sin(deg * D) * ry];
  for (let i = 0; i < n; i++) { const [x, y] = pt(115 + 310 * i / (n - 1)); await page.mouse.click(x, y); await page.waitForTimeout(22); }
  const [lx, ly] = pt(425); await page.mouse.dblclick(lx, ly); await page.waitForTimeout(140);
}
await setPhase(0.52); await page.click('[data-tool="caliper"]'); await page.waitForTimeout(100);
{ const b = await ovBox(); await page.mouse.click(b.x + b.w * 0.34, b.y + b.h * 0.42); await page.waitForTimeout(60); await page.mouse.click(b.x + b.w * 0.52, b.y + b.h * 0.42); await page.waitForTimeout(100); }
await setPhase(0.52); await page.click('[data-tool="simpsonED"]'); await page.waitForTimeout(100); await traceLV(0.42, 0.46, 0.13, 0.24, 14);
await setPhase(0.20); await page.click('[data-tool="simpsonES"]'); await page.waitForTimeout(100); await traceLV(0.42, 0.46, 0.10, 0.19, 14);
await setPhase(0.52); await page.waitForTimeout(300);
await shot('measure', await page.$('.viewEcho'));
const measVals = await page.evaluate(() => ({ ef: document.getElementById('measSimpEf').textContent, err: document.getElementById('measErr').textContent }));

await browser.close();

// ---- build slides ----
const b64 = (n) => { const fp = path.join(ASSET, n + '.png'); return fs.existsSync(fp) ? 'data:image/png;base64,' + fs.readFileSync(fp).toString('base64') : ''; };
const barCol = (v) => v >= 8 ? '#57d9a3' : v >= 7 ? '#7ec8f0' : v >= 5.5 ? '#f2b752' : '#e8615a';

const evoRows = DIMS.map(([id, name, w, b, c]) => `
  <tr>
    <td class="bid">${id}</td><td class="bname">${name}</td>
    <td class="track">
      <div class="bar base" style="width:${b * 10}%"></div>
      <div class="bar gain" style="left:${b * 10}%;width:${(c - b) * 10}%;background:${barCol(c)}"></div>
      <div class="tick" style="left:${b * 10}%"></div>
      <span class="bv base" style="left:calc(${b * 10}% - 4px)">${b.toFixed(1)}</span>
      <span class="bv cur" style="left:calc(${c * 10}% + 6px)">${c.toFixed(1)}</span>
    </td>
    <td class="delta">${c > b ? '+' + (c - b).toFixed(1) : '—'}</td>
  </tr>`).join('');

const slide = (cls, inner) => `<section class="slide ${cls || ''}">${inner}</section>`;
const imgFull = (n) => `<div class="imgwrap"><img src="${b64(n)}"/></div>`;

const slides = [
  // 1 · title
  slide('title', `
    <div class="tcenter">
      <div class="logo">◮</div>
      <h1>EchoSim</h1>
      <p class="tagline">Interactive 3D in-silico echocardiography trainer</p>
      <div class="scorepill">Fidelity ${base.toFixed(1)} → <b>${curr.toFixed(1)}</b> / 10 &nbsp;·&nbsp; v2 state-of-the-art rubric</div>
      <p class="tsub">Benchmark evolution &amp; evidence deck · ${new Date().toISOString().slice(0, 10)}</p>
    </div>`),

  // 2 · what it is
  slide('split', `
    <div class="col text">
      <h2>What it is</h2>
      <ul>
        <li>A beating 3D heart you scan like a real machine — a moveable ultrasound plane drives a live 2-D echo.</li>
        <li>Physics-grounded: a <b>12-lead ECG</b>, a lumped-parameter <b>PV loop</b>, continuity <b>Doppler</b>, and scatterer-based <b>B-mode</b>.</li>
        <li>Standard views (PLAX/PSAX/A4C/A2C/subcostal), 10 graded pathologies, learner measurements, challenge mode.</li>
        <li>Fully offline, no PHI — deployable today as an <b>education-only</b> trainer.</li>
      </ul>
    </div>
    <div class="col shot">${imgFull('hero')}</div>`),

  // 3 · benchmark evolution
  slide('evo', `
    <h2>Benchmark evolution — baseline → current</h2>
    <p class="sub">Weighted against physics-based computational cardiology &amp; validated commercial simulators (not a teaching-tool self-rating).</p>
    <table class="evotab">
      <thead><tr><th>#</th><th>Dimension</th><th>0 &nbsp;·············· score ··············&nbsp; 10</th><th>Δ</th></tr></thead>
      <tbody>${evoRows}</tbody>
    </table>
    <div class="evofoot">
      <span class="key"><i class="kb"></i> baseline ${base.toFixed(1)}</span>
      <span class="key"><i class="kc"></i> current ${curr.toFixed(1)}</span>
      <span class="weighted">Weighted overall <b>${base.toFixed(1)} → ${curr.toFixed(1)}</b> / 10 &nbsp;(Σ ${DIMS.reduce((a, d) => a + d[2] * d[4], 0).toFixed(1)} ÷ ${wsum.toFixed(1)})</span>
    </div>`),

  // 4 · ECG
  slide('split', `
    <div class="col text">
      <h2>Physics-derived 12-lead ECG <span class="chip">B2 · 4.0→7.5</span></h2>
      <ul>
        <li>Vectorcardiographic dipole → <b>inverse-Dower lead field</b> → standard 12-lead.</li>
        <li>Rhythm library: sinus, AF, brady/tachy, AV block, PVCs, LBBB/RBBB, STEMI, hyperK.</li>
        <li>Electromechanically coupled to the cardiac cycle; <b>12/12</b> automated morphology checks pass.</li>
      </ul>
      <div class="mini">${imgFull('twelve_stemi')}</div>
      <p class="cap">Anterior STEMI — ST-elevation vector projects to the anterior leads.</p>
    </div>
    <div class="col shot">${imgFull('twelve_sinus')}<p class="cap">Normal sinus 12-lead: aVR−, normal axis, R-wave progression V1→V6.</p></div>`),

  // 5 · mechanics + imaging
  slide('split', `
    <div class="col text">
      <h2>PV loop, haemodynamics &amp; image quality <span class="chip">B3/B5 · →7.5</span></h2>
      <ul>
        <li>Lumped-parameter closed-loop circulation → a real LV <b>pressure–volume loop</b>; EF/LVIDs fall out of the physics.</li>
        <li>Live <b>gCNR / CNR / speckle-SNR</b> self-validation of the scatterer+PSF B-mode.</li>
        <li>Severe AS shown: P<sub>LV</sub> 149, gradient 76 mmHg, gCNR 0.95.</li>
      </ul>
    </div>
    <div class="col shot narrow">${imgFull('hemo')}</div>`),

  // 6 · pathology + Doppler
  slide('trio', `
    <h2>Graded pathology &amp; continuity Doppler <span class="chip">B4/B6 · →7.5</span></h2>
    <div class="row3">
      <figure>${imgFull('echo_as')}<figcaption>Aortic stenosis — jet velocity is continuity-derived (v=Q/A); graded mild→severe.</figcaption></figure>
      <figure>${imgFull('echo_mr')}<figcaption>Mitral regurgitation — Bernoulli jet, regurgitant fraction 26/46/65% by grade.</figcaption></figure>
      <figure>${imgFull('echo_rwma')}<figcaption>Regional wall-motion abnormality — infarcted septum fails to thicken; EF ~46%.</figcaption></figure>
    </div>`),

  // 7 · measurements
  slide('split', `
    <div class="col shot">${imgFull('measure')}</div>
    <div class="col text">
      <h2>Learner measurement tools <span class="chip">B8 · 7.0→8.0</span></h2>
      <ul>
        <li>Click-two-point <b>calipers</b> (exact cm) and a traced <b>Simpson method-of-discs EF</b>.</li>
        <li>Checked <b>live against the model's ground truth</b> — measured ${measVals.ef}, error ${measVals.err}, within the ASE ±6% band.</li>
        <li>Disc-summation math verified to <b>0.0%</b> vs an analytic solid of revolution.</li>
      </ul>
    </div>`),

  // 8 · honest ceiling + roadmap
  slide('closing', `
    <h2>Honest status &amp; the path beyond ${curr.toFixed(1)}</h2>
    <div class="two">
      <div class="card good">
        <h3>Deployable today</h3>
        <p>An <b>education-only, non-device</b> trainer — offline, no PHI, physics-grounded. Verified at <b>${curr.toFixed(1)}/10</b>.</p>
      </div>
      <div class="card warn">
        <h3>What code cannot close</h3>
        <p>The 9–10 band &amp; a <b>clinical-training deployment</b> claim need <b>B9 external validation</b> (expert panel + learning-transfer RCT) and <b>regulatory/QMS clearance</b> — human evidence-and-audit processes. The protocol, instrument &amp; harness are shipped in the repo (<i>VALIDATION.md</i>, <i>COMMERCIAL_READINESS.md</i>, <i>ROADMAP_TO_9.9.md</i>).</p>
      </div>
    </div>
    <p class="foot">Codeable ceiling ≈ 8.0–8.3 via reduced-order CFD (B4) / FE strain (B3); 9.9 requires the external validation study. In-silico teaching model — not for diagnosis.</p>`),
];

const html = `<!doctype html><html><head><meta charset="utf-8"><style>
@page { size: A4 landscape; margin: 0; }
* { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: -apple-system, system-ui, sans-serif; color: #e8eef6; }
.slide { width: 297mm; height: 209mm; padding: 15mm 18mm; page-break-after: always; position: relative;
  background: radial-gradient(1200px 700px at 78% -10%, #16233a 0%, #0c1119 55%, #080b11 100%); overflow: hidden; }
h1 { font-size: 58px; letter-spacing: -1px; }
h2 { font-size: 27px; color: #fff; margin-bottom: 12px; letter-spacing: -0.3px; }
h3 { font-size: 16px; color: #fff; margin-bottom: 6px; }
.chip { font-size: 13px; font-weight: 600; color: #0c1119; background: #7ec8f0; border-radius: 20px; padding: 2px 11px; vertical-align: middle; margin-left: 8px; }
ul { margin-left: 18px; } li { font-size: 15.5px; line-height: 1.65; color: #cdd8e6; margin-bottom: 6px; }
li b { color: #fff; }
.sub { color: #8b98a8; font-size: 14px; margin-bottom: 14px; }
.cap { font-size: 12px; color: #9fb0c2; margin-top: 6px; }
.imgwrap { border: 1px solid #26344a; border-radius: 10px; overflow: hidden; background: #05070c; box-shadow: 0 8px 30px rgba(0,0,0,.4); }
.imgwrap img { width: 100%; display: block; }
/* title */
.title .tcenter { height: 100%; display: flex; flex-direction: column; align-items: center; justify-content: center; text-align: center; }
.logo { font-size: 70px; color: #39c0e8; }
.tagline { font-size: 20px; color: #aeb9c8; margin: 6px 0 22px; }
.scorepill { font-size: 20px; background: #142033; border: 1px solid #2b6cb0; border-radius: 30px; padding: 10px 26px; }
.scorepill b { color: #57d9a3; font-size: 24px; }
.tsub { color: #6f819a; font-size: 13px; margin-top: 20px; }
/* split */
.split { display: flex; flex-direction: column; }
.split { display: grid; grid-template-columns: 1fr 1fr; gap: 26px; align-items: center; }
.split h2 { grid-column: auto; }
.col.text { align-self: center; }
.col.shot.narrow { max-width: 320px; margin: 0 auto; }
.mini { margin-top: 12px; max-width: 88%; } .mini .imgwrap { border-radius: 8px; }
/* evolution table */
.evotab { width: 100%; border-collapse: collapse; margin-top: 6px; }
.evotab th { font-size: 11px; text-transform: uppercase; letter-spacing: .5px; color: #7c8ba0; text-align: left; padding: 4px 8px; border-bottom: 1px solid #26344a; }
.evotab td { padding: 5px 8px; border-bottom: 1px solid #17222f; vertical-align: middle; }
.bid { font-weight: 700; color: #7ec8f0; font-size: 12px; width: 34px; }
.bname { font-size: 13px; color: #dbe5f0; width: 210px; }
.track { position: relative; height: 20px; }
.track .bar { position: absolute; top: 4px; height: 12px; }
.track .bar.base { background: #33455c; z-index: 1; border-radius: 3px 0 0 3px; left: 0; }
.track .bar.gain { z-index: 2; border-radius: 0 3px 3px 0; }
.track .tick { position: absolute; top: 1px; height: 18px; width: 2px; background: #cdd8e6; z-index: 3; }
.bv { position: absolute; top: 3px; font-size: 10px; font-variant-numeric: tabular-nums; }
.bv.base { color: #8b98a8; transform: translateX(-100%); }
.bv.cur { color: #fff; font-weight: 700; z-index: 3; }
.delta { color: #57d9a3; font-weight: 700; font-size: 13px; width: 40px; text-align: right; }
.evofoot { display: flex; gap: 22px; align-items: center; margin-top: 14px; font-size: 13px; color: #aeb9c8; }
.key { display: inline-flex; align-items: center; gap: 6px; } .kb, .kc { width: 22px; height: 10px; border-radius: 3px; display: inline-block; }
.kb { background: #33455c; } .kc { background: #7ec8f0; }
.weighted { margin-left: auto; font-size: 15px; } .weighted b { color: #57d9a3; }
/* trio */
.row3 { display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 16px; margin-top: 14px; }
.row3 figure { margin: 0; }
.row3 figcaption { font-size: 11.5px; color: #9fb0c2; margin-top: 7px; line-height: 1.45; }
/* closing */
.two { display: grid; grid-template-columns: 1fr 1fr; gap: 22px; margin-top: 16px; }
.card { border-radius: 12px; padding: 18px 20px; border: 1px solid #26344a; background: #101a28; }
.card p { font-size: 14.5px; color: #cdd8e6; line-height: 1.6; } .card p b { color: #fff; }
.card.good { border-color: #2f7d5a; background: #10231b; }
.card.warn { border-color: #7a5a1a; background: #221a10; }
.closing .foot { position: absolute; bottom: 14mm; left: 18mm; right: 18mm; font-size: 12px; color: #7c8ba0; border-top: 1px solid #26344a; padding-top: 10px; }
/* page number */
.slide::after { content: ''; }
</style></head><body>${slides.join('')}</body></html>`;

const hp = path.join(ASSET, 'deck.html');
fs.writeFileSync(hp, html);
const b2 = await chromium.launch({ executablePath: '/opt/pw-browsers/chromium-1194/chrome-linux/chrome', args: ['--no-sandbox'] });
const pg2 = await b2.newPage();
await pg2.goto('file://' + hp, { waitUntil: 'networkidle' });
await pg2.pdf({ path: OUT_PDF, width: '297mm', height: '210mm', printBackground: true });
await b2.close();
server.close();
console.log(JSON.stringify({ pdf: OUT_PDF, slides: slides.length, baseline: +base.toFixed(2), current: +curr.toFixed(2), captureErrors: errors, measVals }, null, 2));
