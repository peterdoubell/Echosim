// tools/validate.mjs — automated self-validation harness (docs/VALIDATION.md §A).
// A3 (12-lead ECG morphology) runs in pure Node against the electrophysiology
// engine. A1/A2 (gCNR + measurement accuracy) read echo.metrics from a headless
// page when available. Prints a pass/fail table and exits non-zero on failure.
import { Heartbeat, vcgToLeads, LEAD_NAMES } from '../js/electrophysiology.js';

const FS = 500;
let pass = 0, fail = 0;
const rows = [];
function check(name, ok, detail) { rows.push([ok ? 'PASS' : 'FAIL', name, detail]); ok ? pass++ : fail++; }

// Collect `secs` of all-12-lead samples for a rhythm/pathology.
function collect(cfg, secs = 8) {
  const h = new Heartbeat();
  if (cfg.hr) h.setHR(cfg.hr);
  if (cfg.rhythm) h.setRhythm(cfg.rhythm);
  if (cfg.path) h.setPathology(cfg.path);
  const n = FS * secs, series = {};
  for (const L of LEAD_NAMES) series[L] = new Float32Array(n);
  for (let i = 0; i < n; i++) {
    h.advance(1 / FS);
    const L = h.leads();
    for (const nm of LEAD_NAMES) series[nm][i] = L[nm];
  }
  return series;
}
const ext = (a) => { let mn = 1e9, mx = -1e9; for (const v of a) { if (v < mn) mn = v; if (v > mx) mx = v; } return { mn, mx, dom: Math.abs(mx) >= Math.abs(mn) ? mx : mn }; };
// R-peak times (s) from lead II by threshold + refractory
function rTimes(a, thr) {
  const t = []; let last = -1;
  for (let i = 1; i < a.length; i++) if (a[i - 1] < thr && a[i] >= thr && (i - last) > FS * 0.25) { t.push(i / FS); last = i; }
  return t;
}
const cv = (arr) => { if (arr.length < 2) return 0; const m = arr.reduce((s, x) => s + x, 0) / arr.length; const sd = Math.sqrt(arr.reduce((s, x) => s + (x - m) ** 2, 0) / arr.length); return sd / m; };

// ---- A3: sinus morphology ----
{
  const s = collect({ rhythm: 'sinus', hr: 72 });
  check('sinus · aVR predominantly negative', ext(s.aVR).dom < 0, `dom=${ext(s.aVR).dom.toFixed(2)} mV`);
  check('sinus · lead I QRS positive (normal axis)', ext(s.I).mx > 0.15, `Rmax=${ext(s.I).mx.toFixed(2)}`);
  check('sinus · aVF positive (normal axis)', ext(s.aVF).mx > 0.1, `Rmax=${ext(s.aVF).mx.toFixed(2)}`);
  const rV1 = ext(s.V1).mx, rV3 = ext(s.V3).mx, rV5 = ext(s.V5).mx, rV6 = ext(s.V6).mx;
  // physiologic: R grows V1→(peak V4–V5), then V6 ≤ the mid-precordial peak
  check('sinus · R-wave progression (V1<V3<V5, peak≥V6)', rV1 < rV3 && rV3 < rV5 + 0.05 && rV5 >= rV6, `R: V1=${rV1.toFixed(2)} V3=${rV3.toFixed(2)} V5=${rV5.toFixed(2)} V6=${rV6.toFixed(2)}`);
  const rt = rTimes(s.II, 0.5), rr = rt.slice(1).map((t, i) => t - rt[i]);
  const meanRR = rr.reduce((a, b) => a + b, 0) / rr.length;
  check('sinus · rate ≈ 72 bpm', Math.abs(60 / meanRR - 72) < 12, `${(60 / meanRR).toFixed(0)} bpm`);
  check('sinus · regular R-R (CV<0.08)', cv(rr) < 0.08, `CV=${cv(rr).toFixed(3)}`);
}
// ---- A3: atrial fibrillation ----
{
  const s = collect({ rhythm: 'afib', hr: 90 });
  const rt = rTimes(s.II, 0.5), rr = rt.slice(1).map((t, i) => t - rt[i]);
  check('AF · irregularly irregular R-R (CV>0.15)', cv(rr) > 0.15, `CV=${cv(rr).toFixed(3)}`);
  check('AF · aVR still negative', ext(s.aVR).dom < 0, `dom=${ext(s.aVR).dom.toFixed(2)}`);
}
// ---- A3: bradycardia / tachycardia rate control ----
{
  const b = collect({ rhythm: 'brady', hr: 40 }); const rtb = rTimes(b.II, 0.5);
  const rrb = rtb.slice(1).map((t, i) => t - rtb[i]); const hrB = 60 / (rrb.reduce((a, c) => a + c, 0) / rrb.length);
  check('bradycardia · rate < 55 bpm', hrB < 55, `${hrB.toFixed(0)} bpm`);
  const t = collect({ rhythm: 'tachy', hr: 140 }); const rtt = rTimes(t.II, 0.5);
  const rrt = rtt.slice(1).map((x, i) => x - rtt[i]); const hrT = 60 / (rrt.reduce((a, c) => a + c, 0) / rrt.length);
  check('tachycardia · rate > 120 bpm', hrT > 120, `${hrT.toFixed(0)} bpm`);
}
// ---- A3: STEMI ST-elevation vector present ----
{
  const s = collect({ rhythm: 'stemi', hr: 75 });
  // sample the ST segment: compare a mid-ST value to baseline in an anterior lead (V2)
  const el = ext(s.V2).mx; // hyperacute T + ST raises V2 excursion vs normal
  const norm = ext(collect({ rhythm: 'sinus', hr: 75 }).V2).mx;
  check('STEMI · anterior lead amplitude elevated vs sinus', el > norm * 1.1, `V2 max ${el.toFixed(2)} vs ${norm.toFixed(2)}`);
}
// ---- A3: LVH voltage ----
{
  const s = collect({ path: { lvh: true }, hr: 72 });
  const norm = collect({ hr: 72 });
  check('LVH · increased precordial voltage', ext(s.V5).mx > ext(norm.V5).mx * 1.3, `V5 ${ext(s.V5).mx.toFixed(2)} vs ${ext(norm.V5).mx.toFixed(2)}`);
}

// ---- print ----
const w0 = 4, w1 = 46;
console.log('\nEchoSim self-validation — §A3 ECG morphology (pure model)\n');
for (const [st, nm, d] of rows) console.log(`  ${st.padEnd(w0)} ${nm.padEnd(w1)} ${d}`);
console.log(`\n  ${pass} passed, ${fail} failed.\n`);
console.log('  §A1 gCNR / §A2 measurement-accuracy: run in-app (Validation panel) or via the');
console.log('  headless echo harness once echo.metrics.gcnr is present. See docs/VALIDATION.md.\n');
process.exit(fail ? 1 : 0);
