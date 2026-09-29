// tools/verify-timing.mjs — cardiac-cycle timing audit (pure model, no browser).
//
// Drives the rhythm engine at several heart rates and reads the valve events off
// the circulation the renderer actually uses (hemodynamics at the warped
// mechanical phase), then checks them against physiology:
//   - the heart-rate warp's breakpoints sit on the solved normal trace, and the
//     warp is the identity at the reference 72 bpm;
//   - LV ejection time follows Weissler (LVET ≈ 413 − 1.7·HR ms) while the
//     isovolumic periods barely change, so diastole absorbs the R-R change;
//   - aortic closure (A2) stays at the end of the T wave (|AVC − T-end| < 30 ms);
//   - the mechanical phase runs on through the electromechanical delay after
//     each QRS (no freeze at phase 0, no skipped late diastole);
//   - severe MR ejects over a near-normal time at a near-normal peak aortic flow;
//   - the right-heart valves run on their own clock: T1 after M1, the pulmonic
//     valve opening before and closing after the aortic (A2-P2 split).
//
// Usage: node tools/verify-timing.mjs      (exit 0 = pass, 1 = fail)
import { Heartbeat, mechPhase, lvetAt, REF_RR, PH_AVO, PH_AVC, PH_MVO } from '../js/electrophysiology.js';
import { hemodynamics, hemoSummary } from '../js/hemodynamics.js';
import { geometryAt } from '../js/cardiac-model.js';

let pass = 0, fail = 0;
const rows = [];
function check(name, ok, detail) { rows.push([ok ? 'PASS' : 'FAIL', name, detail]); ok ? pass++ : fail++; }
const ms = (s) => (s * 1000).toFixed(0);

// ---- valve events of the solved normal trace (reference phase) ----
function traceEvents(path = {}) {
  const N = 4000, ev = {};
  let pAo = hemodynamics(1 - 1 / N, path).Qao > 1, pMv = hemodynamics(1 - 1 / N, path).Qmv > 0.5;
  for (let i = 0; i < N; i++) {
    const p = i / N, h = hemodynamics(p, path), ao = h.Qao > 1, mv = h.Qmv > 0.5;
    if (ao && !pAo) ev.AVO = p; if (!ao && pAo) ev.AVC = p;
    if (mv && !pMv) ev.MVO = p; if (!mv && pMv) ev.MVC = p;
    pAo = ao; pMv = mv;
  }
  return ev;
}
{
  const ev = traceEvents();
  const off = Math.max(Math.abs(ev.AVO - PH_AVO), Math.abs(ev.AVC - PH_AVC), Math.abs(ev.MVO - PH_MVO));
  check('warp breakpoints match the solved trace (±0.01)', off <= 0.01,
    `AVO ${ev.AVO.toFixed(3)}/${PH_AVO}  AVC ${ev.AVC.toFixed(3)}/${PH_AVC}  MVO ${ev.MVO.toFixed(3)}/${PH_MVO}`);
  let dmax = 0;
  for (let t = 0; t < REF_RR; t += 0.001) dmax = Math.max(dmax, Math.abs(mechPhase(t, REF_RR) - Math.min(t / REF_RR, 0.999)));
  check('warp is the identity at 72 bpm', dmax < 2e-3, `max |Δphase| = ${dmax.toExponential(1)}`);
}

// ---- run the rhythm engine and measure per-beat intervals ----
const FS = 2000;
function run(hr, secs = 9) {
  const h = new Heartbeat();
  h.setHR(hr);
  const n = Math.round(FS * secs), t = new Float64Array(n), ph = new Float64Array(n), II = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    h.advance(1 / FS);
    II[i] = h.leads().II;
    ph[i] = h.mechanical().phase;
    t[i] = (i + 1) / FS;
  }
  // R peaks (lead II upstroke + local max)
  const R = [];
  for (let i = 1; i < n - 1; i++) {
    if (II[i - 1] < 0.5 && II[i] >= 0.5) {
      let j = i; while (j + 1 < n && II[j + 1] > II[j]) j++;
      if (!R.length || t[j] - R[R.length - 1] > 0.25) R.push(t[j]);
    }
  }
  // valve events in real time from the circulation at the warped phase
  const ev = { AVO: [], AVC: [], MVO: [], MVC: [] };
  let pAo = false, pMv = true, maxJump = 0, zeroRun = 0, maxZeroRun = 0;
  for (let i = 0; i < n; i++) {
    const s = hemodynamics(ph[i], {}), ao = s.Qao > 1, mv = s.Qmv > 0.5;
    if (i > 0) {
      if (ao && !pAo) ev.AVO.push(t[i]); if (!ao && pAo) ev.AVC.push(t[i]);
      if (mv && !pMv) ev.MVO.push(t[i]); if (!mv && pMv) ev.MVC.push(t[i]);
      let d = Math.abs(ph[i] - ph[i - 1]); d = Math.min(d, 1 - d);
      maxJump = Math.max(maxJump, d);
    }
    // (the very first beat has no previous one to continue, so it starts at 0)
    zeroRun = ph[i] === 0 && t[i] > 1 ? zeroRun + 1 : 0; maxZeroRun = Math.max(maxZeroRun, zeroRun);
    pAo = ao; pMv = mv;
  }
  // per-beat intervals for complete beats (skip the first, warm-up beat)
  const beats = [];
  for (let b = 1; b + 1 < R.length; b++) {
    const r0 = R[b], r1 = R[b + 1];
    const first = (arr, lo, hi) => arr.find((x) => x >= lo && x < hi);
    const mvc = first(ev.MVC, r0 - 0.05, r0 + 0.15), avo = first(ev.AVO, r0, r1);
    const avc = first(ev.AVC, r0, r1), mvo = first(ev.MVO, r0, r1);
    if ([mvc, avo, avc, mvo].some((x) => x == null)) continue;
    // T-wave end by the tangent method on lead II: steepest descent after the T
    // peak, extrapolated to the (isoelectric, zero) baseline
    const i0 = Math.round((r0 + 0.12) * FS) - 1, i1 = Math.round(Math.min(r0 + 0.6, r1 - 0.2) * FS) - 1;
    let ip = i0; for (let i = i0; i <= i1; i++) if (II[i] > II[ip]) ip = i;
    let is = ip, sMin = 0;
    for (let i = ip; i <= i1; i++) { const sl = (II[i + 1] - II[i - 1]) * FS / 2; if (sl < sMin) { sMin = sl; is = i; } }
    const tEnd = t[is] - II[is] / sMin;
    beats.push({ rr: r1 - r0, ivct: avo - mvc, lvet: avc - avo, ivrt: mvo - avc, dia: r1 - mvo, avc: avc - r0, tEnd: tEnd - r0 });
  }
  const mean = (k) => beats.reduce((a, x) => a + x[k], 0) / beats.length;
  const worst = (f) => beats.reduce((a, x) => Math.max(a, f(x)), 0);
  return { beats, mean, worst, maxJump, maxZeroMs: maxZeroRun / FS * 1000 };
}

const res = {};
for (const hr of [50, 72, 110, 120]) {
  const r = res[hr] = run(hr);
  const lvetErr = r.worst((x) => Math.abs(x.lvet - lvetAt(60 / x.rr)));
  check(`${hr} bpm · LVET follows Weissler (±25 ms)`, r.beats.length >= 3 && lvetErr < 0.025,
    `LVET ${ms(r.mean('lvet'))} ms vs ${ms(lvetAt(hr))} ms (worst Δ ${ms(lvetErr)} ms, ${r.beats.length} beats)`);
  const avcErr = r.worst((x) => Math.abs(x.avc - x.tEnd));
  check(`${hr} bpm · A2 at the end of the T wave (<30 ms)`, avcErr < 0.030,
    `AVC R+${ms(r.mean('avc'))} ms, T-end R+${ms(r.mean('tEnd'))} ms (worst Δ ${ms(avcErr)} ms)`);
  check(`${hr} bpm · mechanical phase continuous through the QRS`, r.maxJump < 0.02 && r.maxZeroMs < 5,
    `max step ${r.maxJump.toFixed(3)}, longest hold at 0: ${r.maxZeroMs.toFixed(1)} ms`);
}
{
  const a = res[50], b = res[72], c = res[110];
  const d = (k) => Math.max(Math.abs(a.mean(k) - b.mean(k)), Math.abs(c.mean(k) - b.mean(k)));
  check('IVCT and IVRT nearly rate-independent (±15 ms, 50-110 bpm)', d('ivct') < 0.015 && d('ivrt') < 0.015,
    `IVCT ${ms(a.mean('ivct'))}/${ms(b.mean('ivct'))}/${ms(c.mean('ivct'))} ms, IVRT ${ms(a.mean('ivrt'))}/${ms(b.mean('ivrt'))}/${ms(c.mean('ivrt'))} ms`);
  const share = (a.mean('dia') - c.mean('dia')) / (a.mean('rr') - c.mean('rr'));
  check('diastole absorbs most of the R-R change (>70 %)', share > 0.70, `${(share * 100).toFixed(0)} % of ΔRR (50 → 110 bpm)`);
}

// ---- severe MR: near-normal ejection time and peak aortic flow ----
{
  const n = hemoSummary({}), m = hemoSummary({ mr: true });
  check('severe MR · LVET ≥ 220 ms', m.ejectTime >= 0.22, `${ms(m.ejectTime)} ms (normal ${ms(n.ejectTime)})`);
  check('severe MR · peak aortic flow ≤ 1.2× normal', m.QaoMax <= 1.2 * n.QaoMax, `${m.QaoMax.toFixed(0)} vs ${n.QaoMax.toFixed(0)} mL/s`);
  check('severe MR · regurgitant fraction ≥ 50 %', m.regurgFraction >= 0.5, `RF ${(m.regurgFraction * 100).toFixed(0)} %, EF ${m.EF.toFixed(0)} %`);
}

// ---- right-heart valve timing (reference cycle) ----
{
  const N = 2000, ev = {};
  const open = (v) => v > 0.1;
  let prev = geometryAt(1 - 1 / N, {}).valves;
  let kMax = -1, kAt = 0, kRMax = -1, kRAt = 0;
  for (let i = 0; i < N; i++) {
    const p = i / N, G = geometryAt(p, {}), v = G.valves;
    for (const nm of ['mitral', 'tricuspid', 'aortic', 'pulmonic']) {
      if (open(v[nm]) && !open(prev[nm])) ev[nm + 'O'] = p;
      if (!open(v[nm]) && open(prev[nm])) ev[nm + 'C'] = p;
    }
    if (G.k > kMax) { kMax = G.k; kAt = p; }
    if (G.kRV > kRMax) { kRMax = G.kRV; kRAt = p; }
    prev = v;
  }
  const dt = (a, b) => { let d = (ev[a] - ev[b]) % 1; if (d > 0.5) d -= 1; if (d < -0.5) d += 1; return d * REF_RR; };
  const t1 = dt('tricuspidC', 'mitralC'), pvo = dt('aorticO', 'pulmonicO'), split = dt('pulmonicC', 'aorticC');
  check('T1 follows M1 by 20-30 ms', t1 >= 0.018 && t1 <= 0.032, `${ms(t1)} ms`);
  check('pulmonic valve opens before the aortic (5-20 ms)', pvo >= 0.005 && pvo <= 0.020, `${ms(pvo)} ms`);
  check('A2-P2 split 20-50 ms', split >= 0.020 && split <= 0.050, `${ms(split)} ms`);
  const rvet = dt('pulmonicC', 'pulmonicO'), lvet = dt('aorticC', 'aorticO');
  check('RV ejection longer than LV ejection', rvet > lvet, `RVET ${ms(rvet)} ms vs LVET ${ms(lvet)} ms`);
  check('RV longitudinal shortening peaks after the LV', kRAt > kAt && (kRAt - kAt) * REF_RR < 0.06, `LV k peak @ ${kAt.toFixed(3)}, RV @ ${kRAt.toFixed(3)}`);
}

console.log('\nEchoSim cardiac-cycle timing audit\n');
for (const [st, nm, d] of rows) console.log(`  ${st.padEnd(4)} ${nm.padEnd(56)} ${d}`);
console.log(`\n  ${pass} passed, ${fail} failed.\n`);
process.exit(fail ? 1 : 0);
