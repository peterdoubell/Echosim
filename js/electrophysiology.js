// electrophysiology.js — a physiologically-structured 12-lead ECG.
//
// Instead of a hand-drawn P-QRS-T bump, this derives the ECG from a moving
// cardiac dipole (a vectorcardiogram, VCG): atrial depolarisation, the septal→
// free-wall→basal ventricular activation sequence, and ventricular repolarisation
// are each modelled as time-localised 3D dipole generators along their true
// anatomic axes. The instantaneous heart vector is then projected onto the
// standard leads with the (inverse) Dower lead-field transform to yield all
// 12 leads simultaneously — the same VCG→12-lead mapping used clinically.
//
// A rhythm engine schedules beats (sinus + HRV, AF, brady/tachy, AV block,
// ventricular ectopy) and pathology modifiers reshape the morphology (LVH, LBBB/
// RBBB, STEMI ST-shift, low-voltage effusion, hyperkalaemia). The ventricular
// activation time also drives mechanical contraction (electromechanical coupling),
// so the echo and the ECG are one coupled system.
//
// References: dipole/Gaussian ECG model (McSharry/Sameni); inverse-Dower VCG→12-
// lead transform (Edenbrandt & Pahlm). Amplitudes in mV, times in seconds.

import { clamp } from './mathutils.js';

// Inverse-Dower coefficients: lead = row · [X, Y, Z]  (Frank VCG axes:
// X = right→left, Y = superior→inferior, Z = posterior→anterior).
const DOWER = {
  V1: [-0.515, 0.157, -0.917], V2: [0.044, 0.164, -1.387], V3: [0.882, 0.098, -1.277],
  V4: [1.213, 0.127, -0.601], V5: [1.125, 0.127, -0.086], V6: [0.831, 0.076, 0.230],
  I: [0.632, -0.235, 0.059], II: [0.235, 1.066, -0.132],
};
export const LEAD_NAMES = ['I', 'II', 'III', 'aVR', 'aVL', 'aVF', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6'];

// Project a buffered VCG sample [X,Y,Z] to all 12 leads (for the strip / grid).
export function vcgToLeads(v) { return leadsFromVCG(v, {}); }

// Project a VCG vector [X,Y,Z] to all 12 leads (III/aVR/aVL/aVF from I & II).
function leadsFromVCG(v, out) {
  const dot = (r) => r[0] * v[0] + r[1] * v[1] + r[2] * v[2];
  const I = dot(DOWER.I), II = dot(DOWER.II);
  out.I = I; out.II = II; out.III = II - I;
  out.aVR = -(I + II) / 2; out.aVL = I - II / 2; out.aVF = II - I / 2;
  out.V1 = dot(DOWER.V1); out.V2 = dot(DOWER.V2); out.V3 = dot(DOWER.V3);
  out.V4 = dot(DOWER.V4); out.V5 = dot(DOWER.V5); out.V6 = dot(DOWER.V6);
  return out;
}

// A VCG "wave" = a Gaussian-in-time 3D dipole generator. t0 relative to the R
// peak (s); s = width (s); a = [X,Y,Z] peak dipole (mV-equivalent).
function gauss(t, t0, s) { const d = (t - t0) / s; return Math.exp(-0.5 * d * d); }

// Baseline normal-sinus wave set (adult). QRS is three overlapping generators
// (septal q, main R, terminal s) whose vector loop gives physiologic lead shapes.
function normalWaves() {
  return {
    P:  { t0: -0.16, s: 0.022, a: [0.11, 0.06, 0.04] },   // atrial: left-inferior-anterior, small
    Q:  { t0: -0.020, s: 0.010, a: [-0.10, -0.02, 0.10] },  // septal r→l reversal: rightward-anterior
    R:  { t0: 0.000, s: 0.018, a: [1.30, 0.60, -0.18] },    // main: left-inferior, mildly posterior (peak R ~V4-V5)
    S:  { t0: 0.028, s: 0.014, a: [-0.28, -0.12, -0.22] },  // terminal: right-superior-posterior
    T:  { t0: 0.230, s: 0.070, a: [0.36, 0.26, -0.03] },    // repolarisation: concordant with R
    st: 0.0,                                                // ST-segment offset (STEMI), scalar→vector below
    stVec: [0, 0, 0],
  };
}

// Apply pathology morphology to a wave set (mutates & returns).
function applyMorphology(W, path, rhythm, beatType) {
  // --- ventricular ectopic beat: wide, bizarre, discordant, no P ---
  if (beatType === 'pvc') {
    W.P.a = [0, 0, 0];
    W.R = { t0: 0.0, s: 0.045, a: [-1.1, -0.7, 0.5] };   // wide, opposite axis
    W.S = { t0: 0.06, s: 0.03, a: [0.5, 0.4, 0.2] };
    W.Q.a = [0, 0, 0];
    W.T = { t0: 0.34, s: 0.10, a: [0.8, 0.55, -0.2] };   // large discordant T
    return W;
  }
  // --- rhythm-level P-wave handling ---
  if (rhythm === 'afib') W.P.a = [0, 0, 0];               // no organised P
  // --- conduction: bundle branch blocks widen & re-vector QRS ---
  if (path.lbbb || rhythm === 'lbbb') {
    W.R = { t0: 0.010, s: 0.040, a: [1.5, 0.5, 0.1] };    // wide, notched, no septal q
    W.Q.a = [0, 0, 0];
    W.S = { t0: 0.055, s: 0.028, a: [-0.2, 0.0, -0.5] };
    W.T = { t0: 0.30, s: 0.09, a: [-0.6, -0.35, 0.1] };   // discordant
  }
  if (path.rbbb || rhythm === 'rbbb') {
    W.S = { t0: 0.055, s: 0.030, a: [-0.6, -0.1, 0.6] };  // terminal rightward-anterior (rSR' in V1)
    W.R2 = { t0: 0.045, s: 0.020, a: [0.3, 0.0, 0.7] };   // R' anterior
  }
  // --- LVH (AS / hypertension): tall QRS + strain (down-sloping ST/T) ---
  if (path.lvh || path.aorticStenosis) {
    W.R.a = [W.R.a[0] * 1.7, W.R.a[1] * 1.6, W.R.a[2] * 1.7];
    W.T.a = [-0.35, -0.2, 0.05];                          // T inversion / strain
    W.stVec = [-0.06, -0.04, 0.02];
  }
  // --- dilated cardiomyopathy: often LBBB-like + low-ish voltage ---
  if (path.dilated) {
    W.R.a = [W.R.a[0] * 0.85, W.R.a[1] * 0.85, W.R.a[2] * 0.85];
    W.R.s = 0.030;                                        // mild widening
  }
  // --- pericardial effusion: low voltage (± electrical alternans handled by beat) ---
  if (path.effusion) {
    for (const k of ['P', 'Q', 'R', 'S', 'T']) W[k].a = W[k].a.map((v) => v * 0.5);
  }
  // --- STEMI (explicit rhythm option): anterior ST elevation vector ---
  if (rhythm === 'stemi') {
    W.stVec = [0.5, 0.15, -0.9];                          // elevation projecting to V1-V4
    W.T.a = [0.7, 0.3, -1.0];                             // hyperacute T
  }
  // --- hyperkalaemia: peaked narrow T + widened QRS ---
  if (rhythm === 'hyperk') {
    W.T = { t0: 0.20, s: 0.035, a: [0.9, 0.6, -0.1] };
    W.R.s = 0.030;
    W.P.a = W.P.a.map((v) => v * 0.4);
  }
  return W;
}

// Sum the wave generators into a VCG vector at time `tr` relative to the beat's
// R peak. ST offset is applied as a plateau between S and T.
function vcgAtBeat(tr, W, out) {
  out[0] = out[1] = out[2] = 0;
  const add = (wv) => {
    if (!wv) return;
    const g = gauss(tr, wv.t0, wv.s);
    out[0] += wv.a[0] * g; out[1] += wv.a[1] * g; out[2] += wv.a[2] * g;
  };
  add(W.P); add(W.Q); add(W.R); add(W.R2); add(W.S); add(W.T);
  // ST segment: smooth plateau from ~end-QRS to T onset
  const stW = 0.5 * (1 + Math.tanh((tr - 0.05) / 0.02)) * 0.5 * (1 + Math.tanh((0.18 - tr) / 0.04));
  out[0] += W.stVec[0] * stW; out[1] += W.stVec[1] * stW; out[2] += W.stVec[2] * stW;
  return out;
}

// ---------------------------------------------------------------------------
// Heart-rate warp of the mechanical cycle.
//
// The circulation (hemodynamics.js) is solved once, at the reference cycle of
// 0.833 s (72 bpm), and indexed by a phase whose 0 is the onset of contraction.
// Dividing the time since contraction onset by the R-R stretched every interval
// in proportion: at 50 bpm LV ejection lasted ~425 ms and aortic closure came
// ~110 ms after the end of the T wave; at 120 bpm ejection lasted ~175 ms and
// the valve shut ~75 ms before it. In the heart the isovolumic periods barely
// change with rate, ejection shortens only modestly (Weissler: LVET ≈ 413 −
// 1.7·HR ms) and diastole takes up most of the R-R change, with atrial systole
// locked to the next QRS (a constant PR interval). So the post-QRS time is
// mapped piecewise onto the reference phases of the normal trace:
//   IVCT  (MVC → AVO)   fixed duration
//   LVET  (AVO → AVC)   scaled by LVET(HR) / LVET(72)
//   IVRT  (AVC → MVO)   fixed duration
//   early diastole      absorbs the R-R change (rapid filling + diastasis)
//   atrial systole      fixed ~167 ms before the next contraction; in tachycardia
//                       it shares the short diastole with early filling (E/A
//                       approaching fusion)
// At 72 bpm the map is the identity. tools/verify-timing.mjs checks the
// breakpoints against the solved trace and aortic closure against the T wave.
// ---------------------------------------------------------------------------
export const REF_RR = 0.833;           // reference cycle (s) — hemodynamics CYCLE
export const PH_AVO = 0.072;           // aortic opening (normal trace)
export const PH_AVC = 0.421;           // aortic closure
export const PH_MVO = 0.513;           // mitral opening
export const PH_ATR = 0.80;            // start of atrial systole (A-wave onset)
// Weissler's LV ejection time (s) for a heart rate (bpm), clamped to a sane range.
export function lvetAt(hr) { return clamp(0.413 - 0.0017 * hr, 0.18, 0.36); }
const LVET_REF = lvetAt(60 / REF_RR);

// Reference mechanical phase [0, 1) at time tp (s) since contraction onset
// (QRS + electromechanical delay) in a beat whose R-R is rr (s).
export function mechPhase(tp, rr) {
  let ivct = PH_AVO * REF_RR;
  let ej = (PH_AVC - PH_AVO) * REF_RR * lvetAt(60 / rr) / LVET_REF;
  let ivrt = (PH_MVO - PH_AVC) * REF_RR;
  // very fast rates: keep at least ~28 % of the cycle for filling
  const sys = ivct + ej + ivrt, cap = 0.72 * rr;
  if (sys > cap) { const f = cap / sys; ivct *= f; ej *= f; ivrt *= f; }
  const dia = rr - (ivct + ej + ivrt);
  const atr = Math.min((1 - PH_ATR) * REF_RR, 0.5 * dia);
  const early = dia - atr;
  // piecewise-linear: [segment duration, phase at start, phase at end]
  const t = Math.max(0, tp);
  let t0 = 0;
  const segs = [[ivct, 0, PH_AVO], [ej, PH_AVO, PH_AVC], [ivrt, PH_AVC, PH_MVO],
    [early, PH_MVO, PH_ATR], [atr, PH_ATR, 1]];
  for (const [d, p0, p1] of segs) {
    if (t < t0 + d) return clamp(p0 + (p1 - p0) * (t - t0) / d, 0, 0.999);
    t0 += d;
  }
  return 0.999;
}

// ---------------------------------------------------------------------------
// Rhythm engine — schedules beats and drives the mechanical phase.
// ---------------------------------------------------------------------------
export class Heartbeat {
  constructor() {
    this.hr = 72;
    this.rhythm = 'sinus';   // sinus | afib | brady | tachy | avblock1 | pvc | lbbb | rbbb | stemi | hyperk
    this.path = {};
    this._beats = [];        // recent beat onsets: {t, rr, type, W}
    this._t = 0;
    this._nextAt = 0;
    this._seed = 987654321;
    this._emDelay = 0.04;    // electromechanical delay (s): contraction lags QRS
    this._lastVcg = [0, 0, 0];
    this._leads = {};
    this._genBeat(0);        // seed first beat at t=0
  }
  _rand() { this._seed = (this._seed * 1103515245 + 12345) & 0x7fffffff; return this._seed / 0x7fffffff; }

  _baseRR() {
    const r = this.rhythm;
    let hr = this.hr;
    if (r === 'brady') hr = Math.min(hr, 46);
    if (r === 'tachy') hr = Math.max(hr, 130);
    return 60 / hr;
  }

  // Decide the next beat's RR and type from the current rhythm.
  _scheduleNext(prevType) {
    const base = this._baseRR();
    const r = this.rhythm;
    if (r === 'afib') {
      // irregularly irregular: RR uniformly jittered ±35%
      return { rr: base * (0.65 + 0.7 * this._rand()), type: 'normal' };
    }
    if (r === 'pvc') {
      // sinus with an interpolated PVC roughly every 4th beat + compensatory pause
      if (prevType !== 'pvc' && this._rand() < 0.28) return { rr: base * 0.62, type: 'pvc' };
      if (prevType === 'pvc') return { rr: base * 1.38, type: 'normal' }; // compensatory pause
    }
    // sinus (and morphology-only rhythms): base RR + respiratory sinus arrhythmia
    const rsa = 1 + 0.05 * Math.sin(this._t * 2 * Math.PI * 0.25);
    return { rr: base * rsa, type: 'normal' };
  }

  _genBeat(atT) {
    const prev = this._beats.length ? this._beats[this._beats.length - 1].type : 'normal';
    const sched = this._scheduleNext(prev);
    let W = applyMorphology(normalWaves(), this.path, this.rhythm, sched.type);
    // AV block: prolong PR (shift P earlier relative to R)
    if (this.rhythm === 'avblock1') W.P.t0 = -0.28;
    // QT scales with RR: the T peak moves with sqrt(RR) (Bazett-ish) and the wave
    // narrows with the cube root (Fridericia), so the end of the T wave keeps
    // pace with aortic closure (A2) as the rate changes
    const rrC = clamp(sched.rr, 0.3, 1.5) / 0.857;
    W.T.t0 *= Math.sqrt(rrC);
    W.T.s *= Math.cbrt(rrC);
    const beat = { t: atT, rr: sched.rr, type: sched.type, W };
    this._beats.push(beat);
    if (this._beats.length > 6) this._beats.shift();
    this._nextAt = atT + sched.rr;
    return beat;
  }

  // Advance absolute time and schedule beats as needed.
  advance(dt) {
    this._t += dt;
    while (this._t >= this._nextAt) this._genBeat(this._nextAt);
  }

  // Current 12-lead voltages (mV) at absolute time this._t.
  leads() {
    // sum contributions from recent beats (waves overlap across the R-R boundary)
    const v = [0, 0, 0], tmp = [0, 0, 0];
    for (const b of this._beats) {
      const tr = this._t - (b.t); // beat.t marks the R peak time
      if (tr < -0.35 || tr > 0.6) continue;
      vcgAtBeat(tr, b.W, tmp);
      v[0] += tmp[0]; v[1] += tmp[1]; v[2] += tmp[2];
    }
    // afib fibrillatory baseline in place of P
    if (this.rhythm === 'afib') {
      const f = 0.05 * (Math.sin(this._t * 2 * Math.PI * 7.3) + 0.6 * Math.sin(this._t * 2 * Math.PI * 11.7));
      v[0] += f; v[1] += f * 0.6;
    }
    this._lastVcg = v;
    return leadsFromVCG(v, this._leads);
  }

  // Mechanical beat phase (0..1) for the CURRENT beat, with the electromechanical
  // delay applied, so contraction follows depolarisation. Returns {phase, rr}.
  // The phase indexes the circulation, which is solved once at the reference
  // cycle (72 bpm), so time since the delayed QRS is WARPED onto it (mechPhase)
  // rather than divided by the R-R: systole keeps its physiological length at
  // any rate and diastole absorbs the rest.
  mechanical() {
    // find the beat whose window contains _t (nearest onset at/just before _t);
    // until its contraction starts (the electromechanical delay after the QRS)
    // the previous beat's late diastole — atrial systole, mitral closure — is
    // still playing, so it keeps the phase instead of a 40 ms freeze at 0
    let i = 0;
    for (let j = 0; j < this._beats.length; j++) if (this._beats[j].t <= this._t + 1e-6) i = j;
    let cur = this._beats[i];
    if (this._t - cur.t < this._emDelay && i > 0) cur = this._beats[i - 1];
    const phase = mechPhase(this._t - cur.t - this._emDelay, cur.rr);
    return { phase, rr: cur.rr, beatType: cur.type };
  }

  setHR(hr) { this.hr = hr; }
  setRhythm(r) { this.rhythm = r; }
  setPathology(p) { this.path = p || {}; }
  get vcg() { return this._lastVcg; }
}
