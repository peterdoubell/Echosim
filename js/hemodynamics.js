// hemodynamics.js — a lumped-parameter closed-loop circulation model.
//
// Replaces the old kinematic "scale the cavity by a phase curve" hack with a
// real, self-consistent circulation: a time-varying-elastance left ventricle
// (E(t)·(V−V0)) coupled to a two-element Windkessel systemic afterload (R,C),
// a mitral inflow resistance fed by a left-atrial/venous filling pressure, and
// an aortic valve that opens only when Plv > Pao. Integrating this to its limit
// cycle yields a physiological LV pressure–volume loop, from which EDV / ESV /
// stroke volume / EF and the instantaneous transvalvular gradients emerge —
// rather than being hand-set. Pathology shifts real parameters (contractility,
// valve resistance, afterload, filling pressure) so the loop moves the way a
// diseased heart's loop actually moves.
//
// Everything is pressure in mmHg, volume in mL, time in s, flow in mL/s,
// elastance in mmHg/mL, compliance in mL/mmHg, resistance in mmHg·s/mL. The
// heavy integration runs ONCE per pathology and is memoised; per-frame callers
// (geometryAt) hit a cheap sampled/interpolated trace.
//
// References: Suga–Sagawa time-varying elastance; Stergiopulos two-element
// Windkessel; standard closed-loop lumped circulation (CircAdapt-style, reduced).

import { clamp } from './mathutils.js';

const CYCLE = 0.833;          // reference cycle length (s) at ~72 bpm
// Ventricular activation is ASYMMETRIC: contraction to peak elastance takes
// longer than the relaxation that follows, and end-systole sits at the peak. A
// symmetric bump over 0..0.40 put peak elastance at tau 0.20, which ended
// ejection there and gave an LV ejection time of only ~125 ms (normal 280-320)
// with a peak aortic flow of ~1190 mL/s (normal ~450) — the stroke volume was
// right but it was expelled far too fast, so the LV then sat frozen at ESV for
// the rest of systole. TC puts peak elastance (end-systole, aortic closure) at
// ~333 ms post-QRS. TR sets the isovolumic pressure fall: 0.12 gives an IVRT of
// ~85 ms (normal 70-100) with mitral opening at phase ~0.50; a slower decay
// (0.22) stretched IVRT to ~150 ms, the impaired-relaxation pattern.
const TC = 0.40;              // contraction: activation rises to peak (end-systole)
const TR = 0.14;              // relaxation: activation decays back to zero
const NSAMPLE = 240;          // samples stored per output cycle
const STEPS = 3000;           // integration steps per cycle
const NCYCLE = 10;            // cycles run to reach the limit cycle

// Ventricular activation e(τ) in [0,1]: a raised-sine rise to peak elastance at
// τ = TC (end-systole), then a faster raised-cosine decay to zero at TC+TR, zero
// through the rest of diastole. τ is the cardiac phase in [0,1) (0 == onset of
// systole). Both joins are C1 (zero slope at τ=0, τ=TC and τ=TC+TR).
function activation(tau) {
  if (tau < TC) { const u = tau / TC, sn = Math.sin(Math.PI * 0.5 * u), cs = 0.5 * (1 - Math.cos(Math.PI * u)); const w = 0.4; return sn * (1 - w) + cs * w; }
  if (tau < TC + TR) { const c = Math.cos(Math.PI * 0.5 * (tau - TC) / TR); return c * c; }
  return 0;
}

// Left-atrial pressure source: a resting reservoir pressure plus a late-diastolic
// a-wave (atrial kick ~τ 0.9). Elevated resting LAP models MS / MR congestion.
function atrialP(tau, p) {
  // the a-wave peaks just before the QRS (atrial systole ends at mitral closure,
  // ~40-60 ms before ejection); wrap-safe so its tail continues across tau = 0
  let dd = tau - (0.95); dd -= Math.round(dd);
  const d = dd / 0.045;
  const aWave = Math.exp(-d * d);
  return p.Pla0 + p.aWave * aWave;
}

// Per-pathology circulation parameters. Normal values are calibrated to
// EDV≈121, ESV≈47, EF≈61%, aortic ~125/76, LVEDP≈8.
function circParams(path) {
  const p = {
    Emax: 2.7,   // end-systolic elastance (contractility)
    Emin: 0.085, // diastolic (passive) elastance / stiffness
    V0: 10,      // unstressed LV volume
    Pla0: 8.0,   // resting left-atrial / filling pressure (raised → stronger E wave)
    aWave: 4.4,  // atrial-kick pressure bump (tuned → E-dominant normal filling, E/A≈1.1)
    Rmv: 0.010,  // mitral (inflow) resistance
    // Rao 0.018 gave a peak LV→aorta gradient of ~14 mmHg in a NORMAL heart,
    // which on echo reads as mild aortic stenosis; a normal valve is 2-5 mmHg.
    // AS overrides Rao explicitly per grade below, so this only affects normals.
    Rao: 0.006,
    // (aortic-valve resistance, ∝ 1/AVA)
    Rsys: 0.95,  // systemic vascular resistance
    Cao: 1.50,   // arterial (Windkessel) compliance (sets pulse pressure)
    Pven: 4,     // downstream systemic venous pressure
    Rmr: Infinity, // mitral regurgitant-orifice resistance (∞ = competent)
  };
  // graded severity (mild/moderate/severe) for valve lesions; default severe so
  // presets that set only the flag keep their original "severe by design" behaviour.
  const g = gradeOf(path);
  if (path.aorticStenosis) {                                      // small AVA → LV pressure overload
    p.Rao = ({ mild: 0.068, moderate: 0.125, severe: 0.24 })[g];
    p.Emax = ({ mild: 3.85, moderate: 4.3, severe: 5.15 })[g];   // compensatory, well above normal (2.7)
    // compensated AS keeps a normal-to-high blood pressure (~120-140/70-80): the
    // systemic vascular resistance is raised, so the LV (not the aorta) carries
    // the valve gradient — LV systolic pressure ~190-210 in severe AS
    p.Rsys *= 1.3;
  }
  if (path.lvh && !path.dilated) { p.Emin = 0.12; }               // stiff, hypertrophied → higher LVEDP
  if (path.dilated) { p.Emax = 0.60; p.Emin = 0.060; p.V0 = 30; p.Rsys = 1.30; p.Cao = 1.4; p.Pla0 = 13; }
  // DCM: annular dilatation and leaflet tethering leave a central FUNCTIONAL MR
  // (moderate: regurgitant fraction ~30 %) even with structurally normal leaflets
  if (path.dilated && !path.mr) { p.Rmr = 1.2; }
  if (path.mr) {
    // Regurgitant orifice + congested LA, plus the chronic volume-overload
    // REMODELLING of primary MR: the LV dilates (unstressed volume V0 rises) so it
    // is not emptied to nothing through the low-resistance leak — EF stays
    // high-normal (~65-75 %) with a dilated, hyperdynamic ventricle. Calibrated to
    // ASE grading: regurgitant fraction ~24 / 35 / 54 %, volume ~22 / 40 / 80 mL.
    // Chronic MR keeps the forward stroke volume by ejecting a larger total volume
    // over a near-normal time, not by emptying faster: with normal contractility
    // the unloaded LV reached its end-systolic volume ~190 ms into severe-MR
    // ejection at a peak aortic flow ~1.3x normal, so the cusps opened for only
    // two-thirds of systole. The (subclinically) lower end-systolic elastance of
    // chronic volume overload gives LVET ~230 ms at a near-normal peak aortic flow.
    p.Rmr = ({ mild: 1.6, moderate: 0.9, severe: 0.38 })[g];
    p.Pla0 = ({ mild: 9, moderate: 11, severe: 14 })[g];
    p.V0 = ({ mild: 20, moderate: 34, severe: 40 })[g];
    p.Emax = Math.min(p.Emax, ({ mild: 2.7, moderate: 2.4, severe: 1.8 })[g]);
    p.aWave = 6;
  }
  if (path.mitralStenosis) {                                      // narrow inflow, high LAP
    p.Rmv = ({ mild: 0.043, moderate: 0.078, severe: 0.128 })[g];
    p.Pla0 = ({ mild: 12, moderate: 16, severe: 20 })[g];
    p.aWave = 8;
  }
  if (path.rvpo) {
    // pulmonary hypertension: the failing right heart delivers less to the left
    // (low LA filling pressure) and the flattened, leftward-shifted septum stiffens
    // LV filling (ventricular interdependence) — the D-shaped LV is underfilled
    p.Pla0 = Math.min(p.Pla0, 7); p.Emin *= 1.3;            // LV EDV ~90 mL
  }
  if (path.rwma) { p.Emax = Math.min(p.Emax, 1.8); }              // regional dysfunction → low global contractility
  return p;
}

// Mitral valve area (cm^2) per grade: the ONE area the leaflet geometry (planimetry)
// is drawn to. circParams' Rmv/Pla0 are tuned so the Hakki area of the solved
// circulation agrees with it (verify-anatomy checks this).
export const MS_AREA = { mild: 1.8, moderate: 1.3, severe: 0.95 };

// Normalise a pathology's severity grade to 'mild' | 'moderate' | 'severe'.
// Default 'severe' preserves the original single-severity presets.
export function gradeOf(path = {}) {
  const gr = path.grade;
  return (gr === 'mild' || gr === 'moderate' || gr === 'severe') ? gr : 'severe';
}

// Integrate the closed loop to its limit cycle and return a sampled trace of one
// cycle plus scalar summary (EDV/ESV/SV/EF/pressures/gradients/regurg fraction).
function simulate(path) {
  const p = circParams(path);
  const dt = CYCLE / STEPS;
  let Vlv = 120, Pao = 80, Qao = 0;
  // LVOT/aortic inertance (mmHg·s²/mL) and aortic characteristic impedance (mmHg·s/mL)
  const L = 0.0015, Zc = 0.08;

  const tr = {
    phase: new Float64Array(NSAMPLE), Plv: new Float64Array(NSAMPLE),
    Pao: new Float64Array(NSAMPLE), Vlv: new Float64Array(NSAMPLE),
    Qmv: new Float64Array(NSAMPLE), Qao: new Float64Array(NSAMPLE),
    Qmr: new Float64Array(NSAMPLE), Pla: new Float64Array(NSAMPLE),
  };

  for (let cyc = 0; cyc < NCYCLE; cyc++) {
    const record = cyc === NCYCLE - 1;
    let sIdx = 0;
    for (let s = 0; s < STEPS; s++) {
      const tau = s / STEPS;
      const a = activation(tau);
      const E = p.Emin + (p.Emax - p.Emin) * a;
      const Plv = E * (Vlv - p.V0);
      const Pla = atrialP(tau, p);
      // valves (diode + resistance); mitral opens when LA>LV (diastole)
      const Qmv = Pla > Plv ? (Pla - Plv) / p.Rmv : 0;
      // aortic outflow has INERTIA (blood mass in the LVOT / proximal aorta) and
      // meets the aorta's characteristic impedance: the flow accelerates over
      // ~80-100 ms to a mid-systolic peak instead of jumping to its maximum the
      // instant the valve opens. Semi-implicit (stable) update; the valve closes
      // when the forward flow decelerates to zero.
      if (Qao > 0 || Plv > Pao) {
        Qao = (Qao + dt / L * (Plv - Pao)) / (1 + dt * (p.Rao + Zc) / L);
        if (Qao < 0) Qao = 0;
      }
      const Qmr = (p.Rmr !== Infinity && Plv > Pla) ? (Plv - Pla) / p.Rmr : 0;
      const Qsys = (Pao - p.Pven) / p.Rsys;
      Vlv += (Qmv - Qao - Qmr) * dt;
      Pao += ((Qao - Qsys) / p.Cao) * dt;
      if (Vlv < 0) Vlv = 0;

      if (record) {
        // write this step into every sample slot it covers
        const target = Math.floor((s + 1) * NSAMPLE / STEPS);
        while (sIdx < target && sIdx < NSAMPLE) {
          tr.phase[sIdx] = sIdx / NSAMPLE;
          tr.Plv[sIdx] = Plv; tr.Pao[sIdx] = Pao + Zc * Qao; tr.Vlv[sIdx] = Vlv;   // Pao: ascending-aortic pressure
          tr.Qmv[sIdx] = Qmv; tr.Qao[sIdx] = Qao; tr.Qmr[sIdx] = Qmr; tr.Pla[sIdx] = Pla;
          sIdx++;
        }
      }
    }
  }

  // summary scalars from the recorded cycle
  let EDV = -1e9, ESV = 1e9, PlvSys = -1e9, PaoSys = -1e9, PaoDia = 1e9;
  let grad = 0, fwd = 0, regurg = 0, tMin = 0, PlaAtSys = 0, gradMV = 0, mvSum = 0, mvN = 0;
  let QaoMax = 0, QmvMax = 0, QmrMax = 0, iQmax = 0, iOpen = -1, iClose = -1;
  const dtS = CYCLE / NSAMPLE;
  for (let i = 0; i < NSAMPLE; i++) {
    if (tr.Vlv[i] > EDV) EDV = tr.Vlv[i];
    if (tr.Vlv[i] < ESV) { ESV = tr.Vlv[i]; tMin = tr.phase[i]; }
    if (tr.Plv[i] > PlvSys) { PlvSys = tr.Plv[i]; PlaAtSys = tr.Pla[i]; }
    if (tr.Pao[i] > PaoSys) PaoSys = tr.Pao[i];
    if (tr.Pao[i] < PaoDia) PaoDia = tr.Pao[i];
    if (tr.Qao[i] > 1 && tr.Qao[i] * p.Rao > grad) grad = tr.Qao[i] * p.Rao;   // transvalvular drop
    if (tr.Qao[i] > QaoMax) { QaoMax = tr.Qao[i]; iQmax = i; }
    if (tr.Qao[i] > 1) { if (iOpen < 0) iOpen = i; iClose = i; }
    if (tr.Qmv[i] > 1) {                                   // transmitral LA-LV gradient during filling
      const dm = tr.Pla[i] - tr.Plv[i];
      if (dm > gradMV) gradMV = dm;
      mvSum += dm; mvN++;
    }
    if (tr.Qmv[i] > QmvMax) QmvMax = tr.Qmv[i];
    if (tr.Qmr[i] > QmrMax) QmrMax = tr.Qmr[i];
    fwd += tr.Qao[i] * dtS;
    regurg += tr.Qmr[i] * dtS;
  }
  const SV = EDV - ESV;
  const sum = {
    EDV, ESV, SV, EF: EDV > 0 ? (SV / EDV) * 100 : 0,
    PlvSys, PaoSys, PaoDia, PlaAtSys, gradient: grad,
    // transmitral gradient (mmHg): peak and diastolic mean, and the Hakki valve area
    // (cardiac output L/min over sqrt(mean gradient)) — what continuity/PHT would report
    gradientMV: gradMV, meanGradMV: mvN ? mvSum / mvN : 0,
    mvaHakki: mvN && mvSum > 0 ? (SV * 60 / CYCLE / 1000) / Math.sqrt(mvSum / mvN) : 0,
    forwardSV: fwd, regurgVol: regurg,
    regurgFraction: SV > 0 ? clamp(regurg / SV, 0, 1) : 0,
    tMinVol: tMin,
    QaoMax, QmvMax, QmrMax, // per-cycle peak flows, for jet-velocity envelopes
    accelTime: (iQmax - iOpen) * dtS, ejectTime: (iClose - iOpen + 1) * dtS,
    Pla0: p.Pla0,           // chronic LA pressure (drives LA remodelling in the anatomy)
  };
  return { tr, sum, params: p };
}

// ---- memoisation ---------------------------------------------------------
// The integration is heavy but depends only on the pathology flags, so cache by
// a signature of the mechanically-relevant flags (identical to how echo.js
// memoises its fixed-phase geometries). Bounded: one entry per distinct case.
const _cache = new Map();
function sig(path) {
  return [path.aorticStenosis, path.lvh, path.dilated, path.mr, path.mitralStenosis,
    path.rwma ? 1 : 0, path.rvpo].map((v) => (v ? 1 : 0)).join('') + ':' + gradeOf(path);
}
function solved(path) {
  const key = sig(path);
  let r = _cache.get(key);
  if (!r) { r = simulate(path); _cache.set(key, r); }
  return r;
}

// Reusable singleton so the ~per-frame call in geometryAt allocates nothing.
const _state = {
  phase: 0, Plv: 0, Pao: 0, Vlv: 0, Pla: 0, Qmv: 0, Qao: 0, Qmr: 0,
  rho: 1, kv: 0, sum: null,
};

// Instantaneous circulation state at a cardiac phase (linear-interpolated from
// the cached limit-cycle trace). `rho` = Vlv/EDV (cavity-volume ratio that drives
// the geometry); `kv` = (EDV−Vlv)/(EDV−ESV) normalised contraction fraction.
export function hemodynamics(phase, path = {}) {
  const { tr, sum } = solved(path);
  const x = ((phase % 1) + 1) % 1 * NSAMPLE;
  const i0 = Math.floor(x) % NSAMPLE;
  const i1 = (i0 + 1) % NSAMPLE;
  const f = x - Math.floor(x);
  const lerp = (arr) => arr[i0] + (arr[i1] - arr[i0]) * f;
  _state.phase = phase;
  _state.Plv = lerp(tr.Plv); _state.Pao = lerp(tr.Pao); _state.Vlv = lerp(tr.Vlv);
  _state.Pla = lerp(tr.Pla); _state.Qmv = lerp(tr.Qmv); _state.Qao = lerp(tr.Qao);
  _state.Qmr = lerp(tr.Qmr);
  _state.rho = sum.EDV > 0 ? clamp(_state.Vlv / sum.EDV, 0.05, 1) : 1;
  _state.kv = sum.SV > 0 ? clamp((sum.EDV - _state.Vlv) / sum.SV, 0, 1.2) : 0;
  _state.sum = sum;
  return _state;
}

// Scalar summary of the loop (EDV/ESV/SV/EF/pressures/gradient/regurg fraction).
export function hemoSummary(path = {}) { return solved(path).sum; }

// Mitral-inflow E/A ratio from the modelled diastolic flow: peak early (E, rapid
// filling just after mitral opening) over peak late (A, atrial kick). Normal filling
// is E-dominant (≈ 0.8–1.5); it falls below 1 with impaired relaxation. Memoised.
let _eaPath = null, _eaVal;
export function eaRatio(path = {}) {
  if (path === _eaPath) return _eaVal;
  const { tr } = solved(path);
  let e = 0, a = 0;
  for (let i = 0; i < NSAMPLE; i++) {
    const ph = i / NSAMPLE, q = tr.Qmv[i];
    if (ph >= 0.36 && ph < 0.62) { if (q > e) e = q; }
    else if (ph >= 0.80) { if (q > a) a = q; }
  }
  _eaVal = a > 0 ? e / a : null;
  _eaPath = path;
  return _eaVal;
}

// A sampled PV loop for one cardiac cycle: [{V, P}, …] (V in mL, P in mmHg),
// for the UI to draw the pressure–volume loop. n defaults to the native trace.
export function pvLoop(path = {}, n = 120) {
  const { tr } = solved(path);
  const out = [];
  for (let i = 0; i < n; i++) {
    const idx = Math.round(i * NSAMPLE / n) % NSAMPLE;
    out.push({ V: tr.Vlv[idx], P: tr.Plv[idx], phase: tr.phase[idx] });
  }
  return out;
}
