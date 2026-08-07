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
const TSYS = 0.40;            // systolic activation fraction of the cycle
const NSAMPLE = 240;          // samples stored per output cycle
const STEPS = 3000;           // integration steps per cycle
const NCYCLE = 10;            // cycles run to reach the limit cycle

// Ventricular activation e(τ) in [0,1]: a raised-sine bump over systole, zero
// through diastole. τ is the cardiac phase in [0,1) (0 == onset of systole).
function activation(tau) {
  if (tau >= TSYS) return 0;
  const s = Math.sin(Math.PI * tau / TSYS);
  return s * s;
}

// Left-atrial pressure source: a resting reservoir pressure plus a late-diastolic
// a-wave (atrial kick ~τ 0.9). Elevated resting LAP models MS / MR congestion.
function atrialP(tau, p) {
  const d = (tau - 0.90) / 0.045;
  const aWave = Math.exp(-d * d);
  return p.Pla0 + p.aWave * aWave;
}

// Per-pathology circulation parameters. Normal values are calibrated to
// EDV≈120, ESV≈50, EF≈58%, aortic 120/78, LVEDP≈8.
function circParams(path) {
  const p = {
    Emax: 3.0,   // end-systolic elastance (contractility)
    Emin: 0.085, // diastolic (passive) elastance / stiffness
    V0: 10,      // unstressed LV volume
    Pla0: 8.0,   // resting left-atrial / filling pressure (raised → stronger E wave)
    aWave: 4.4,  // atrial-kick pressure bump (tuned → E-dominant normal filling, E/A≈1.1)
    Rmv: 0.010,  // mitral (inflow) resistance
    Rao: 0.018,  // aortic-valve resistance (∝ 1/AVA)
    Rsys: 1.15,  // systemic vascular resistance
    Cao: 1.6,    // arterial (Windkessel) compliance
    Pven: 4,     // downstream systemic venous pressure
    Rmr: Infinity, // mitral regurgitant-orifice resistance (∞ = competent)
  };
  // graded severity (mild/moderate/severe) for valve lesions; default severe so
  // presets that set only the flag keep their original "severe by design" behaviour.
  const g = gradeOf(path);
  if (path.aorticStenosis) {                                      // small AVA → LV pressure overload
    p.Rao = ({ mild: 0.052, moderate: 0.090, severe: 0.14 })[g];
    p.Emax = ({ mild: 3.2, moderate: 3.4, severe: 3.6 })[g];
  }
  if (path.lvh && !path.dilated) { p.Emin = 0.12; }               // stiff, hypertrophied → higher LVEDP
  if (path.dilated) { p.Emax = 0.60; p.Emin = 0.060; p.V0 = 30; p.Rsys = 1.30; p.Cao = 1.4; p.Pla0 = 13; }
  if (path.mr) {                                                  // regurgitant orifice + congested LA
    p.Rmr = ({ mild: 0.90, moderate: 0.33, severe: 0.11 })[g];
    p.Pla0 = ({ mild: 8, moderate: 9, severe: 10 })[g];
    p.aWave = 6;
  }
  if (path.mitralStenosis) {                                      // narrow inflow, high LAP
    p.Rmv = ({ mild: 0.045, moderate: 0.065, severe: 0.090 })[g];
    p.Pla0 = ({ mild: 11, moderate: 14, severe: 16 })[g];
    p.aWave = 8;
  }
  if (path.rwma) { p.Emax = Math.min(p.Emax, 1.8); }              // regional dysfunction → low global contractility
  return p;
}

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
  let Vlv = 120, Pao = 80;

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
      const Qao = Plv > Pao ? (Plv - Pao) / p.Rao : 0;
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
          tr.Plv[sIdx] = Plv; tr.Pao[sIdx] = Pao; tr.Vlv[sIdx] = Vlv;
          tr.Qmv[sIdx] = Qmv; tr.Qao[sIdx] = Qao; tr.Qmr[sIdx] = Qmr; tr.Pla[sIdx] = Pla;
          sIdx++;
        }
      }
    }
  }

  // summary scalars from the recorded cycle
  let EDV = -1e9, ESV = 1e9, PlvSys = -1e9, PaoSys = -1e9, PaoDia = 1e9;
  let grad = 0, fwd = 0, regurg = 0, tMin = 0, PlaAtSys = 0;
  let QaoMax = 0, QmvMax = 0, QmrMax = 0;
  const dtS = CYCLE / NSAMPLE;
  for (let i = 0; i < NSAMPLE; i++) {
    if (tr.Vlv[i] > EDV) EDV = tr.Vlv[i];
    if (tr.Vlv[i] < ESV) { ESV = tr.Vlv[i]; tMin = tr.phase[i]; }
    if (tr.Plv[i] > PlvSys) { PlvSys = tr.Plv[i]; PlaAtSys = tr.Pla[i]; }
    if (tr.Pao[i] > PaoSys) PaoSys = tr.Pao[i];
    if (tr.Pao[i] < PaoDia) PaoDia = tr.Pao[i];
    if (tr.Qao[i] > 1 && tr.Plv[i] - tr.Pao[i] > grad) grad = tr.Plv[i] - tr.Pao[i];
    if (tr.Qao[i] > QaoMax) QaoMax = tr.Qao[i];
    if (tr.Qmv[i] > QmvMax) QmvMax = tr.Qmv[i];
    if (tr.Qmr[i] > QmrMax) QmrMax = tr.Qmr[i];
    fwd += tr.Qao[i] * dtS;
    regurg += tr.Qmr[i] * dtS;
  }
  const SV = EDV - ESV;
  const sum = {
    EDV, ESV, SV, EF: EDV > 0 ? (SV / EDV) * 100 : 0,
    PlvSys, PaoSys, PaoDia, PlaAtSys, gradient: grad,
    forwardSV: fwd, regurgVol: regurg,
    regurgFraction: SV > 0 ? clamp(regurg / SV, 0, 1) : 0,
    tMinVol: tMin,
    QaoMax, QmvMax, QmrMax, // per-cycle peak flows, for jet-velocity envelopes
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
    path.rwma ? 1 : 0].map((v) => (v ? 1 : 0)).join('') + ':' + gradeOf(path);
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
