// hemodynamics.js — a lumped-parameter closed-loop circulation model.
//
// Replaces the old kinematic "scale the cavity by a phase curve" hack with a
// real, self-consistent circulation: a time-varying-elastance left ventricle
// (E(t)·(V−V0)) coupled to a two-element Windkessel systemic afterload (R,C),
// a mitral inflow (resistance + inertance) fed by a compliant, contracting left
// atrium that fills from the pulmonary veins, and
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

// Left-atrial activation (atrial systole): a Gaussian bump peaking ~80 ms before
// the QRS, so that the inertial A-wave flow it drives has decelerated to zero
// just as the ventricle contracts (mitral closure at the QRS, T1 ~25 ms later);
// wrap-safe so its tail continues across tau = 0. It raises the LA elastance, so
// the a-wave and the A wave come out of the LA volume it holds.
function atrialAct(tau) {
  let dd = tau - (0.905); dd -= Math.round(dd);
  const d = dd / 0.045;
  return Math.exp(-d * d);
}

// Per-pathology circulation parameters. Normal values are calibrated to
// EDV≈121, ESV≈47, EF≈61%, aortic ~125/76, LVEDP≈8.
function circParams(path) {
  const p = {
    Emax: 2.7,   // end-systolic elastance (contractility)
    Emin: 0.085, // diastolic (passive) elastance / stiffness
    V0: 10,      // unstressed LV volume
    // Left atrium: a compliant chamber, Pla = Ela(t)·(Vla − V0la), filled from the
    // pulmonary veins (a venous pressure behind a small resistance Rpv) and
    // emptied through the mitral valve. Ela rises by aWave·ElaMin at atrial
    // systole (the booster pump), so the a-wave, the reservoir filling in systole,
    // the MR v-wave and the E-wave deceleration all come from the LA volume.
    // Pla0 is the chronic MEAN LA pressure: the pulmonary venous pressure behind
    // Rpv is set, cycle by cycle, to Pla0 + Rpv x (mean vein inflow), so the mean
    // LAP of the solved loop is Pla0 while its waveform (a, x, v, y) is free.
    Pla0: 8.0,
    Rpv: 0.040,  // pulmonary-vein -> LA resistance (phasic S/D vein flow, a still diastasis)
    ElaMin: 0.20, // passive LA elastance (stiff in MS, compliant in chronic MR)
    aWave: 1.1,  // atrial activation: Ela rises to ElaMin·(1 + aWave) at the A wave
    V0la: 10,    // unstressed LA volume
    // Annular descent: as the LV long axis shortens in systole the mitral annulus
    // is pulled toward the apex and the LA lengthens, so the atrium gains capacity
    // (its reservoir phase, and the systolic pulmonary-vein S wave). laDesc is the
    // capacity gained at a normal end-systole (annulus area x MAPSE, ~8 cm2 x 1.5 cm);
    // it follows the LV long-axis shortening, rho^0.24 as in the geometry.
    laDesc: 20,
    Rmv: 0.009,  // mitral (inflow) resistance
    // Mitral inertance (blood mass across the annulus): the E wave accelerates over
    // ~80-100 ms instead of jumping to its peak, overshoots, decelerates to a
    // near-still diastasis, then the A wave.
    Lmv: 0.0009,
    // Stenotic mitral ORIFICE (MS): a Bernoulli loss, dP = Kmv·Q² (= 4v² with
    // v = Q / area), on top of Rmv. With it the transmitral velocity falls
    // linearly as the LA empties (the "ski slope"), at a rate set by the orifice
    // area and the chamber elastances — the pressure half-time.
    Kmv: 0,
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
  if (path.dilated) {
    p.Emax = 0.60; p.V0 = 30; p.Rsys = 1.30; p.Cao = 1.4;
    // advanced DCM: a raised filling pressure into a stiffer, overfilled LV and a
    // failing, dilated atrium (weak booster) give the RESTRICTIVE pattern that goes
    // with an LVEF ~25 % and functional MR: a tall E wave that decelerates fast
    // (DT < 150 ms) and a small A (E/A > 2), not impaired relaxation.
    p.Pla0 = 18; p.Emin = 0.075;
    p.ElaMin = 0.35; p.V0la = 40; p.aWave = 0.2; p.laDesc = 8;
  }
  // DCM: annular dilatation and leaflet tethering leave a central FUNCTIONAL MR
  // (moderate: regurgitant fraction ~30 %) even with structurally normal leaflets
  if (path.dilated && !path.mr) { p.Rmr = 1.05; }
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
    // The LA takes the regurgitant volume in systole: a dilated chamber (V0la up)
    // whose pressure climbs to a v-wave (severe: >= 25 mmHg) that blunts the
    // late-systolic MR velocity, then empties fast into the LV (a tall E wave,
    // E-dominant inflow); its booster function wanes as it dilates.
    p.Rmr = ({ mild: 1.6, moderate: 0.85, severe: 0.33 })[g];
    p.Pla0 = ({ mild: 9, moderate: 12, severe: 17 })[g];
    p.V0 = ({ mild: 20, moderate: 34, severe: 40 })[g];
    p.Emax = Math.min(p.Emax, ({ mild: 2.7, moderate: 2.4, severe: 1.8 })[g]);
    p.ElaMin = ({ mild: 0.25, moderate: 0.32, severe: 0.4 })[g];
    p.V0la = ({ mild: 25, moderate: 38, severe: 50 })[g];
    p.aWave = ({ mild: 0.9, moderate: 0.7, severe: 0.5 })[g];
    p.laDesc = ({ mild: 18, moderate: 14, severe: 10 })[g];
  }
  if (path.mitralStenosis) {                                      // narrow inflow, high LAP
    // a Bernoulli orifice (effective area ~0.68 x the planimetered MS_AREA, the
    // contraction of the funnel-shaped jet) behind a dilated, STIFF atrium filled
    // through congested pulmonary veins: the LA rises to a v-wave in systole and
    // empties slowly through the whole of diastole, so the gradient (and the CW
    // velocity) decays along a ski slope whose pressure half-time follows the
    // area; the fibrosed atrium keeps only a weak booster (reduced A bump). The
    // rigid rheumatic annulus barely descends.
    const aEff = 0.684 * MS_AREA[g];
    p.Kmv = 1 / (2500 * aEff * aEff);
    p.Rmv = 0.005;
    p.Pla0 = ({ mild: 12, moderate: 16, severe: 20 })[g];
    p.ElaMin = ({ mild: 0.4, moderate: 0.6, severe: 0.8 })[g];
    p.V0la = ({ mild: 50, moderate: 60, severe: 70 })[g];
    p.Rpv = ({ mild: 0.2, moderate: 0.35, severe: 0.5 })[g];
    p.laDesc = ({ mild: 12, moderate: 8, severe: 5 })[g];
    p.aWave = 0.6;
  }
  if (path.rvpo) {
    // pulmonary hypertension: the failing right heart delivers less to the left
    // (low LA filling pressure) and the flattened, leftward-shifted septum stiffens
    // LV filling (ventricular interdependence) — the D-shaped LV is underfilled
    p.Pla0 = Math.min(p.Pla0, 7); p.Emin *= 1.3;            // LV EDV ~90 mL
  }
  if (path.rwma) { p.Emax = Math.min(p.Emax, 1.8); }              // regional dysfunction → low global contractility
  Object.assign(p, shuntParams(path));
  return p;
}

// ---- intracardiac shunts (self-contained) ---------------------------------
// Secundum ASD and VSD as Bernoulli orifices between the modelled LEFT heart and
// a prescribed RIGHT-heart pressure (the right heart is not lumped-modelled):
//   Q = A_eff x 50 sqrt(|dP|) mL/s  (v = sqrt(dP/4) m/s = 50 sqrt(dP) cm/s)
// ASD: LA -> RA through the cycle, driven by the small LA-RA difference (a few
// mmHg: ~1 m/s), largest in late systole (LA v wave) and early diastole, with a
// dip as both atria contract. VSD: LV -> RV, a fast jet in systole that follows
// the LV-RV gradient (RV systolic ~25 mmHg), a trickle in diastole. Either way
// the pulmonary circuit carries the shunt on top of the systemic output, so
// Qp:Qs = (forward SV + shunt volume) / forward SV; the LA keeps its chronic mean
// pressure (Pla0) because the pulmonary venous return rises with the shunt.
function shuntParams(path) {
  return {
    asdA: path.asd ? 0.62 : 0,     // effective ASD orifice (cm2): the ~1.4 x 0.7 cm oval defect x Cd
    vsdA: path.vsd ? 0.12 : 0,     // restrictive VSD vena contracta (cm2)
  };
}
// right-atrial pressure over the cycle (mmHg): mean ~3.5, a wave just BEFORE the
// left atrial one (the RA is activated first), the c wave as the tricuspid
// bulges back at the QRS, x descent in systole, v wave at
// end-systole, y descent once the tricuspid opens
function praAt(tau) {
  const g = (c, w) => { let d = tau - c; d -= Math.round(d); return Math.exp(-(d / w) * (d / w)); };
  return 3.5 + 4.2 * g(0.875, 0.045) + 2.4 * g(0.0, 0.035) - 1.6 * g(0.22, 0.1) + 1.6 * g(0.46, 0.06) - 1.4 * g(0.6, 0.07);
}
// right-ventricular pressure (mmHg): ~25/4, systole on the left heart's clock
function prvAt(tau, a) { return 4 + 22 * Math.min(1, a * 1.15); }
const orifice = (A, dP) => (A > 0 ? A * 50 * dP / Math.sqrt(Math.abs(dP) + 0.05) : 0);

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

// shunt-jet velocity (m/s, signed) from the pressure difference across it (Bernoulli)
const shuntV = (dP) => (dP >= 0 ? 1 : -1) * Math.sqrt(Math.abs(dP) / 4);

// Integrate the closed loop to its limit cycle and return a sampled trace of one
// cycle plus scalar summary (EDV/ESV/SV/EF/pressures/gradients/regurg fraction).
function simulate(path) {
  const p = circParams(path);
  const dt = CYCLE / STEPS;
  let Vlv = 120, Pao = 80, Qao = 0, Vla = 50, Qmv = 0;
  let edvPrev = 120, edvRun = 0;      // last cycle's EDV (annular-descent reference)
  let Ppv = p.Pla0 + p.Rpv * 90, qpvSum = 0;   // pulmonary venous pressure (see Pla0)
  const descNorm = 1 - Math.pow(0.4, 0.24);
  // LVOT/aortic inertance (mmHg·s²/mL) and aortic characteristic impedance (mmHg·s/mL)
  const L = 0.0015, Zc = 0.08;

  const tr = {
    phase: new Float64Array(NSAMPLE), Plv: new Float64Array(NSAMPLE),
    Pao: new Float64Array(NSAMPLE), Vlv: new Float64Array(NSAMPLE),
    Qmv: new Float64Array(NSAMPLE), Qao: new Float64Array(NSAMPLE),
    Qmr: new Float64Array(NSAMPLE), Pla: new Float64Array(NSAMPLE),
    Vla: new Float64Array(NSAMPLE), Qpv: new Float64Array(NSAMPLE),
    Qasd: new Float64Array(NSAMPLE), Qvsd: new Float64Array(NSAMPLE),
    dPasd: new Float64Array(NSAMPLE), dPvsd: new Float64Array(NSAMPLE),
  };

  for (let cyc = 0; cyc < NCYCLE; cyc++) {
    const record = cyc === NCYCLE - 1;
    let sIdx = 0;
    if (cyc > 0) { edvPrev = edvRun; Ppv = p.Pla0 + p.Rpv * qpvSum / STEPS; }
    edvRun = 0; qpvSum = 0;
    for (let s = 0; s < STEPS; s++) {
      const tau = s / STEPS;
      const a = activation(tau);
      const E = p.Emin + (p.Emax - p.Emin) * a;
      const Plv = E * (Vlv - p.V0);
      if (Vlv > edvRun) edvRun = Vlv;
      const desc = p.laDesc * Math.max(0, 1 - Math.pow(Math.min(1, Vlv / edvPrev), 0.24)) / descNorm;
      const Pla = p.ElaMin * (1 + p.aWave * atrialAct(tau)) * Math.max(0, Vla - p.V0la - desc);
      // mitral valve: a diode with resistance AND inertance (semi-implicit, like
      // the aortic outflow below); it opens when LA > LV and shuts as the
      // forward flow decelerates to zero
      if (Qmv > 0 || Pla > Plv) {
        Qmv = (Qmv + dt / p.Lmv * (Pla - Plv)) / (1 + dt * (p.Rmv + p.Kmv * Qmv) / p.Lmv);
        if (Qmv < 0) Qmv = 0;
      }
      // pulmonary-vein inflow (bidirectional: atrial systole pushes a little back)
      const Qpv = (Ppv - Pla) / p.Rpv;
      qpvSum += Qpv;
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
      const dPasd = p.asdA ? Pla - praAt(tau) : 0, dPvsd = p.vsdA ? Plv - prvAt(tau, a) : 0;
      const Qasd = orifice(p.asdA, dPasd), Qvsd = orifice(p.vsdA, dPvsd);
      Vlv += (Qmv - Qao - Qmr - Qvsd) * dt;
      Vla += (Qpv + Qmr - Qmv - Qasd) * dt;
      Pao += ((Qao - Qsys) / p.Cao) * dt;
      if (Vlv < 0) Vlv = 0;

      if (record) {
        // write this step into every sample slot it covers
        const target = Math.floor((s + 1) * NSAMPLE / STEPS);
        while (sIdx < target && sIdx < NSAMPLE) {
          tr.phase[sIdx] = sIdx / NSAMPLE;
          tr.Plv[sIdx] = Plv; tr.Pao[sIdx] = Pao + Zc * Qao; tr.Vlv[sIdx] = Vlv;   // Pao: ascending-aortic pressure
          tr.Qmv[sIdx] = Qmv; tr.Qao[sIdx] = Qao; tr.Qmr[sIdx] = Qmr; tr.Pla[sIdx] = Pla;
          tr.Vla[sIdx] = Vla; tr.Qpv[sIdx] = Qpv;
          tr.Qasd[sIdx] = Qasd; tr.Qvsd[sIdx] = Qvsd; tr.dPasd[sIdx] = dPasd; tr.dPvsd[sIdx] = dPvsd;
          sIdx++;
        }
      }
    }
  }

  // summary scalars from the recorded cycle
  let EDV = -1e9, ESV = 1e9, PlvSys = -1e9, PaoSys = -1e9, PaoDia = 1e9;
  let grad = 0, fwd = 0, regurg = 0, tMin = 0, PlaAtSys = 0, gradMV = 0, mvSum = 0, mvN = 0;
  let QaoMax = 0, QmvMax = 0, QmrMax = 0, iQmax = 0, iOpen = -1, iClose = -1, mrGradMax = 0;
  let shunt = 0, vAsdMax = 0, vAsdMin = 0, vVsdMax = 0;
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
    if (tr.Qmr[i] > 0.5 && tr.Plv[i] - tr.Pla[i] > mrGradMax) mrGradMax = tr.Plv[i] - tr.Pla[i];
    fwd += tr.Qao[i] * dtS;
    regurg += tr.Qmr[i] * dtS;
    shunt += (tr.Qasd[i] + tr.Qvsd[i]) * dtS;
    vAsdMax = Math.max(vAsdMax, shuntV(tr.dPasd[i])); vAsdMin = Math.min(vAsdMin, shuntV(tr.dPasd[i]));
    vVsdMax = Math.max(vVsdMax, shuntV(tr.dPvsd[i]));
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
    mrGradMax,              // peak LV-LA gradient while regurgitating (MR jet Vmax = sqrt(/4))
    accelTime: (iQmax - iOpen) * dtS, ejectTime: (iClose - iOpen + 1) * dtS,
    Pla0: p.Pla0,           // chronic LA pressure (drives LA remodelling in the anatomy)
    // intracardiac shunt (ASD / VSD): net left-to-right volume per beat, the
    // pulmonary-to-systemic flow ratio and the shunt-jet velocities (m/s, + = L->R)
    shuntVol: shunt, QpQs: fwd > 0 ? (fwd + shunt) / fwd : 1,
    vAsdPeak: vAsdMax, vAsdMin, vVsdPeak: vVsdMax,
  };
  Object.assign(sum, laMetrics(tr), inflowMetrics(tr));
  return { tr, sum, params: p };
}

// LA pressure and volume over the cycle: mean LAP, the systolic v-wave peak, the
// a-wave peak, and the reservoir (max), pre-atrial-contraction and minimum volumes
// with the passive (conduit) and active (booster) emptying fractions.
function laMetrics(tr) {
  let pSum = 0, vWave = 0, aWaveP = 0, vMax = -1e9, vMin = 1e9;
  for (let i = 0; i < NSAMPLE; i++) {
    const ph = tr.phase[i];
    pSum += tr.Pla[i];
    if (ph < 0.56 && tr.Pla[i] > vWave) vWave = tr.Pla[i];
    if ((ph >= 0.84 || ph < 0.03) && tr.Pla[i] > aWaveP) aWaveP = tr.Pla[i];
    if (tr.Vla[i] > vMax) vMax = tr.Vla[i];
    if (tr.Vla[i] < vMin) vMin = tr.Vla[i];
  }
  const vPreA = tr.Vla[Math.round(0.84 * NSAMPLE)];
  return {
    PlaMean: pSum / NSAMPLE, PlaVWave: vWave, PlaAWave: aWaveP,
    LAVmax: vMax, LAVmin: vMin, LAVpreA: vPreA,
    laPassiveEF: vMax > 0 ? (vMax - vPreA) / vMax : 0,
    laActiveEF: vPreA > 0 ? (vPreA - vMin) / vPreA : 0,
    laTotalEF: vMax > 0 ? (vMax - vMin) / vMax : 0,
  };
}

// Mitral-inflow (E/A) envelope metrics read off the modelled transmitral flow.
// The Doppler velocity is the flow scaled by a constant orifice area (and in MS
// the E-wave peak IS the Bernoulli peak gradient), so the timing and shape ratios
// below are what a PW/CW trace of the inflow shows:
//   eAccel  time from mitral opening to the E peak (s)
//   eDT     deceleration time: E peak to where the deceleration slope (fitted
//           between 85 % and 45 % of E) reaches baseline (s)
//   pht     pressure half-time: E peak to E/sqrt(2) on the velocity (s), the
//           early (steep) slope extrapolated if the A wave arrives first
//   diastasis  lowest flow between the E and A waves as a fraction of E
//   eaFlow  E / A peak-flow ratio; aPresent: an A-wave bump above the pre-A flow
function inflowMetrics(tr) {
  const dtS = CYCLE / NSAMPLE, Q = tr.Qmv;
  let iO = -1;
  for (let i = Math.round(0.3 * NSAMPLE); i < NSAMPLE; i++) if (Q[i] > 0.5) { iO = i; break; }
  if (iO < 0) return {};
  const iAon = Math.round(0.86 * NSAMPLE);
  let iE = iO;
  for (let i = iO; i < iAon; i++) { if (Q[i] > Q[iE]) iE = i; else if (Q[i] < 0.7 * Q[iE]) break; }
  const E = Q[iE];
  const cross = (frac) => {          // first time after the E peak the flow drops below frac·E
    for (let i = iE; i < iAon; i++) if (Q[i] <= frac * E) {
      const f = (Q[i - 1] - frac * E) / Math.max(1e-9, Q[i - 1] - Q[i]);
      return (i - 1 + f) * dtS;
    }
    return null;
  };
  const tE = iE * dtS, t85 = cross(0.85), t45 = cross(0.45);
  let eDT = null;
  if (t85 != null && t45 != null && t45 > t85) eDT = (t85 - tE) + 0.85 * (t45 - t85) / 0.40;
  else if (t85 != null) {             // slow (ski-slope) decay: slope from the peak to the pre-A flow
    const qA = Q[iAon - 1], tA = (iAon - 1) * dtS;
    eDT = E > qA ? (tA - tE) * E / (E - qA) : null;
  }
  let pht = cross(Math.SQRT1_2);
  if (pht != null) pht -= tE;
  else {
    const qA = Q[iAon - 1], tA = (iAon - 1) * dtS;
    pht = E > qA ? (tA - tE) * (1 - Math.SQRT1_2) * E / (E - qA) : null;
  }
  // diastasis: the flow nadir between the E and A waves (up to the A-wave upstroke)
  let dMin = Infinity, q75 = 0;
  for (let i = iE; i < Math.round(0.90 * NSAMPLE); i++) dMin = Math.min(dMin, Q[i]);
  for (let i = Math.round(0.75 * NSAMPLE); i <= Math.round(0.82 * NSAMPLE); i++) q75 = Math.max(q75, Q[i]);
  let A = 0;
  for (let i = iAon - 4; i < NSAMPLE + 8; i++) A = Math.max(A, Q[i % NSAMPLE]);
  const preA = Q[iAon - 6];
  return {
    mvoPhase: iO / NSAMPLE, eAccel: (iE - iO) * dtS, eDT, pht,
    diastasis: E > 0 ? dMin / E : 0, lateDiastolic: E > 0 ? q75 / E : 0, eaFlow: A > 0 ? E / A : null,
    aPresent: A > preA * 1.15 + 2,
  };
}

// ---- memoisation ---------------------------------------------------------
// The integration is heavy but depends only on the pathology flags, so cache by
// a signature of the mechanically-relevant flags (identical to how echo.js
// memoises its fixed-phase geometries). Bounded: one entry per distinct case.
const _cache = new Map();
function sig(path) {
  return [path.aorticStenosis, path.lvh, path.dilated, path.mr, path.mitralStenosis,
    path.rwma ? 1 : 0, path.rvpo, path.asd, path.vsd].map((v) => (v ? 1 : 0)).join('') + ':' + gradeOf(path);
}
function solved(path) {
  const key = sig(path);
  let r = _cache.get(key);
  if (!r) { r = simulate(path); _cache.set(key, r); }
  return r;
}

// Reusable singleton so the ~per-frame call in geometryAt allocates nothing.
const _state = {
  phase: 0, Plv: 0, Pao: 0, Vlv: 0, Pla: 0, Qmv: 0, Qao: 0, Qmr: 0, vAsd: 0, vVsd: 0, Qasd: 0, Qvsd: 0,
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
  _state.Qasd = lerp(tr.Qasd); _state.Qvsd = lerp(tr.Qvsd);
  _state.vAsd = shuntV(lerp(tr.dPasd)); _state.vVsd = shuntV(lerp(tr.dPvsd));
  _state.rho = sum.EDV > 0 ? clamp(_state.Vlv / sum.EDV, 0.05, 1) : 1;
  _state.kv = sum.SV > 0 ? clamp((sum.EDV - _state.Vlv) / sum.SV, 0, 1.2) : 0;
  _state.sum = sum;
  return _state;
}

// Scalar summary of the loop (EDV/ESV/SV/EF/pressures/gradient/regurg fraction).
export function hemoSummary(path = {}) { return solved(path).sum; }

// LA phasic function relative to the normal heart: c = conduit (passive) and
// b = booster (active) emptying fractions of this circulation over the normal
// ones. The anatomy scales its atrial reservoir/conduit/booster curve by these,
// so the drawn LA empties as this case's LA does (a stiff MS atrium barely
// moves; a failing DCM atrium loses its booster). {c: 1, b: 1} for a normal heart.
const _laCache = new Map();
export function laPhasic(path = {}) {
  const key = sig(path);
  let r = _laCache.get(key);
  if (!r) {
    const s = solved(path).sum, n = solved({}).sum;
    const rel = (a, b) => clamp(b > 0 ? a / b : 1, 0.15, 1.4);
    r = key === sig({}) ? { c: 1, b: 1 } : { c: rel(s.laPassiveEF, n.laPassiveEF), b: rel(s.laActiveEF, n.laActiveEF) };
    _laCache.set(key, r);
  }
  return r;
}

// Mitral-inflow E/A ratio from the modelled diastolic flow: peak early (E, rapid
// filling just after mitral opening, which the mitral inertance puts ~100 ms
// later) over peak late (A, atrial kick). Normal filling is E-dominant (~1-2); it
// falls below 1 with impaired relaxation and exceeds 2 when a high LA pressure
// meets a stiff LV (restrictive). Memoised.
let _eaPath = null, _eaVal;
export function eaRatio(path = {}) {
  if (path === _eaPath) return _eaVal;
  const { tr } = solved(path);
  let e = 0, a = 0;
  for (let i = 0; i < NSAMPLE; i++) {
    const ph = i / NSAMPLE, q = tr.Qmv[i];
    if (ph >= 0.36 && ph < 0.72) { if (q > e) e = q; }
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
