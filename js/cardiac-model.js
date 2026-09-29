// cardiac-model.js
// A dependency-free analytic model of the human heart used as the single source
// of truth for BOTH the 3D anatomy view and the 2D echo sector renderer.
//
// Everything lives in "heart space" (centimetres):
//   +y = towards the base (top),  -y = towards the apex (bottom)
//   +x = patient's left  (the LV / systemic side)
//   +z = anterior (towards the chest wall / transducer)
//
// Anatomy comes from anatomy.js: a landmark-based signed-distance heart at adult
// reference dimensions (LV, RV crescent + RVOT, atria on their annuli with a
// shared septum, aortic root / arch / descending aorta, PA, coronary sinus,
// pericardium, liver). This module adds the valves (thin hinged leaflets placed
// on those landmarks), the subvalvular apparatus and a piece-wise bulk-velocity
// blood-flow field, and couples everything to the lumped-parameter circulation.

import { clamp, smoothstep, pulse, unit } from './mathutils.js';
import { anatomyParams, atrialFill, bodyClassify, visceralDist, lumenDist, mitralLift, buildRwma, rwmaWeight, BODY, CFG, LM, BODY_AX } from './anatomy.js';
import { hemodynamics, hemoSummary, pvLoop, gradeOf, eaRatio, MS_AREA } from './hemodynamics.js';

// Re-export the lumped-parameter circulation API so echo.js / the UI can read the
// PV loop and instantaneous haemodynamics through cardiac-model (the single model
// entry point for the renderers) without importing hemodynamics.js directly.
export { hemodynamics, hemoSummary, pvLoop, gradeOf, eaRatio };

export const TISSUE = {
  OUTSIDE: 0,
  MYO: 1,        // myocardium (heart muscle)
  LV: 2,         // left-ventricular blood pool
  RV: 3,         // right-ventricular blood pool
  LA: 4,         // left atrium
  RA: 5,         // right atrium
  AORTA: 6,      // aortic root / ascending aorta lumen
  VALVE: 7,      // valve leaflet tissue
  PERICARDIUM: 8,// pericardial fluid (effusion)
  LUNG: 9,       // surrounding soft tissue / shadowing
  PERI_LINE: 10, // parietal pericardium / organ capsule (bright fibrous line)
  LIVER: 11,     // liver parenchyma (subcostal window)
  VEIN: 12,      // extracardiac venous blood (IVC, hepatic veins)
  FAT: 13,       // epicardial fat in the atrioventricular groove
  VWALL: 14,     // great-vessel wall
};

// Which flow compartment a velocity sample belongs to (for Doppler).
export const FLOW = {
  NONE: 0,
  MITRAL_IN: 1,   // LA -> LV  (diastolic inflow, towards apex)
  TRICUSPID_IN: 2,// RA -> RV
  LVOT: 3,        // LV -> aorta (systolic ejection)
  RVOT: 4,        // RV -> pulmonary artery
  MR_JET: 5,      // mitral regurgitation LV -> LA (systole)
  AS_JET: 6,      // aortic stenosis high-velocity jet
  VSD_JET: 7,     // ventricular septal defect LV -> RV
  ASD_JET: 8,     // atrial septal defect LA -> RA
  TR_JET: 9,      // tricuspid regurgitation RV -> RA
  PV_FLOW: 10,    // pulmonary-vein flow into the LA (reverses in severe MR)
};

// ---------------------------------------------------------------------------
// Baseline geometry (end-diastolic) — coarse ONE-ellipsoid-per-chamber proxies.
//
// SOURCE OF TRUTH: anatomy.js CFG (imported) is the authoritative geometry for
// TISSUE classification and the 2D/3D anatomy. BASE below is only the coarse
// ellipsoid proxy that the flow field, chamber labels and the numeric
// Measurements panel sample. The LV/LA/RA proxies are derived directly from CFG
// (so LVIDd/EF and Doppler gating match the visual anatomy — no dual-source
// drift), while the RV and aorta are single ellipsoids that bound the multi-lobe
// RV crescent / aortic root for cheap "am I inside this chamber?" flow gating.
// Radii are semi-axes of each ellipsoid.
// ---------------------------------------------------------------------------
// The LV proxy is the prolate ellipsoid with the cavity's length and equatorial
// radius (so 2·r[0] is LVIDd and the apex-anchored top tracks the annulus: MAPSE);
// the RV proxy bounds the crescent body. The atrial proxies are rebuilt per
// phase from the live anatomy (they sit on their moving annuli).
const BASE = {
  lv: { c: CFG.lv.c, r: CFG.lv.r, wall: CFG.lv.wall },   // LVIDd = 2*r[0] = 4.9 cm
  rv: { c: CFG.rvBody.c, r: CFG.rvBody.r, wall: CFG.rvWall },
  // aortic-root proxy (sinuses) for LVOT / aorta flow gating and the "Ao" label
  aorta: { c: [LM.A[0] + LM.U_AO[0] * 0.9, LM.A[1] + LM.U_AO[1] * 0.9, LM.A[2] + LM.U_AO[2] * 0.9], r: [1.6, 1.6, 1.6] },
};

// Valve definitions (end-diastolic): annulus centre + normal (pointing
// downstream->upstream for the AV valves, i.e. toward the atrium; along the
// outflow for the semilunar valves), annulus radius and when it opens. All come
// from the anatomy landmarks. The tricuspid annulus sits 0.85 cm APICAL to the
// mitral — the normal septal-leaflet offset that identifies the RV on A4C.
const _tvN0 = unit([LM.T[0] + 2.3, LM.T[1] - (CFG.rvBody.c[1] - CFG.rvBody.r[1] + 0.6), LM.T[2] - 1.9]);
const VALVES = {
  mitral:    { c: LM.M, n: unit([0.02, 1, -0.08]), r: LM.MV_R, opensInDiastole: true },
  tricuspid: { c: LM.T, n: _tvN0, r: LM.TV_R, opensInDiastole: true },
  aortic:    { c: LM.A, n: LM.U_AO, r: LM.AO.annR, opensInDiastole: false },
  pulmonic:  { c: LM.PV, n: LM.U_PA, r: LM.PV_R, opensInDiastole: false },
};

// Cross-file contract: expose the valve plane definitions (centre c, normal n,
// annulus radius r, opening phase) so heart3d.js can import them instead of
// duplicating the numbers. Do NOT mutate — treat as read-only.
export const VALVE_DEFS = VALVES;

// ---------------------------------------------------------------------------
// Leaflet geometry. Real AV/semilunar leaflets are thin curtains that hinge at
// the annulus and swing as flaps — NOT a flat disc that irises open. Each valve
// is modelled as a funnel of leaflet shells that extend DOWNSTREAM (`flow`) from
// the annulus: at the hinge (axial 0) the shell sits at the annulus radius; at
// the coaptation depth `cd` it tapers to the free-edge radius `open·r`. So when
// the valve shuts (open→0) the free edges meet on the axis and the leaflets tent
// into the downstream chamber (the mitral systolic dome / aortic diastolic
// closure line); when it opens they swing back toward the walls, leaving a
// central orifice of radius open·r (the diastolic "parallel leaflets" / systolic
// aortic "box"). The mitral is bileaflet-asymmetric: the anterior leaflet (facing
// the aortomitral curtain / LVOT) is longer and drapes deeper than the posterior.
function _perp(v, axis) {                     // component of v perpendicular to unit `axis`
  const d = v[0] * axis[0] + v[1] * axis[1] + v[2] * axis[2];
  return unit([v[0] - d * axis[0], v[1] - d * axis[1], v[2] - d * axis[2]]);
}
function _cross(a, b) {
  return unit([a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]);
}
// Orthonormal in-plane basis (u,w) perpendicular to a unit axis `f` — used to give
// the semilunar valves an angular coordinate for their three cusps/commissures.
function _basis(f) {
  const ref = Math.abs(f[0]) < 0.9 ? [1, 0, 0] : [0, 0, 1];
  const u = _perp(ref, f);
  return { u, w: _cross(f, u) };
}
const _mFlow = unit([LM.apex[0] - LM.M[0], LM.apex[1] - LM.M[1], LM.apex[2] - LM.M[2]]); // mitral: into the LV, at the apex
const _tvFlow = [-_tvN0[0], -_tvN0[1], -_tvN0[2]]; // tricuspid: into the RV, at the RV apex
const _mAP = _perp([VALVES.aortic.c[0] - VALVES.mitral.c[0],
                    VALVES.aortic.c[1] - VALVES.mitral.c[1],
                    VALVES.aortic.c[2] - VALVES.mitral.c[2]], _mFlow);
// tricuspid septal-leaflet axis: toward the interventricular septum / mitral
const _tvAP = _perp([VALVES.mitral.c[0] - VALVES.tricuspid.c[0],
                     VALVES.mitral.c[1] - VALVES.tricuspid.c[1],
                     VALVES.mitral.c[2] - VALVES.tricuspid.c[2]], _tvFlow);
const _aoFlow = LM.U_AO;                      // aortic: up the root
const _puFlow = LM.U_PA;                      // pulmonic: up the PA
const _aoB = _basis(_aoFlow), _puB = _basis(_puFlow);
const LEAFLETS = {
  // AV valves (mitral/tricuspid): flow = downstream leaflet-extension axis; ap =
  // in-plane axis toward the longer (anterior) leaflet; com = intercommissural axis
  // (posterior scallops); cd = coaptation depth; asym = anterior/posterior ratio.
  mitral:    { flow: _mFlow, cd: 1.75, thick: 0.075, thickP: 0.09,
               ap: _mAP, com: _cross(_mFlow, _mAP), scallops: 3, asym: 1.55 },
  // tricuspid: the septal leaflet (toward the IVS, +ap) is the shortest and least
  // mobile; the anterior/lateral leaflet is longer — so asym < 1 shortens +ap.
  tricuspid: { flow: _tvFlow, cd: 1.35, thick: 0.08, ap: _tvAP, com: _cross(_tvFlow, _tvAP), asym: 0.7 },
  // Semilunar valves (aortic/pulmonic): three cusps. u/w give the angular frame;
  // cusps=3 triggers the trilobed (triangular) orifice + commissural coaptation
  // seams that read as the short-axis "Mercedes" Y when the cusps shut.
  // aortic cusps share the anatomy's sinus frame (PSAX-AV screen-right / anterior):
  // cusp centres at 90 (right-coronary, anterior), 210 (non-coronary) and 330 deg
  // (left-coronary), so each cusp sits in its own sinus of Valsalva and the
  // commissures (30 / 150 / 270 deg) fall between the sinuses
  aortic:    { flow: _aoFlow, cd: 0.95, thick: 0.065, ap: null, asym: 1,
               cusps: 3, u: LM.E_SCR, w: LM.E_ANT, commOff: Math.PI / 6 },
  pulmonic:  { flow: _puFlow, cd: 0.9,  thick: 0.065, ap: null, asym: 1,
               cusps: 3, u: _puB.u, w: _puB.w, commOff: 0.0 },
};

// The three aortic cusps' identities + centroid anchors (in heart space), for the
// aortic-valve short-axis labels. Identity is assigned by ANATOMY, not a fixed
// order: the right-coronary cusp faces the RVOT/pulmonary valve (anterior), the
// non-coronary cusp faces the interatrial septum / right heart, and the
// left-coronary cusp is the remaining (leftward, toward the LA). This stays
// correct regardless of the valve's angular offset.
export const AORTIC_CUSPS = (() => {
  const v = VALVES.aortic, lf = LEAFLETS.aortic, rc = 0.6;
  const cusps = [0, 1, 2].map((k) => {
    const phi = lf.commOff + Math.PI / 3 + k * 2 * Math.PI / 3;
    const cs = Math.cos(phi), sn = Math.sin(phi);
    return { p: [
      v.c[0] + (lf.u[0] * cs + lf.w[0] * sn) * rc,
      v.c[1] + (lf.u[1] * cs + lf.w[1] * sn) * rc,
      v.c[2] + (lf.u[2] * cs + lf.w[2] * sn) * rc] };
  });
  const dot = (p, ref) => (p[0] - v.c[0]) * (ref[0] - v.c[0]) +
    (p[1] - v.c[1]) * (ref[1] - v.c[1]) + (p[2] - v.c[2]) * (ref[2] - v.c[2]);
  // RCC: cusp facing the RVOT/pulmonary valve (anterior)
  let rcc = 0; for (let k = 1; k < 3; k++) if (dot(cusps[k].p, VALVES.pulmonic.c) > dot(cusps[rcc].p, VALVES.pulmonic.c)) rcc = k;
  const rest = [0, 1, 2].filter((k) => k !== rcc);
  // NCC: of the remaining two, the one facing the tricuspid/right heart (IAS)
  const ncc = dot(cusps[rest[0]].p, VALVES.tricuspid.c) > dot(cusps[rest[1]].p, VALVES.tricuspid.c) ? rest[0] : rest[1];
  const lcc = rest.find((k) => k !== ncc);
  cusps[rcc].name = 'RCC'; cusps[ncc].name = 'NCC'; cusps[lcc].name = 'LCC';
  return cusps;
})();

// ---------------------------------------------------------------------------
// Cardiac timing. phase in [0,1); 0 == onset of systole (QRS / mitral close).
// ---------------------------------------------------------------------------
// Atrial contraction (atrial kick) late in diastole.
function atrialKick(phase) {
  return pulse(0.86, 1.0, phase);
}

// Right-heart timing. The right ventricle is not separately lumped-modelled, so it
// replays the left-ventricular trace on its own clock instead of sharing the
// left-heart valve events: the tricuspid closes ~25 ms after the mitral (T1 after
// M1), the pulmonic valve opens ~10 ms before the aortic (shorter RV isovolumic
// contraction against the low pulmonary pressure) and closes ~30 ms after it (P2
// after A2: the longer RV ejection into the compliant pulmonary bed), and the
// tricuspid opens ~10 ms before the mitral (RV IVRT ~35 ms). Each row maps an LV
// reference phase (where the right-heart event happens) onto the LV-trace phase
// whose state the right heart shows then. At the reference 0.833 s cycle: 0.030
// = 25 ms, 0.012 = 10 ms, 0.036 = 30 ms.
export const RV_EVENTS = [
  [0.026, -0.004],   // TVC  <- MVC
  [0.060, 0.072],    // PVO  <- AVO
  [0.457, 0.421],    // PVC (P2) <- AVC (A2)
  [0.501, 0.513],    // TVO  <- MVO
];
export function rvPhase(phase) {
  let p = ((phase % 1) + 1) % 1;
  if (p < RV_EVENTS[0][0]) p += 1;
  for (let i = 0; i < RV_EVENTS.length; i++) {
    const a = RV_EVENTS[i], b = i + 1 < RV_EVENTS.length ? RV_EVENTS[i + 1] : [RV_EVENTS[0][0] + 1, RV_EVENTS[0][1] + 1];
    if (p < b[0]) { const q = a[1] + (b[1] - a[1]) * (p - a[0]) / (b[0] - a[0]); return ((q % 1) + 1) % 1; }
  }
  return phase;
}

// Valve opening 0..1, driven by the modelled transvalvular FLOW so the leaflets
// open exactly when blood starts to cross and shut when it stops (the valve
// events — AVO, AVC, MVO, MVC — and the isovolumic periods all come from the
// circulation). An AV valve swings wide on the E wave, drifts half-shut in
// diastasis and reopens on the A wave; a semilunar valve opens fully in early
// ejection and floats toward closure as flow decelerates. The right-heart valves
// are passed the right heart's own snapshot (rvPhase), so they open and close
// on its clock.
function valveOpening(name, hd, path) {
  const s = hd.sum;
  let o;
  if (VALVES[name].opensInDiastole) {
    const q = hd.Qmv;
    o = q > 1 ? 0.3 + 0.7 * clamp(q / (0.55 * s.QmvMax), 0, 1) : 0.02;
  } else {
    const q = hd.Qao;
    o = q > 1 ? 0.35 + 0.65 * clamp(q / (0.35 * s.QaoMax), 0, 1) : 0.02;
  }
  // pathology tweaks
  if (path) {
    if (name === 'aortic' && path.aorticStenosis) o *= 0.28;   // restricted opening
    if (name === 'mitral' && path.mitralStenosis) o *= 0.35;
  }
  return clamp(o, 0, 1);
}

// ---------------------------------------------------------------------------
// Effective orifice areas (cm^2) and phase durations used to turn the modelled
// stroke volume / flows into physiologic jet PEAK velocities via the continuity
// equation (v = Q / A). Making the Doppler velocities fall out of the modelled
// SV and orifice area (rather than hand-set constants) is what makes severity
// track physiology: AS = same SV forced through a small AVA -> high velocity.
// ---------------------------------------------------------------------------
const AVA_NORM = 3.0, AVA_AS = 0.6;   // aortic (LVOT/valve) effective orifice area
const MVA_NORM = 5.0;                 // mitral effective orifice area (MS: hemodynamics.js MS_AREA)
const ET_EJECT = 0.33, T_FILL = 0.37; // systolic ejection / diastolic filling time (s)
const K_AO = 1.56, K_MV = 2.2;        // peak/mean profile factors (peak ≈ k·mean)

// Build the per-frame haemodynamic snapshot the flow field reads. Returns a
// FRESH object (never the hemodynamics() singleton) so two geometryAt() results
// — e.g. echo's end-diastolic and end-systolic geometries — never alias.
function hemoFrame(hd, path, rvHd) {
  const s = hd.sum;
  // effective orifice area follows the severity grade (continuity: same SV through
  // a smaller area → higher velocity), so the Doppler peak and BSE band both track it.
  const g = gradeOf(path);
  const AVA = path.aorticStenosis ? ({ mild: 0.88, moderate: 0.71, severe: 0.53 })[g] : AVA_NORM;
  // peak jet speeds (m/s): continuity for forward flows, Bernoulli (ΔP=4v²) for MR
  // AS: the Doppler peak IS the modelled peak LV->Ao gradient (Bernoulli), so the
  // spectral trace, the reported gradient and the severity grade cannot disagree
  const vAoPeak = path.aorticStenosis ? Math.sqrt(Math.max(0, s.gradient) / 4)
    : K_AO * s.forwardSV / (AVA * ET_EJECT) / 100;
  // MS: like AS, the Doppler E-wave peak IS the modelled peak LA-LV gradient
  // (Bernoulli), so the trace, the mean gradient, the LAP and the planimetered
  // MVA (MS_AREA) all describe the same valve
  const vMitPeak = path.mitralStenosis ? Math.sqrt(Math.max(0, s.gradientMV) / 4)
    : K_MV * s.forwardSV / (MVA_NORM * T_FILL) / 100;
  const vMRPeak = Math.sqrt(Math.max(0, s.PlvSys - s.PlaAtSys) / 4);
  // TR jet peak (m/s) encodes RV systolic pressure: high with pressure overload
  // (pulmonary hypertension), moderate for isolated TR. Feeds the PASP estimate.
  const vTRPeak = path.tr ? (path.rvpo ? 4.3 : 3.0) : 0;
  // instantaneous 0..1 envelopes from the modelled flow waveform (carry E/A,
  // ejection and regurgitant timing) so velocityAt only multiplies scalars.
  const ejAo = s.QaoMax > 0 ? clamp(hd.Qao / s.QaoMax, 0, 1) : 0;
  const ejMv = s.QmvMax > 0 ? clamp(hd.Qmv / s.QmvMax, 0, 1) : 0;
  const ejMr = s.QmrMax > 0 ? clamp(hd.Qmr / s.QmrMax, 0, 1) : 0;
  // right-heart envelopes (tricuspid inflow, RV outflow) on the right heart's clock
  const r = rvHd || hd;
  const ejTv = s.QmvMax > 0 ? clamp(r.Qmv / s.QmvMax, 0, 1) : 0;
  const ejPv = s.QaoMax > 0 ? clamp(r.Qao / s.QaoMax, 0, 1) : 0;
  return {
    Plv: hd.Plv, Pao: hd.Pao, Vlv: hd.Vlv, Pla: hd.Pla,
    Qao: hd.Qao, Qmv: hd.Qmv, Qmr: hd.Qmr, rho: hd.rho, kv: hd.kv,
    EDV: s.EDV, ESV: s.ESV, SV: s.SV, EF: s.EF, forwardSV: s.forwardSV,
    regurgFraction: s.regurgFraction, gradient: s.gradient,
    PlvSys: s.PlvSys, PaoSys: s.PaoSys, PaoDia: s.PaoDia,
    vAoPeak, vMitPeak, vMRPeak, vTRPeak, ejAo, ejMv, ejMr, ejTv, ejPv,
  };
}

// ---------------------------------------------------------------------------
// Live geometry for a given phase + pathology. Returns per-chamber ellipsoids
// (centre + radii) plus wall thicknesses, so both renderers stay consistent.
// ---------------------------------------------------------------------------
// Reference (normal) end-diastolic volume, for scaling remodelled ventricles.
let _edvRef = null;
const edvRef = () => (_edvRef != null ? _edvRef : (_edvRef = hemoSummary({}).EDV));

export function geometryAt(phase, path = {}) {
  const kick = atrialKick(phase);

  // PHYSICS COUPLING: the lumped-parameter circulation solves the LV pressure–
  // volume loop. hd.rho = Vlv/EDV is the ABSOLUTE cavity-volume ratio (1 at
  // end-diastole, ESV/EDV at end-systole), so it carries the true ejection
  // fraction (low in DCM/ischemia). The cavity semi-axes are scaled so their
  // product tracks rho EXACTLY (sShort²·sLong = rho), with the long axis
  // shortening less than the short axis; hence the measured EF/LVIDs fall out of
  // the modelled volume rather than a prescribed phase curve. hd.kv (normalised)
  // still times the wall-thickening. The right heart reads the same trace on its
  // own clock (rvPhase); copied out first, as hemodynamics() returns a singleton.
  const phaseRV = rvPhase(phase);
  const hdR = hemodynamics(phaseRV, path);
  const rvHd = { Qmv: hdR.Qmv, Qao: hdR.Qao, kv: hdR.kv, sum: hdR.sum };
  const kvR = rvHd.kv, kR = clamp(kvR, 0, 1);
  const hd = hemodynamics(phase, path);
  const rho = hd.rho;                 // Vlv / EDV (absolute cavity-volume ratio)
  const kv = hd.kv;                   // normalised contraction fraction (timing)
  // wall-thickening timing follows the modelled LV volume, so the walls, valves
  // and flows all share the circulation's event times (the RV: kvR / kR above)
  const k = clamp(kv, 0, 1);

  // Remodelling from the circulation: the cavity is sized to the modelled EDV
  // relative to normal, dilating more in the short axis than the long (a
  // dilated ventricle becomes more spherical: radial^2 * long = volume ratio);
  // the LA enlarges with the chronic LA pressure (MR, MS, DCM).
  const vr = hd.sum.EDV / edvRef();
  const lvScale = Math.pow(vr, 0.4), lvLong = Math.pow(vr, 0.2);
  const laScale = 1 + 0.045 * Math.max(0, (hd.sum.Pla0 || 8) - 8);
  let lvWallMul = 1.0;
  if (path.dilated) { lvWallMul = 0.82; }                    // thinned walls
  if (path.lvh || path.aorticStenosis) { lvWallMul = 1.55; } // concentric hypertrophy (IVSd ~1.4)

  // Long axis shortens less than the short axis. aLong = 0.24 gives, for the
  // normal EF, MAPSE ~1.55 cm, GLS ~-19 % and fractional shortening ~28 % on
  // the 8.3 cm ventricle — all mid-normal.
  const aLong = 0.24;
  const sLong = Math.pow(rho, aLong);
  const sShort = Math.pow(rho, (1 - aLong) / 2); // sShort²·sLong = rho (volume-exact)
  const lvR = [
    BASE.lv.r[0] * lvScale * sShort,
    BASE.lv.r[1] * lvLong * sLong,
    BASE.lv.r[2] * lvScale * sShort,
  ];
  // apex-anchored longitudinal contraction: the LV centre shifts apically as the
  // long axis shortens (apex fixed), so the proxy's top tracks the mitral annulus
  // and descends the full shortening (MAPSE). Volume-preserving.
  const lvR1d = BASE.lv.r[1] * lvLong;          // end-diastolic long semi-axis
  const lvApexY = BASE.lv.c[1] - BASE.lv.r[1] - (lvLong - 1) * 2 * BASE.lv.r[1]; // base fixed, apex moves
  const lvCy = lvApexY + lvR1d * sLong;         // anchored centre y
  // wall thickens as the short axis shrinks (approx. muscle-volume conservation)
  const lvWall = BASE.lv.wall * lvWallMul * clamp(1 / sShort, 1, 1.7);

  const rvR = BASE.rv.r.map((r, i) => {
    const shorten = i === 1 ? 0.12 : 0.26;
    return r * (1 - kvR * shorten);
  });
  const rvWall = BASE.rv.wall * (1 + kR * 0.5);

  // anatomical SDF bundle (built once here so the chordae apparatus below can
  // hang off the same apex-anchored papillary tips + live valve planes).
  const A = anatomyParams(k, kick, path, { sShort, sLong, lvWall, lvScale, lvLong, laScale, kRV: kR, phaseRV }, phase);

  // atria: coarse axis-aligned proxies of the live anatomical atria (which sit on
  // their moving annuli and follow the reservoir/conduit/booster curve), so flow
  // gating, labels and measurements stay locked to the rendered atrial wall.
  const laR = [A.la.r1, A.la.rl, A.la.r2];
  const raR = [A.ra.r1, A.ra.rl, A.ra.r2];
  const valves = {
    mitral: valveOpening('mitral', hd, path),
    tricuspid: valveOpening('tricuspid', rvHd, path),
    aortic: valveOpening('aortic', hd, path),
    pulmonic: valveOpening('pulmonic', rvHd, path),
  };

  return {
    phase, k, kick, kv, phaseRV, kRV: kR,
    // instantaneous lumped-parameter circulation state (PV loop, pressures,
    // stroke volume, EF, jet peak velocities) — read by the echo measurements,
    // the flow field and the UI PV-loop panel. A fresh snapshot (not the
    // hemodynamics() singleton) so distinct geometries never alias.
    hemo: hemoFrame(hd, path, rvHd),
    // Coarse per-chamber ellipsoid proxies — used by the velocity field, the
    // on-image labels and the LVIDd/EF measurements. Tissue classification and
    // the 3D surface use the anatomical SDF (G.A) below, not these.
    lv: { c: [BASE.lv.c[0], lvCy, BASE.lv.c[2]], r: lvR, wall: lvWall },
    rv: { c: BASE.rv.c, r: rvR, wall: rvWall },
    la: { c: A.la.c, r: laR },
    ra: { c: A.ra.c, r: raR },
    aorta: BASE.aorta,
    // anatomical SDF parameter bundle — classify() and the 3D marching-cubes
    // surface both sample this for a topologically correct heart. The same
    // volume-exact LV scaling (sShort/sLong) is threaded through so the SDF
    // cavity, the coarse proxy above and the measured EF all track the PV loop.
    A,
    // subvalvular apparatus: chordae tendineae fanning from the papillary tips
    // to the live leaflet free edges (so they go slack / taut with the leaflets).
    chordae: buildChordae({ A, valves, path }),
    valves,
    path,
  };
}

// ---------------------------------------------------------------------------
// Chordae tendineae. Thin fibrous strands from each papillary-muscle tip to the
// mitral leaflet free edges — echo shows them as fine linear echoes tethering
// the leaflets into the LV cavity (and, tellingly, they go SLACK in diastole and
// pull TAUT in systole). Modelled as a small fan of capsule segments so the
// tissue sampler can paint them over the LV blood pool. Built once per frame off
// the apex-anchored papillary tips (A.pap[i].b) and the mitral coaptation zone.
function buildChordae(G) {
  const A = G.A, path = G.path, segs = [];
  // Each papillary muscle tethers the halves of BOTH leaflets on its own side of
  // the valve (anterolateral PM -> the anterolateral halves, posteromedial PM ->
  // the posteromedial halves), inserting along the free edge.
  const tipsFor = (name, pb, sides) => {
    const lf = LEAFLETS[name], vc = A.valves[name].c;
    const s = Math.sign((pb[0] - vc[0]) * lf.com[0] + (pb[1] - vc[1]) * lf.com[1] + (pb[2] - vc[2]) * lf.com[2]) || 1;
    const out = [];
    for (const sgn of sides) for (const q of [0.3, 0.72]) out.push(avTipWorld(name, G, sgn, s * q));
    return out;
  };
  // mitral: both LV papillary tips -> anterior/posterior leaflet free edges
  const rMit = path.mitralStenosis ? 0.07 : 0.035;  // fine strands; thickened, fused in MS
  for (const p of A.pap) for (const q of tipsFor('mitral', p.b, [1, -1])) segs.push({ a: p.b, b: q, r: rMit });
  // tricuspid: the RV anterior papillary (fused to the moderator band) -> the
  // anterior/posterior leaflet free edge, completing the septum→band→PM→leaflet chain
  if (A.rvPap) {
    for (const q of tipsFor('tricuspid', A.rvPap.b, [-1])) segs.push({ a: A.rvPap.b, b: q, r: 0.033, tv: true });
  }
  return segs;
}

// Global longitudinal strain (GLS): peak systolic shortening of the LV long axis,
// as a percentage (negative = shortening). Derived from the modelled long-axis
// length over the cycle — the same quantity speckle-tracking reports. Normal is
// more negative than about −18 %; it falls toward zero as contractility drops
// (DCM) and drops regionally in RWMA. Memoised per pathology object.
let _glsPath = null, _glsVal = null;
export function longitudinalStrain(path = {}) {
  if (path === _glsPath && _glsVal !== null) return _glsVal;
  const N = 120;
  let Lmax = -1e9, Lmin = 1e9;
  for (let i = 0; i < N; i++) {
    const L = 2 * geometryAt(i / N, path).lv.r[1]; // apex→base long-axis length proxy
    if (L > Lmax) Lmax = L;
    if (L < Lmin) Lmin = L;
  }
  _glsVal = Lmax > 0 ? (Lmin - Lmax) / Lmax * 100 : null;
  _glsPath = path;
  return _glsVal;
}

// ---------------------------------------------------------------------------
// Regional longitudinal strain — AHA 17-segment "bullseye" (speckle-tracking).
// Peak systolic longitudinal strain per segment (%, negative = shortening).
// Normal is ~-20% and roughly uniform; global contractility scales every segment
// (DCM low, LVH mildly low), and a regional wall-motion abnormality drives the
// affected wall's segments toward zero. Segment angles use the model convention
// (lateral +x = 0, anterior +z = +pi/2, septal -x = pi, inferior -z = -pi/2), the
// same frame as the RWMA blend, so the cold segment matches the 2-D akinetic wall.
// ---------------------------------------------------------------------------
const NORMAL_SEG_STRAIN = -20;
// [name, angle, level] level: 0 basal, 1 mid, 2 apical, 3 apex. Six walls at 60°.
const AHA_WALLS6 = [
  ['anterior', Math.PI / 2], ['anteroseptal', Math.PI * 5 / 6], ['inferoseptal', -Math.PI * 5 / 6],
  ['inferior', -Math.PI / 2], ['inferolateral', -Math.PI / 6], ['anterolateral', Math.PI / 6],
];
const AHA_WALLS4 = [ // apical level: 4 walls
  ['anterior', Math.PI / 2], ['septal', Math.PI], ['inferior', -Math.PI / 2], ['lateral', 0],
];
const AHA_SEGMENTS = [
  ...AHA_WALLS6.map(([n, a]) => ({ name: 'basal ' + n, angle: a, level: 0 })),
  ...AHA_WALLS6.map(([n, a]) => ({ name: 'mid ' + n, angle: a, level: 1 })),
  ...AHA_WALLS4.map(([n, a]) => ({ name: 'apical ' + n, angle: a, level: 2 })),
  { name: 'apex', angle: 0, level: 3 },
];
// how "affected" (0..1) a segment is by the RWMA spec (anatomy.js buildRwma):
// the same territory weights the 2-D wall blend uses, so the cold segment sits
// where the akinetic wall is
function segAffected(seg, rw) {
  if (!rw) return 0;
  if (rw.apical) {
    return rw.sev * (seg.level >= 2 ? 1 : seg.level === 1 ? 0.35 : 0.05);
  }
  if (rw.lad) {                                        // longitudinal LAD territory
    const f = [0, 0.45, 1, 1][seg.level], cap = [0, 0, 0.7, 1][seg.level];
    return rw.sev * rwmaWeight(rw, seg.angle, f, cap);
  }
  if (seg.level === 3) return rw.sev * 0.3;            // apex: partial
  const levelF = seg.level === 2 ? 0.8 : 1;            // apical territory slightly less
  return rw.sev * rwmaWeight(rw, seg.angle, 0, 0) * levelF;
}

// Returns { segments: [{name, angle, level, strain}], gls } for a pathology.
let _rsPath = null, _rsVal = null;
export function regionalStrain(path = {}) {
  if (path === _rsPath && _rsVal) return _rsVal;
  let gf = 1;                                           // global contractility factor
  if (path.dilated) gf = 0.42;
  else if (path.aorticStenosis && path.lvh) gf = 0.85;
  else if (path.lvh) gf = 0.9;
  const rw = buildRwma(path);
  const segments = AHA_SEGMENTS.map((s) => {
    const w = segAffected(s, rw);
    return { name: s.name, angle: s.angle, level: s.level, strain: NORMAL_SEG_STRAIN * gf * (1 - w) };
  });
  const gls = segments.reduce((a, s) => a + s.strain, 0) / segments.length;
  _rsVal = { segments, gls };
  _rsPath = path;
  return _rsVal;
}

// squared normalised distance from ellipsoid centre (inside if <= 1)
function ellip(px, py, pz, c, r) {
  const dx = (px - c[0]) / r[0];
  const dy = (py - c[1]) / r[1];
  const dz = (pz - c[2]) / r[2];
  return dx * dx + dy * dy + dz * dz;
}

// Module-scoped singleton result for classify()'s non-valve returns. Callers
// read {tissue, echo} immediately, so returning the same object every call is
// safe and avoids a per-sample allocation. NOTE: valveTissue() writes into a
// SEPARATE singleton (_valveResult) so the two paths can never alias.
const _clsResult = { tissue: TISSUE.OUTSIDE, echo: 0 };
function cls(tissue, echo) {
  _clsResult.tissue = tissue;
  _clsResult.echo = echo;
  return _clsResult;
}

// ---------------------------------------------------------------------------
// classify(): the workhorse. Given a point + a precomputed geometry, decide
// which structure occupies it and return an echogenicity (0..1 brightness).
// ---------------------------------------------------------------------------
const BODY_TO_TISSUE = {
  [BODY.MYO]: TISSUE.MYO,
  [BODY.LV]: TISSUE.LV,
  [BODY.RV]: TISSUE.RV,
  [BODY.LA]: TISSUE.LA,
  [BODY.RA]: TISSUE.RA,
  [BODY.AO]: TISSUE.AORTA,
  [BODY.PAP]: TISSUE.MYO,
  [BODY.OUT]: TISSUE.OUTSIDE,
  [BODY.PERI]: TISSUE.PERI_LINE,
  [BODY.LIVER]: TISSUE.LIVER,
  [BODY.VESSELWALL]: TISSUE.VWALL,
  [BODY.FAT]: TISSUE.FAT,
  [BODY.VEIN]: TISSUE.VEIN,
  [BODY.LUNG]: TISSUE.LUNG,
};

// Pericardial effusion depth (cm) at a point: fluid is DEPENDENT, so it pools
// posteriorly (the patient lies supine / left-lateral for echo) and is thin
// anteriorly. Distance is measured from the LV centre along the body-posterior
// axis. The circumferential maximum (~1.8 cm posteriorly) is a moderate-large
// effusion; anteriorly it tapers to ~0.3 cm.
const _P = BODY_AX.P;
function effusionGap(px, py, pz, A) {
  const c = A.lv.c;
  const dep = (px - c[0]) * _P[0] + (py - c[1]) * _P[1] + (pz - c[2]) * _P[2];
  return 0.3 + 1.5 * smoothstep(-3.5, 2.5, dep);
}

export function classify(px, py, pz, G, path = {}) {
  const A = G.A;

  // --- pericardial effusion: a DEPENDENT echo-free layer between the
  // epicardium and the parietal pericardium (which it pushes outward) — so it
  // stops at the pericardial reflections, sits ANTERIOR to the descending aorta,
  // and is widest posteriorly. ---
  let periOff = 0;
  if (path.effusion) {
    const gap = effusionGap(px, py, pz, A);
    // taper to nothing toward the reflections at the great-vessel roots (no flat cut-off)
    const g = gap * (1 - smoothstep(A.la.c[1] - 0.2, A.la.c[1] + 1.4, py));
    const de = visceralDist(px, py, pz, A); // signed distance to the visceral envelope (epicardium + CS in its fat)
    if (de > 0.02 && de < g) {
      // fluid only displaces the pericardial space itself — never a vessel, lung or liver
      const b0 = bodyClassify(px, py, pz, A, g);
      if (b0.code === BODY.OUT || b0.code === BODY.PERI || b0.code === BODY.FAT) return cls(TISSUE.PERICARDIUM, 0.02);
    }
    periOff = g;
  }

  // --- valves (thin leaflets) tested first so they read over the blood pool ---
  const vtis = valveTissue(px, py, pz, G, path);
  if (vtis) return vtis;

  // --- chordae tendineae: fine strands tethering the AV leaflets to the papillary
  // muscles. Painted over the ventricular blood pool; taut (brighter) in systole,
  // slack in diastole — mitral off the LV papillaries, tricuspid off the RV. ---
  if (G.chordae) {
    let best = 1e9, bestTv = false;
    for (let i = 0; i < G.chordae.length; i++) {
      const c = G.chordae[i];
      const d = segDist(px, py, pz, c.a, c.b) - c.r;
      if (d < best) { best = d; bestTv = !!c.tv; }
    }
    if (best < 0) {
      const opening = G.valves ? (bestTv ? G.valves.tricuspid : G.valves.mitral) : 0;
      const taut = 1 - opening;                            // 1 shut … 0 open
      return cls(TISSUE.VALVE, 0.5 + 0.12 * taut);         // fine, faint strands
    }
  }

  // --- body: blood pools / myocardium / great-vessel walls from the SDF ---
  const body = bodyClassify(px, py, pz, A, periOff);
  // membranous septum: a thin fibrous segment of the septum just below the right/
  // non-coronary commissure (the perimembranous VSD site; a hole there when the
  // defect is perimembranous)
  if (body.code === BODY.MYO) {
    const ms = membSite(G);
    if (ms.inside(px, py, pz)) {
      // (perimembranous VSD: a ~0.9 cm hole through the membrane along the shunt axis)
      if (path.vsd === 'perimembranous' && offAxisDist(px, py, pz, ms.c, ms.rad) < 0.45) return cls(TISSUE.LV, 0.02);
      if (lumenDist(px, py, pz, A, 'LV') > MEMB_T) return cls(TISSUE.RV, 0.03);   // the RV side reaches the membrane
      return cls(TISSUE.VWALL, 0.55);
    }
  }
  // VSD: a real ~1 cm defect through the muscular septum along the shunt axis
  if (path.vsd && path.vsd !== 'perimembranous' && body.code === BODY.MYO && offAxisDist(px, py, pz, VSD_CORE, AX_VSD) < VSD_RADIUS &&
      Math.abs((px - VSD_CORE[0]) * AX_VSD[0] + (py - VSD_CORE[1]) * AX_VSD[1] + (pz - VSD_CORE[2]) * AX_VSD[2]) < 1.6)
    return cls(TISSUE.LV, 0.02);
  return cls(BODY_TO_TISSUE[body.code], body.echo);
}

// Unsigned distance from a point to the segment a—b (chordae are capsules).
function segDist(px, py, pz, a, b) {
  const bx = b[0] - a[0], by = b[1] - a[1], bz = b[2] - a[2];
  const px0 = px - a[0], py0 = py - a[1], pz0 = pz - a[2];
  const bb = bx * bx + by * by + bz * bz;
  let t = bb > 1e-9 ? (px0 * bx + py0 * by + pz0 * bz) / bb : 0;
  t = t < 0 ? 0 : t > 1 ? 1 : t;
  const dx = px0 - bx * t, dy = py0 - by * t, dz = pz0 - bz * t;
  return Math.hypot(dx, dy, dz);
}

// Inside the aortic root / ascending aorta lumen (anatomical SDF, behind a cheap
// bound around the root so it costs nothing for most samples).
const _AO_C = [LM.A[0] + LM.U_AO[0] * 1.6, LM.A[1] + LM.U_AO[1] * 1.6, LM.A[2] + LM.U_AO[2] * 1.6];
function aortaLumen(px, py, pz, G) {
  const dx = px - _AO_C[0], dy = py - _AO_C[1], dz = pz - _AO_C[2];
  if (dx * dx + dy * dy + dz * dz > 16) return false;
  return lumenDist(px, py, pz, G.A, 'AOROOT') < 0;
}

// Separate singleton from _clsResult so a valve hit and a subsequent non-valve
// classify() return can never share (alias) the same object.
// `calc` flags calcified / rheumatic leaflets (AS cusps, MS leaflets): the
// renderer keys acoustic shadowing and the bright calcific look on it, never on
// brightness, so normal valves cannot shadow.
const _valveResult = { tissue: TISSUE.VALVE, echo: 0, calc: 0 };

// ---------------------------------------------------------------------------
// Atrioventricular valves: hinged leaflets on a D-shaped, saddle-shaped annulus.
// In the valve frame (u = antero-posterior, +anterior leaflet side; c =
// intercommissural; a = axial, +downstream) each leaflet is a surface swept
// along c: a BODY hinging at the annulus plus a COAPTATION segment. Shut, the
// bodies reach inward to meet along the coaptation line (the short-axis
// "smile"), tented ~0.5 cm into the ventricle, and the coaptation segments lie
// against each other (~0.4 cm of apposition); open, the bodies swing toward the
// walls — the anterior mitral leaflet to near the septum — leaving the
// fish-mouth orifice between the tips. Leaflet lengths follow from the annulus:
// anterior mitral ~2.2 cm, posterior ~1.4 cm (scalloped P1-P3), shorter toward
// the commissures. Dimensions: mitral annulus ~2.8 (AP) x 3.2 cm (IC), saddle
// height ~0.6 cm; tricuspid ~3.0 x 3.4 cm with a short septal leaflet (+u side).
const AVL = {
  mitral: { Rap: 1.38, Ric: 1.62, saddle: 0.3, coapt: 0.36, ovl: 0.45, tent: 0.45,
            openA: 1.72, openP: 1.25, thickA: 0.075, thickP: 0.085, dShape: true, scallops: 3 },
  tricuspid: { Rap: 1.5, Ric: 1.72, saddle: 0.15, coapt: 0.62, ovl: 0.4, tent: 0.35,
               openA: 1.35, openP: 1.4, thickA: 0.07, thickP: 0.07, dShape: false, scallops: 0 },
};
// Leaflet polyline in the (u, a) plane at intercommissural fraction q, for side
// sgn (+1 anterior leaflet hinged at uA, -1 posterior hinged at uP). Writes the
// hinge, body end and tip into _lp. `mod` carries pathology: {ms, prolapse}.
const _lp = { uh: 0, ah: 0, ub: 0, ab: 0, ut: 0, at: 0, len: 1 };
function avLeafletPoly(P, sgn, q, o, sc, mod) {
  const q2 = q * q > 1 ? 1 : q * q;
  const half = Math.sqrt(1 - q2);
  const Rap = P.Rap * sc;
  const uA = P.dShape ? Rap * 0.97 * Math.sqrt(Math.max(0, 1 - q2 * q2 * q2)) : Rap * half;
  const uP = -Rap * half;
  const w = uA - uP;
  const ah = -P.saddle * (1 - 2 * q2);                      // saddle: atrial at A/P, low at the commissures
  const uc = uP + P.coapt * w;
  // functional (tethered) regurgitation: the papillary muscles pull the
  // coaptation point further into the ventricle (DCM mitral tenting)
  const tent = P.tent + (mod && mod.tent ? mod.tent : 0);
  const ac = ah + tent;
  const uh = sgn > 0 ? uA : uP;
  const reach = sgn > 0 ? uA - uc : uc - uP;               // horizontal reach to the coaptation line
  const Lb = Math.sqrt(reach * reach + tent * tent);
  const ovl = P.ovl * half + 0.08;
  let oo = o;
  if (mod && mod.ms && q2 > 0.2) oo *= 0.12;                // rheumatic commissural fusion
  const ac0 = Math.atan2(tent, reach);                      // shut body angle (tenting)
  const ao = sgn > 0 ? P.openA : P.openP;                   // fully-open body angle
  let alpha = ac0 + (ao - ac0) * oo;
  let beta = Math.PI / 2 * (1 - oo) + alpha * oo;           // coaptation segment: along the axis when shut
  if (mod && mod.ms) {
    // doming: the body bellies into the LV while the tethered tip is held back —
    // the diastolic "hockey-stick" anterior leaflet of rheumatic mitral stenosis
    alpha = ac0 + (1.25 - ac0) * oo;
    beta = ac0 + 0.25 * oo;
  }
  if (mod && mod.prolapse > 0 && sgn < 0 && q2 < 0.16 && o < 0.5) {
    // P2 prolapse: the middle posterior scallop billows back above the annular
    // plane in systole and fails to meet the anterior leaflet (primary MR)
    const k = mod.prolapse * (1 - o * 2) * (1 - q2 / 0.16);
    alpha = alpha * (1 - k) + (-0.55) * k;
    beta = beta * (1 - k) + (-0.2) * k;
  }
  _lp.uh = uh; _lp.ah = ah;
  _lp.ub = uh - sgn * Lb * Math.cos(alpha); _lp.ab = ah + Lb * Math.sin(alpha);
  _lp.ut = _lp.ub - sgn * ovl * Math.cos(beta); _lp.at = _lp.ab + ovl * Math.sin(beta);
  if (mod && mod.ms) {
    // planimetry-true orifice: the fused tips sit on an ellipse of AP width msGap
    // (commissures fused beyond |q| ≈ 0.45), opening with the transmitral flow
    const on = o / 0.35 > 1 ? 1 : o / 0.35;
    const e = 1 - q2 / 0.2;
    const hg = e > 0 ? 0.5 * mod.msGap * sc * Math.sqrt(e) * on : 0;
    // rheumatic doming ("hockey stick"): the fused tip sits at the orifice edge,
    // well into the LV, and the body bellies convexly toward the LV between hinge
    // and tip instead of running straight (leaflets never cross the orifice)
    const Lt = Lb + ovl;
    _lp.ut = uc + sgn * hg;
    const reachT = Math.abs(_lp.ut - uh);
    _lp.at = ah + Math.sqrt(Math.max(0.3, Lt * Lt * 0.8 - reachT * reachT)) * (0.35 + 0.65 * on);
    const bulge = (sgn > 0 ? 0.32 : 0.18) * on;               // anterior leaflet domes most
    _lp.ub = uh + (_lp.ut - uh) * 0.5 - sgn * bulge * 0.4;
    _lp.ab = ah + (_lp.at - ah) * 0.62 + bulge;
  }
  if (mod && mod.gap && o < 0.5) {
    // malcoaptation: with a dilated annulus the leaflets no longer meet — the
    // tips stop short of each other, leaving the regurgitant orifice visible
    const g = mod.gap * (1 - 2 * o) * half;
    _lp.ut += sgn * g * 0.5; _lp.ub += sgn * g * 0.25;
  }
  _lp.len = Lb + ovl;
  return _lp;
}
// distance from (u, a) to segment (u0, a0)-(u1, a1); also sets _segT (0..1 along it)
let _segT = 0;
function seg2d(u, a, u0, a0, u1, a1) {
  const du = u1 - u0, da = a1 - a0;
  const L2 = du * du + da * da || 1e-9;
  let t = ((u - u0) * du + (a - a0) * da) / L2;
  t = t < 0 ? 0 : t > 1 ? 1 : t;
  _segT = t;
  const x = u - (u0 + du * t), y = a - (a0 + da * t);
  return Math.sqrt(x * x + y * y);
}
// Is (u, a, c) inside a leaflet? Returns the leaflet-length fraction (0 hinge ..
// 1 tip) of the hit, or -1.
function avLeafletHit(P, u, a, c, o, sc, mod, thickMul) {
  const Ric = P.Ric * sc;
  const q = c / Ric;
  if (q <= -1.02 || q >= 1.02) return -1;
  for (let sgn = 1; sgn >= -1; sgn -= 2) {
    const L = avLeafletPoly(P, sgn, q, o, sc, mod);
    // posterior mitral scallops (P1/P2/P3): the free edge is notched between them
    let edgeCut = 1;
    if (sgn < 0 && P.scallops) edgeCut = 0.9 + 0.1 * Math.abs(Math.cos(P.scallops * Math.PI * 0.5 * (q + 1)));
    const d1 = seg2d(u, a, L.uh, L.ah, L.ub, L.ab);
    const t1 = _segT;
    const d2 = seg2d(u, a, L.ub, L.ab, L.ut, L.at);
    const t2 = _segT;
    const bodyFrac = 1 - (P.ovl + 0.08) / L.len;
    const f = d1 <= d2 ? t1 * bodyFrac : bodyFrac + t2 * (1 - bodyFrac);
    if (f > edgeCut) continue;
    const th0 = (sgn > 0 ? P.thickA : P.thickP) * thickMul;
    // thin membrane tapering to the free edge; rheumatic tips are bulbous
    const th = mod && mod.ms ? th0 * (1.1 + 0.9 * f) : th0 * (1.05 - 0.4 * f);
    if ((d1 < d2 ? d1 : d2) <= th) return f;
  }
  return -1;
}
// Mitral pathology modifiers: rheumatic stenosis (commissural fusion, doming),
// and a posterior (P2) prolapse whose extent grows with the MR grade.
const PROLAPSE = { mild: 0, moderate: 0.55, severe: 1 };
// tricuspid malcoaptation gap (cm) in TR, by grade
const TV_MOD = { mild: { gap: 0.1 }, moderate: { gap: 0.25 }, severe: { gap: 0.45 } };
const _modCache = new WeakMap();
function avMod(path) {
  let m = _modCache.get(path);
  if (!m) {
    // MS orifice: fused commissures leave an ellipse ~1.45 cm wide (|q| < 0.45);
    // its AP opening is set so the planimetered area matches the grade's MVA.
    const msA = MS_AREA[gradeOf(path)] || MS_AREA.severe;
    m = { ms: !!path.mitralStenosis, msGap: 4 * msA / (Math.PI * 1.45), prolapse: path.mr ? PROLAPSE[gradeOf(path)] : 0,
      tent: path.dilated && !path.mr ? 0.55 : 0 };
    _modCache.set(path, m);
  }
  return m;
}
// World-space free-edge (tip) point of a leaflet at intercommissural fraction q.
function avTipWorld(name, G, sgn, q) {
  const P = AVL[name], lf = LEAFLETS[name], vc = G.A.valves[name].c;
  const o = G.valves[name], sc = name === 'mitral' ? G.A.lv.mvR / LM.MV_R : 1;
  const mod = name === 'mitral' ? avMod(G.path || {}) : null;
  const L = avLeafletPoly(P, sgn, q, o, sc, mod);
  const c = q * P.Ric * sc;
  const t = [
    vc[0] + lf.ap[0] * L.ut + lf.flow[0] * L.at + lf.com[0] * c,
    vc[1] + lf.ap[1] * L.ut + lf.flow[1] * L.at + lf.com[1] * c,
    vc[2] + lf.ap[2] * L.ut + lf.flow[2] * L.at + lf.com[2] * c,
  ];
  // the mitral apparatus rides the angle-dependent annular excursion
  if (name === 'mitral') t[1] += mitralLift(t[0], t[1], t[2], G.A);
  return t;
}

function valveTissue(px, py, pz, G, path) {
  // apex-anchored contraction descends the AV valve planes with the annulus: a
  // valve at end-diastolic height y maps to apexY + (y−apexY)·lsy (same material
  // map as the LV wall). Mitral follows the LV; aortic sits at the LV base too.
  const live = G.A.valves;
  for (const name in VALVES) {
    const v = VALVES[name];
    const lf = LEAFLETS[name];
    const vc = live[name].c;
    // vector from the live annulus centre to the sample, and its split into the
    // leaflet axial (downstream) component `a` and the in-plane radial vector.
    const dx = px - vc[0], dz = pz - vc[2];
    // (the mitral hinges follow the angle-dependent annular excursion)
    const dy = py - vc[1] - (name === 'mitral' ? mitralLift(px, py, pz, G.A) : 0);
    if (dx * dx + dy * dy + dz * dz > 9) continue;       // cheap reject (>3 cm away)
    const f = lf.flow;
    const a = dx * f[0] + dy * f[1] + dz * f[2];          // axial: 0 at hinge, +downstream
    const rx = dx - a * f[0], ry = dy - a * f[1], rz = dz - a * f[2];
    const rad = Math.hypot(rx, ry, rz);                   // in-plane radius from axis
    if (rad > v.r + (AVL[name] ? 0.9 : 0.2)) continue;   // outside the annulus footprint
    const open = G.valves[name];                          // 0 shut … 1 fully open
    if (AVL[name]) {
      // hinged bileaflet AV valve (see avLeafletPoly)
      if (a < -0.9 || a > 3.2) continue;
      const u = rx * lf.ap[0] + ry * lf.ap[1] + rz * lf.ap[2];
      const c = rx * lf.com[0] + ry * lf.com[1] + rz * lf.com[2];
      const ms = name === 'mitral' && !!path.mitralStenosis;
      const mod = name === 'mitral' ? avMod(path) : (path.tr ? TV_MOD[gradeOf(path)] : null);
      // tricuspid annular dilatation in TR / pulmonary hypertension (~20 %)
      const sc = name === 'mitral' ? G.A.lv.mvR / LM.MV_R : (path.tr || path.rvpo ? 1.2 : 1);
      const f = avLeafletHit(AVL[name], u, a, c, open, sc, mod, ms ? 1.9 : 1);
      if (f < 0) continue;
      _valveResult.tissue = TISSUE.VALVE;
      _valveResult.echo = ms ? 0.98 : 0.8 + 0.06 * (1 - f);
      _valveResult.calc = ms ? 1 : 0;
      return _valveResult;
    }
    // Semilunar cusps: three curtains hinging at the annulus.
    let cd = lf.cd, thick = lf.thick;
    if (lf.ap) {
      const side = (rx * lf.ap[0] + ry * lf.ap[1] + rz * lf.ap[2]) / (rad + 1e-6);
      cd *= 1 + (lf.asym - 1) * Math.max(0, side);        // anterior (+ap) drapes deeper
      if (side < 0) {                                     // posterior leaflet
        thick = lf.thickP || thick;                       // shorter but a touch thicker
        if (lf.com && lf.scallops) {                      // P1/P2/P3 scalloped free edge
          const comCoord = (rx * lf.com[0] + ry * lf.com[1] + rz * lf.com[2]) / (v.r + 1e-6);
          cd *= 1 + 0.12 * Math.cos(lf.scallops * Math.PI * comCoord) * -side;
        }
      }
    }
    if (a < -0.05 || a > cd) continue;                    // only the leaflet span, hinge→edge
    const t = a <= 0 ? 0 : a / cd;                        // 0 at hinge, 1 at free edge
    // Free-edge orifice radius. AV valves iris symmetrically; the semilunar valves
    // are TRILOBED — three cusps whose orifice bulges toward the commissures and
    // pinches at the cusp centres (the triangular systolic aortic orifice).
    let edgeR = open * v.r, kh = 0, cptr = 0, isCusp = lf.cusps === 3;
    if (isCusp) {
      const ang = Math.atan2(rx * lf.w[0] + ry * lf.w[1] + rz * lf.w[2],
                             rx * lf.u[0] + ry * lf.u[1] + rz * lf.u[2]);
      const k = 3 * (ang - lf.commOff);
      kh = Math.atan2(Math.sin(k), Math.cos(k));          // angle to nearest commissure ×3
      cptr = 0.5 * (Math.cos(k) + 1);                     // 1 at commissures, 0 at cusp centres
      edgeR = open * v.r * (0.58 + 0.42 * cptr);          // triangular orifice
    }
    // Curved cusp: an S-profile (smoothstep) so it leaves the annulus nearly
    // parallel to the outflow axis, bellies, then curves in to the free edge —
    // a real doming/tenting curtain, not a straight cone.
    // A closing semilunar cusp is a hammock: it hugs the sinus wall and turns
    // in to the centre only near its free edge, so a short-axis cut above the
    // coaptation line shows the sinus wall, not a concentric ring of cusp.
    let s = t * t * (3 - 2 * t);
    if (isCusp) s += (t * t * t - s) * (1 - open);
    const shellR = v.r + (edgeR - v.r) * s;               // curved radius at this depth
    // Thin membrane tapering to a fine free edge (thickest at the annular base).
    let th = thick * (0.45 + 0.55 * (1 - t) * (1 - t));
    // Nodulus of Arantius: a fibrous thickening at the CENTRE of each semilunar
    // cusp's free edge, where the three cusps meet — bulk up the free-edge tip at
    // the cusp centre so the coaptation shows the three little nodules.
    let nod = 0;
    if (isCusp && t > 0.7) nod = 0.11 * ((t - 0.7) / 0.3) * (1 - cptr);
    th += nod;
    // calcific AS: thick (3-5 mm), nodular cusps — calcium masses concentrated at the
    // cusp bodies and bases rather than a thin bright membrane
    if (name === 'aortic' && path.aorticStenosis) {
      const lump = 0.5 + 0.5 * Math.sin(rx * 7.1 + ry * 5.3 + rz * 6.7);
      th = th * 2.6 + 0.06 + 0.06 * lump * (1 - t * 0.5);
    }
    let hit = Math.abs(rad - shellR) <= th;
    // Semilunar commissural coaptation seams: three radial lines meeting centrally
    // as the cusps shut — the short-axis "Mercedes" Y, and the closure line in LAX.
    if (!hit && isCusp) {
      const closed = 1 - open;
      const angDist = Math.abs(kh) / 3;                   // angular distance to a commissure
      if (closed > 0.25 && a > 0.6 * cd && rad < v.r * 0.9 && angDist * rad < 0.07) hit = true;
    }
    if (!hit) continue;
    let echo = 0.8 + 0.05 * (1 - t) + (nod > 0.03 ? 0.04 : 0); // nodule reads a touch brighter
    const calc = (name === 'aortic' && path.aorticStenosis) || (name === 'mitral' && path.mitralStenosis);
    if (calc) echo = 0.98; // calcified / restricted → bright, thickened
    _valveResult.tissue = TISSUE.VALVE;
    _valveResult.echo = echo;
    _valveResult.calc = calc ? 1 : 0;
    return _valveResult;
  }
  return null;
}

// ---------------------------------------------------------------------------
// Flow-field helpers. A jet/tract is a stream about an axis, not a chamber
// fill: speed is a fast on-axis core that falls off with a Gaussian in the
// perpendicular (off-axis) direction and decays along the direction travelled.
// ---------------------------------------------------------------------------

// Signed pulmonary-venous flow velocity (m/s): positive = forward, toward the
// left atrium; negative = retrograde, back up the vein. Normal adults show S >= D;
// the small late-diastolic 'a' reversal is the atrium contracting against open
// vein mouths. Severe MR inverts the systolic component (systolic flow reversal).
function pvFlowSpeed(phase, path) {
  const s = pulse(0.05, 0.35, phase);        // systolic (S) wave
  const d = pulse(0.50, 0.72, phase);        // early-diastolic (D) wave
  const a = pulse(0.86, 1.0, phase);         // atrial reversal
  const severeMR = !!path.mr && gradeOf(path) === 'severe';
  const sComp = severeMR ? -0.55 * s : 0.48 * s;
  return sComp + 0.42 * d - 0.20 * a;
}

// Off-axis Gaussian falloff: 1 on the centreline, ->0 as perpendicular distance
// d grows past the stream half-width w.
function gauss(d, w) {
  const t = d / w;
  return Math.exp(-t * t);
}

// Perpendicular distance from point (px,py,pz) to the line through anchor `a`
// with UNIT direction `u`. The signed distance travelled along the axis is
// stashed in _axis.along (module scope) to avoid allocating a return tuple.
const _axis = { along: 0 };
function offAxisDist(px, py, pz, a, u) {
  const dx = px - a[0], dy = py - a[1], dz = pz - a[2];
  const along = dx * u[0] + dy * u[1] + dz * u[2];
  _axis.along = along;
  const rx = dx - along * u[0], ry = dy - along * u[1], rz = dz - along * u[2];
  return Math.hypot(rx, ry, rz);
}

// Precomputed unit flow axes (constant in heart space). Anchors reference the
// valve centres in VALVES; jet cores are given explicitly below.
// All derived from the anatomy landmarks so the flow follows the rendered heart.
const _sub3 = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const _dot3 = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const LA_MID = [LM.LA_FLOOR0[0] + LM.A_LA[0] * 2.3, LM.LA_FLOOR0[1] + LM.A_LA[1] * 2.3, LM.LA_FLOOR0[2] + LM.A_LA[2] * 2.3];
const RA_MID = [LM.RA_FLOOR0[0] + LM.A_RA[0] * 2.1, LM.RA_FLOOR0[1] + LM.A_RA[1] * 2.1, LM.RA_FLOOR0[2] + LM.A_RA[2] * 2.1];
const AX_MITRAL_IN = _mFlow;                    // mitral annulus -> LV apex
const AX_TRICUSPID_IN = _tvFlow;                // tricuspid annulus -> RV apex
const AX_LVOT = unit(_sub3(LM.A, LM.LVOT0));    // LVOT -> aorta
const AX_RVOT = unit(_sub3(LM.PV, LM.RVOT[1])); // infundibulum -> pulmonary artery
const AX_MR = unit(_sub3(LA_MID, LM.M));        // mitral -> into the LA
// posterior-leaflet (P2) prolapse drives an ECCENTRIC jet away from the flail
// leaflet: anteriorly, hugging the anterior LA wall behind the aortic root
// (functional MR in DCM stays central)
const AX_MR_ECC = (() => {
  const toAo = unit(_sub3(LM.A, LM.M));
  const t = _sub3(toAo, [AX_MR[0] * _dot3(toAo, AX_MR), AX_MR[1] * _dot3(toAo, AX_MR), AX_MR[2] * _dot3(toAo, AX_MR)]);
  const tn = unit(t);
  return unit([AX_MR[0] + 0.75 * tn[0], AX_MR[1] + 0.75 * tn[1], AX_MR[2] + 0.75 * tn[2]]);
})();
const AX_TR = unit(_sub3(RA_MID, LM.T));        // tricuspid -> into the RA
const AX_VSD = unit([-1, 0.1, 0.15]);           // LV -> RV across the septum
const AX_ASD = LM.IAS_N;                        // LA -> RA across the septum
// mid muscular septum, where the mid-papillary PSAX and the A4C planes cross, so the
// defect is seen in both standard views
const VSD_CORE = [-2.4, -3.6, LM.A4C_Z(-2.4, -3.6)];
const ASD_CORE = LM.IAS.fossaC;      // fossa ovalis (secundum ASD)

// Membranous septum: the thin fibrous part of the septum just below the commissure
// between the right-coronary and non-coronary cusps (the one that faces the right
// heart), where the LVOT wall meets the RA/RV. The muscular septum is ~1.3 cm there;
// an ellipsoid pit opens the RV/RA side of it down to a MEMB_T membrane on the LV
// side. Follows the aortic root, so it moves with the valve plane. Memoised per
// anatomy object.
export const MEMB_T = 0.12;                                        // membrane thickness (cm)
const MEMB_R = [1.3, 0.8, 0.6];                             // pit semi-axes: radial (toward the right heart), tangential, axial (cm)
const _membCache = new WeakMap();
export function membSite(G) {
  let s = _membCache.get(G.A);
  if (s) return s;
  const A = G.A.valves.aortic.c, U = LM.U_AO;
  const toT = _sub3(G.A.valves.tricuspid.c, A);
  let rad = null, best = -1e9;
  for (const deg of [30, 150, 270]) {                       // the three commissures
    const a = deg * Math.PI / 180;
    const r = [0, 1, 2].map((i) => LM.E_SCR[i] * Math.cos(a) + LM.E_ANT[i] * Math.sin(a));
    const d = _dot3(r, toT);
    if (d > best) { best = d; rad = r; }
  }
  const tan = unit([U[1] * rad[2] - U[2] * rad[1], U[2] * rad[0] - U[0] * rad[2], U[0] * rad[1] - U[1] * rad[0]]);
  const c = [0, 1, 2].map((i) => A[i] + rad[i] * 2.3 - U[i] * 0.5);
  s = {
    c, rad, tan,
    inside(px, py, pz) {
      const q0 = px - c[0], q1 = py - c[1], q2 = pz - c[2];
      const ur = (q0 * rad[0] + q1 * rad[1] + q2 * rad[2]) / MEMB_R[0];
      const ut = (q0 * tan[0] + q1 * tan[1] + q2 * tan[2]) / MEMB_R[1];
      const ua = (q0 * U[0] + q1 * U[1] + q2 * U[2]) / MEMB_R[2];
      return ur * ur + ut * ut + ua * ua < 1;
    },
  };
  _membCache.set(G.A, s);
  return s;
}
// Regurgitant-jet spread regions (ellipsoid centres): MR reaches deep into the
// LA body, TR into the RA body. Used to bound each jet's colour to its receiving
// atrium so no mosaic leaks across the septum into the other side of the heart.
const MR_REGION = LA_MID;  // LA body (mitral regurgitation target)
const TR_REGION = RA_MID;  // RA body (tricuspid regurgitation target)
// which side of the interatrial septum a point lies on (<0 = LA side)
const iasSide = (px, py, pz) => (px - LM.IAS_P[0]) * LM.IAS_N[0] + (py - LM.IAS_P[1]) * LM.IAS_N[1] + (pz - LM.IAS_P[2]) * LM.IAS_N[2];

// Diastolic mitral-inflow vortex. A pure swirl about a fixed axis is
// divergence-free (∇·v = 0), so superimposing it on the inflow jet keeps the
// field mass-consistent while reproducing the physiologic filling vortex that
// curls behind the anterior mitral leaflet on colour Doppler. Axis ~ anterior-
// posterior so the recirculation reads in the A4C / PLAX planes.
const VORTEX_AXIS = unit([0, 0, 1]);
const VORTEX_CORE = [
  VALVES.mitral.c[0] + AX_MITRAL_IN[0] * 1.5,
  VALVES.mitral.c[1] + AX_MITRAL_IN[1] * 1.5,
  VALVES.mitral.c[2] + AX_MITRAL_IN[2] * 1.5,
];
const VORTEX_R = 1.3;          // core radius (cm)
const VORTEX_STRENGTH = 0.45;  // swirl gain (m/s per cm, before the Gaussian roll-off)

// Module-scoped singleton velocity result + best-tracking flag. velocityAt()
// returns THIS object (or null) every call; callers read it immediately.
const _velResult = { vx: 0, vy: 0, vz: 0, speed: 0, flow: FLOW.NONE, turbulent: false };
let _velHasBest = false;

// Lifted out of velocityAt so no closure is allocated per sample; state (the
// current best) lives in the module-scoped singleton above. Keeps the fastest
// contribution. dir components are expected pre-normalised (unit axes).
function considerFlow(inside, dirx, diry, dirz, speed, flow, turbulent) {
  if (!inside || speed <= 0.01) return;
  if (_velHasBest && speed <= _velResult.speed) return;
  _velResult.vx = dirx * speed;
  _velResult.vy = diry * speed;
  _velResult.vz = dirz * speed;
  _velResult.speed = speed;
  _velResult.flow = flow;
  _velResult.turbulent = turbulent;
  _velHasBest = true;
}

// Shared regurgitant-jet contribution (MR / TR are geometrically identical: a
// systolic vena-contracta core firing across a shut AV valve into the receiving
// atrium, plus an upstream ~1/r PISA flow-convergence hemisphere on the
// ventricular side). `downstreamGate` is the already-evaluated "sample is inside
// the jet's receiving-atrium region" test; `upstreamInside` is the ventricle
// containment for PISA. `p` carries the per-lesion tuning (peak velocity, width,
// decay, PISA gains, turbulence thresholds).
function regurgJet(px, py, pz, phase, valveC, axis, downstreamGate, upstreamInside, p, env, peak) {
  // env / peak: when supplied (MR) the jet timing comes from the modelled
  // regurgitant flow waveform and the peak from the modelled LV-LA gradient
  // (Bernoulli); otherwise fall back to the analytic pulse + tuned peak (TR).
  // (TR spans ventricular systole incl. the isovolumic periods: TV closure at the
  // QRS to TV opening at ~0.50; MR passes its modelled envelope instead)
  const jet = env != null ? env : pulse(0.02, 0.48, phase);
  const pk = peak != null ? peak : p.peak;
  // downstream core: fast on-axis, Gaussian off-axis, decaying along travel
  const perp = offAxisDist(px, py, pz, valveC, axis);
  const along = _axis.along;
  const w = p.w0 + p.wSlope * clamp(along, 0, p.wClamp);
  const shape = gauss(perp, w) * Math.exp(-clamp(along, 0, p.decayClamp) / p.decay);
  const speed = jet * pk * shape;
  considerFlow(downstreamGate, axis[0], axis[1], axis[2], speed, p.flow, jet > p.coreTurb);
  // proximal isovelocity surface area (PISA): on the ventricular (upstream)
  // side the flow converges radially toward the orifice, speeding up (~1/r).
  const dOrif = Math.hypot(px - valveC[0], py - valveC[1], pz - valveC[2]);
  if (along < 0 && upstreamInside && dOrif < p.pisaR) {
    const inv = 1 / (dOrif + 1e-4);
    const conv = jet * clamp(p.pisaGain / (dOrif + 0.2), 0, p.pisaMax);
    considerFlow(true, (valveC[0] - px) * inv, (valveC[1] - py) * inv, (valveC[2] - pz) * inv,
      conv, p.flow, jet > p.pisaTurb && dOrif < 0.6);
  }
}

// Per-lesion tuning for regurgJet. Core widths are kept broad enough that the
// finite echo sampling grid reliably hits the high-velocity core, so the tracked
// spectral/Peak-V reads the true SEVERE value (~5.5 MR) rather than under-reading
// a sub-pixel vena contracta; the jet still stays strictly within the LA.
const MR_JET_PARAMS = { w0: 1.0, wSlope: 0.12, wClamp: 5, decayClamp: 9, decay: 9.0,
  peak: 5.7, flow: FLOW.MR_JET, coreTurb: 0.2, pisaR: 1.6, pisaGain: 0.6, pisaMax: 4, pisaTurb: 0.25 };
const TR_JET_PARAMS = { w0: 1.05, wSlope: 0.12, wClamp: 4, decayClamp: 6, decay: 5.0,
  peak: 3.2, flow: FLOW.TR_JET, coreTurb: 0.2, pisaR: 1.5, pisaGain: 0.5, pisaMax: 3, pisaTurb: 0.25 };

// Radius (cm) of the modelled VSD channel. A ~1 cm defect with a ~4.5 m/s jet is a
// RESTRICTIVE defect: the model has no L-R shunt volume (no Qp:Qs, LA/LV overload).
export const VSD_RADIUS = 0.5;

// Shunt lesions' badge: size claims must follow the geometry and the circulation.
// A 'large' VSD is non-restrictive (radius >= 0.9 cm) with Qp:Qs >= 1.5, which this
// model does not represent, so its defects are labelled restrictive / secundum.
export function shuntLabel(path = {}) {
  if (path.vsd) {
    return { severity: 'restrictive shunt', label: path.vsd === 'perimembranous' ? 'restrictive perimembranous VSD' : 'restrictive muscular VSD' };
  }
  if (path.asd) return { severity: 'shunt', label: 'secundum ASD' };
  return null;
}

// What a CW sweep through a stenotic jet seen in the plane reads: the modelled
// peak (the vena contracta), not whatever a 2-D sample happens to land on. `peak`
// is the cycle peak (m/s) and `env` its instantaneous 0..1 envelope; peak = 0 when
// `flow` is not a stenotic jet of this case.
const _cw = { peak: 0, env: 0 };
export function stenosisCw(flow, H, path = {}) {
  _cw.peak = 0; _cw.env = 0;
  if (!H) return _cw;
  if (flow === FLOW.AS_JET && path.aorticStenosis) { _cw.peak = H.vAoPeak; _cw.env = H.ejAo; }
  else if (flow === FLOW.MITRAL_IN && path.mitralStenosis) { _cw.peak = H.vMitPeak; _cw.env = H.ejMv; }
  return _cw;
}

// ---------------------------------------------------------------------------
// Blood-flow velocity field (cm/s) at a point & phase, for colour Doppler.
// Returns {vx,vy,vz, speed, flow, turbulent} (a reused singleton) or null.
// speed in m/s magnitude. Every compartment is a tapering stream: away from a
// jet/tract the returned velocity is null / near-zero.
// ---------------------------------------------------------------------------
export function velocityAt(px, py, pz, G, path = {}) {
  const phase = G.phase;
  const H = G.hemo;   // per-frame haemodynamic scalars (jet peaks + flow envelopes)
  _velHasBest = false;

  const inLV = ellip(px, py, pz, G.lv.c, G.lv.r) <= 1;
  const inRV = ellip(px, py, pz, G.rv.c, G.rv.r) <= 1;
  const inLA = ellip(px, py, pz, G.la.c, G.la.r) <= 1;
  const inRA = ellip(px, py, pz, G.ra.c, G.ra.r) <= 1;
  const inAo = aortaLumen(px, py, pz, G);

  // --- diastolic inflow (E + A): a cone from the annulus that narrows and
  //     decays toward the apex. STRICTLY gated to diastole (phase ~0.42-1.0);
  //     the AV valves are shut in systole so there must be no apex-directed
  //     ventricular colour then. The E/A pulse windows already lie inside
  //     diastole (E 0.50-0.66, A 0.84-1.0); the gate makes that guarantee hard.
  if (phase >= 0.42) {
    // mitral inflow -> LV apex. Peak velocity comes from the modelled inflow via
    // continuity (v = Q/A): H.vMitPeak already carries the E/A magnitude and MS
    // scaling (small MVA -> high, prolonged velocity); H.ejMv is the instantaneous
    // E-then-A flow envelope from the circulation model.
    const mPerp = offAxisDist(px, py, pz, VALVES.mitral.c, AX_MITRAL_IN);
    const mAlong = _axis.along;
    const mW = 1.25 * (1 - 0.5 * clamp(mAlong / 4.5, 0, 1)); // narrows toward apex
    const mShape = gauss(mPerp, mW) * Math.exp(-clamp(mAlong, 0, 6) / 4.5);
    const mSpeed = H.vMitPeak * H.ejMv * mShape;
    if (inLV && mAlong > -0.3 && mSpeed > 0.01) {
      // superimpose the divergence-free filling vortex: a pure swirl about
      // VORTEX_AXIS (tangential = axis × r_perp) whose magnitude rises then
      // decays with radius (Lamb-Oseen-like), so the inflow curls behind the
      // anterior leaflet instead of reading as a straight jet.
      const rx = px - VORTEX_CORE[0], ry = py - VORTEX_CORE[1], rz = pz - VORTEX_CORE[2];
      const axl = rx * VORTEX_AXIS[0] + ry * VORTEX_AXIS[1] + rz * VORTEX_AXIS[2];
      const qx = rx - axl * VORTEX_AXIS[0], qy = ry - axl * VORTEX_AXIS[1], qz = rz - axl * VORTEX_AXIS[2];
      const rr = Math.hypot(qx, qy, qz);
      const tx = VORTEX_AXIS[1] * qz - VORTEX_AXIS[2] * qy;
      const ty = VORTEX_AXIS[2] * qx - VORTEX_AXIS[0] * qz;
      const tz = VORTEX_AXIS[0] * qy - VORTEX_AXIS[1] * qx;
      const tlen = Math.hypot(tx, ty, tz) + 1e-6;
      const vth = VORTEX_STRENGTH * rr * Math.exp(-(rr / VORTEX_R) * (rr / VORTEX_R)) * H.ejMv;
      const vx = AX_MITRAL_IN[0] * mSpeed + (tx / tlen) * vth;
      const vy = AX_MITRAL_IN[1] * mSpeed + (ty / tlen) * vth;
      const vz = AX_MITRAL_IN[2] * mSpeed + (tz / tlen) * vth;
      const sp = Math.hypot(vx, vy, vz);
      considerFlow(true, vx / sp, vy / sp, vz / sp, sp, FLOW.MITRAL_IN,
        path.mitralStenosis && mSpeed > 1.5);
    }

    // tricuspid inflow -> RV apex (right heart not lumped-modelled; the filling
    // envelope on the right heart's clock, at a lower normal velocity)
    const tPerp = offAxisDist(px, py, pz, VALVES.tricuspid.c, AX_TRICUSPID_IN);
    const tAlong = _axis.along;
    const tW = 1.35 * (1 - 0.5 * clamp(tAlong / 4.5, 0, 1));
    const tShape = gauss(tPerp, tW) * Math.exp(-clamp(tAlong, 0, 6) / 4.5);
    const tSpeed = 0.55 * H.ejTv * tShape;
    considerFlow(inRV && tAlong > -0.3, AX_TRICUSPID_IN[0], AX_TRICUSPID_IN[1], AX_TRICUSPID_IN[2],
      tSpeed, FLOW.TRICUSPID_IN, false);
  }

  // --- systolic ejection: a THIN column hugging the outflow centreline only,
  //     NOT a whole-ventricle fill. Present only during ejection (0.04-0.34).
  if (H.ejAo > 0.001) {
    // LVOT -> aorta (centreline anchored at the aortic valve). Clip the stream
    // to the aortic (downstream, +along) side of the valve plane so it reads as
    // ONE directed jet aimed into the aorta rather than two symmetric lobes
    // flanking the annulus (the off-axis Gaussian alone is symmetric in `along`).
    const lPerp = offAxisDist(px, py, pz, VALVES.aortic.c, AX_LVOT);
    const lAlong = _axis.along;
    // LVOT/aortic-valve jet peak from continuity (v = SV / (AVA·ET)): with the
    // SAME stroke volume forced through a reduced AVA, AS auto-scales to ~4 m/s
    // and normal sits ~1 m/s — severity now follows physiology, not a constant.
    const lSpeed = H.vAoPeak * H.ejAo * gauss(lPerp, 0.7);
    considerFlow(((inLV && py > G.lv.c[1]) || inAo) && lAlong > -0.25,
      AX_LVOT[0], AX_LVOT[1], AX_LVOT[2],
      lSpeed, path.aorticStenosis ? FLOW.AS_JET : FLOW.LVOT, path.aorticStenosis && H.ejAo > 0.3);
  }
  // RV ejection: on the right heart's clock it starts ~10 ms before and ends
  // ~30 ms after the aortic, so it has its own gate
  if (H.ejPv > 0.001) {
    // RVOT -> pulmonary artery (centreline anchored at the pulmonic valve).
    // Tightened to the true outflow centreline: a narrow Gaussian, clipped to
    // the downstream (+along) side and to the high subpulmonary tract (py > 0.8)
    // so no stray systolic colour blob appears mid-RV near the tricuspid valve.
    const rPerp = offAxisDist(px, py, pz, VALVES.pulmonic.c, AX_RVOT);
    const rAlong = _axis.along;
    const rSpeed = H.ejPv * 0.9 * gauss(rPerp, 0.55); // systolic timing from the RV ejection envelope
    considerFlow(rAlong > -2.4 && rAlong < 3.5 && rPerp < 1.3, AX_RVOT[0], AX_RVOT[1], AX_RVOT[2],
      rSpeed, FLOW.RVOT, false);
  }

  // --- pathological jets: a fast vena-contracta core tapering laterally and
  //     distally (Gaussian off-axis + travel decay). ---
  if (path.mr || path.dilated) {
    // mitral regurgitation (primary, or functional in DCM): systolic jet LV -> LA. The region reaches from the
    // mitral annulus DEEP into the LA body, with a long travel-decay so the
    // brightest turbulent mosaic sits well within the atrium (mid-LA). The
    // downstream colour is confined STRICTLY to the LA side of the crux: the
    // sample must be inside the LA lumen (or, in the annulus->LA transition,
    // above the annulus AND on the LV/LA side of the septum) and must never fall
    // inside the RV/RA proxies — so no blue/mosaic leaks into the RV or septum.
    const region = ellip(px, py, pz, MR_REGION, [2.3, 3.0, 2.0]) <= 1;
    const aboveMV = _dot3(_sub3([px, py, pz], G.A.valves.mitral.c), VALVES.mitral.n) > 0.1;
    const downstream = region && !inRV && iasSide(px, py, pz) < 0 && (inLA || aboveMV);
    // envelope + peak from the modelled regurgitant flow / LV-LA gradient: the
    // regurgitant fraction sets how much colour fills the LA; the systolic
    // LV-LA pressure difference sets the ~5 m/s vena-contracta velocity.
    regurgJet(px, py, pz, phase, VALVES.mitral.c, path.mr ? AX_MR_ECC : AX_MR, downstream, inLV, MR_JET_PARAMS,
      H.ejMr, H.vMRPeak);
  }
  if (path.tr) {
    // tricuspid regurgitation: systolic jet RV -> RA, firing into the RA cavity.
    // The TR peak velocity encodes RV systolic pressure: with pressure overload
    // (pulmonary hypertension) the jet runs fast (~4.3 m/s → PASP ~80 mmHg via
    // 4v²+RAP); isolated TR at normal PA pressure runs ~3.0 m/s.
    const trPeak = path.rvpo ? 4.3 : 3.0;
    const region = ellip(px, py, pz, TR_REGION, [2.3, 2.8, 2.0]) <= 1;
    const aboveTV = _dot3(_sub3([px, py, pz], G.A.valves.tricuspid.c), VALVES.tricuspid.n) > 0.1;
    const downstream = region && iasSide(px, py, pz) > 0 && (inRA || aboveTV);
    regurgJet(px, py, pz, G.phaseRV != null ? G.phaseRV : phase, VALVES.tricuspid.c, AX_TR, downstream, inRV, TR_JET_PARAMS, null, trPeak);
  }
  if (path.vsd) {
    // ventricular septal defect: LV -> RV across septum, systole, turbulent core
    const jet = pulse(0.02, 0.44, phase);
    const peri = path.vsd === 'perimembranous';           // else the mid-muscular defect
    const ms = peri ? membSite(G) : null;
    const core = peri ? ms.c : VSD_CORE, ax = peri ? ms.rad : AX_VSD;
    const region = ellip(px, py, pz, core, peri ? [1.6, 1.6, 1.6] : [1.0, 1.6, 1.2]) <= 1;
    const perp = offAxisDist(px, py, pz, core, ax);
    const along = _axis.along;
    const shape = gauss(perp, 0.55) * Math.exp(-clamp(Math.abs(along), 0, 3) / 2.2);
    const speed = jet * 4.5 * shape;                      // restrictive VSD peak ~4.5 m/s
    considerFlow(region, ax[0], ax[1], ax[2], speed, FLOW.VSD_JET, jet > 0.15);
  }
  // --- pulmonary-vein flow -------------------------------------------------
  // Forward (vein -> LA) in systole (S wave, driven by atrial relaxation and
  // annular descent) and again in early diastole (D wave, once the mitral opens),
  // with a small RETROGRADE 'a' wave as the atrium contracts against the open
  // vein mouths. In SEVERE mitral regurgitation the systolic component blunts and
  // then REVERSES — regurgitant volume is driven back up the veins — which is an
  // ASE-specific criterion for severe MR and one of the most useful confirmatory
  // signs at the bedside. Gated to the posterior LA so the four vein tests stay
  // off the hot path everywhere else.
  const A = G.A;
  if (A && A.pv && pz < A.la.c[2] + 1.2) {
    const dcx = px - A.la.c[0], dcy = py - A.la.c[1], dcz = pz - A.la.c[2];
    if (dcx * dcx + dcy * dcy + dcz * dcz < 16) {          // within 4 cm of the LA centre
      const sp = pvFlowSpeed(phase, path);
      if (sp !== 0) {
        const veins = A.pv.veins;
        for (let i = 0; i < veins.length; i++) {
          const v = veins[i];
          const perp = offAxisDist(px, py, pz, v.b, v.ax); // ax precomputed per frame
          if (perp > v.r * 1.8) continue;                  // not in this vein's stream
          const s = sp * gauss(perp, v.r * 0.9);
          const dir = s >= 0 ? 1 : -1;                     // negative = retrograde
          considerFlow(true, v.ax[0] * dir, v.ax[1] * dir, v.ax[2] * dir,
            s < 0 ? -s : s, FLOW.PV_FLOW, false);
        }
      }
    }
  }
  if (path.asd) {
    // atrial septal defect: LA -> RA, low-velocity, mostly all-cycle
    const region = ellip(px, py, pz, ASD_CORE, [1.2, 1.3, 1.2]) <= 1;
    const perp = offAxisDist(px, py, pz, ASD_CORE, AX_ASD);
    const along = _axis.along;
    const shape = gauss(perp, 0.7) * Math.exp(-clamp(Math.abs(along), 0, 3) / 2.5);
    considerFlow(region, AX_ASD[0], AX_ASD[1], AX_ASD[2], 0.9 * shape, FLOW.ASD_JET, false);
  }

  return _velHasBest ? _velResult : null; // may be null (no meaningful flow)
}
