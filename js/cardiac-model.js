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
import { anatomyParams, rightHeartLoad, hepaticVeinSegs, ivcOfRap, atrialFill, bodyClassify, visceralDist, lumenDist, mitralLift, tricuspidLift, buildRwma, rwmaWeight, BODY, CFG, LM, BODY_AX } from './anatomy.js';
import { hemodynamics, hemoSummary, pvLoop, gradeOf, eaRatio, MS_AREA, laPhasic } from './hemodynamics.js';
import { DCM_SEVERITIES, dcmSeverityOf } from './scenarios.js';

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
  HV_FLOW: 11,    // hepatic-vein flow toward the RA (reverses in systole in severe TR)
  AORTA: 12,      // ascending / arch / descending aortic flow beyond the root
  LAA_FLOW: 13,   // left-atrial-appendage emptying (atrial systole) / filling (early systole)
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
  aortic:    { flow: _aoFlow, cd: 0.95, thick: 0.05, ap: null, asym: 1,
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

// Right-heart timing. The right heart is solved in the same closed loop as the
// left (hemodynamics.js: its own RV elastance, activated RV_TIMING.delay after
// the LV and relaxing on its own time; RA, tricuspid, pulmonary valve and PA), so
// its valve events (T1 after M1, P2 after A2, a shorter isovolumic relaxation),
// its volume and its flows come straight from the solved Qtv, Qpa and Vrv at the
// cycle phase. (Review 12: the earlier replay of the LEFT ventricle's trace on a
// shifted clock is retired.) What remains is phaseRV, the phase at which the
// atrial reservoir / conduit / booster curve, the TR-jet envelope and the hepatic
// veins (all drawn on left-heart timing) are read for the right heart: a warp
// that puts the left heart's AV-valve closure and opening onto the SOLVED
// tricuspid closure and opening of this circulation (rvPhase(phase, path)).
export function rvPhase(phase, path = null) {
  const p = ((phase % 1) + 1) % 1;
  if (!path) return p;
  const w = openTrace(path).warp;
  // anchors: [right-heart phase, left-heart phase] at AV-valve closure and opening
  const a0 = w[0], a1 = w[1];
  let x = p - a0[0]; x -= Math.floor(x);
  const span1 = ((a1[0] - a0[0]) % 1 + 1) % 1, span2 = 1 - span1;
  const l1 = ((a1[1] - a0[1]) % 1 + 1) % 1, l2 = 1 - l1;
  const q = x < span1 ? a0[1] + x * l1 / span1 : a1[1] + (x - span1) * l2 / span2;
  return ((q % 1) + 1) % 1;
}

// Valve opening 0..1, driven by the modelled transvalvular FLOW so the leaflets
// open when blood starts to cross and shut as it stops (the valve events — AVO,
// AVC, MVO, MVC — and the isovolumic periods all come from the circulation). An
// AV valve swings wide on the E wave, drifts toward closure in diastasis and
// reopens on the A wave; a semilunar valve opens fully in early ejection and
// floats toward closure as flow decelerates. The motion is continuous, never a
// one-frame step: a smooth map of the flow read slightly ahead (the leaflets
// answer the pressure gradient, which leads the inertial flow — further on the
// way down, so they are all but shut as forward flow ends), then a first-order
// lag. Closure takes ~20-40 ms, and the leaflets cross 10 % open within a few
// ms of the flow-defined events (AVO, AVC, MVO, MVC). Pre-computed once
// per circulation as a periodic trace indexed by phase; the right-heart valves
// read their own solved flows (tricuspid Qtv, pulmonary Qpa).
const OPEN_N = 960, OPEN_CYCLE = 0.833;                   // samples per cycle; reference cycle (s)
const OPEN_TAU_OPEN = 0.008, OPEN_TAU_CLOSE = 0.008, OPEN_LEAD_O = 0.007, OPEN_LEAD = 0.011;   // s
// (OPEN_Q_AO, the aortic cusps on the aortic flow: 0.2 of peak. Review 12 rounded
// the aortic ejection, whose deceleration now ends steeply: at 0.12 the cusps
// shut in ~20 ms, at the bottom of the 20-40 ms range, and only after the flow
// had stopped. Floating toward closure a little earlier in the deceleration,
// as the sinus vortices push them, they shut over ~27 ms as the flow ends.)
const OPEN_Q_AV = 0.2, OPEN_Q_SL = 0.12, OPEN_Q_AO = 0.2;   // flow (x peak) for a fully open valve
const _openCache = new WeakMap();
function openTrace(path) {
  const sum = hemoSummary(path);
  let tr = _openCache.get(sum);
  if (tr) return tr;
  const qmv = new Float32Array(OPEN_N), qao = new Float32Array(OPEN_N), qtv = new Float32Array(OPEN_N), qpa = new Float32Array(OPEN_N);
  let qtvMax = 0, qpaMax = 0;
  for (let i = 0; i < OPEN_N; i++) {
    const h = hemodynamics(i / OPEN_N, path);
    qmv[i] = h.Qmv; qao[i] = h.Qao; qtv[i] = h.Qtv; qpa[i] = h.Qpa;
    qtvMax = Math.max(qtvMax, h.Qtv); qpaMax = Math.max(qpaMax, h.Qpa);
  }
  const L = Math.round(OPEN_LEAD / OPEN_CYCLE * OPEN_N), Lo = Math.round(OPEN_LEAD_O / OPEN_CYCLE * OPEN_N), dt = OPEN_CYCLE / OPEN_N;
  const build = (q, qFull) => {
    const tgt = new Float32Array(OPEN_N), out = new Float32Array(OPEN_N);
    for (let i = 0; i < OPEN_N; i++) tgt[i] = 0.02 + 0.98 * smoothstep(0, qFull, Math.min(q[(i + Lo) % OPEN_N], q[(i + L) % OPEN_N]));
    let o = tgt[OPEN_N - 1];
    for (let pass = 0; pass < 2; pass++) for (let i = 0; i < OPEN_N; i++) {   // (2nd pass: periodic steady state)
      o += (tgt[i] - o) * (1 - Math.exp(-dt / (tgt[i] > o ? OPEN_TAU_OPEN : OPEN_TAU_CLOSE)));
      out[i] = o;
    }
    return out;
  };
  tr = { av: build(qmv, OPEN_Q_AV * sum.QmvMax), sl: build(qao, OPEN_Q_AO * sum.QaoMax),
    tv: build(qtv, OPEN_Q_AV * qtvMax), pv: build(qpa, OPEN_Q_SL * qpaMax) };
  // AV-valve events (10 % open): closure nearest the QRS, opening nearest
  // mid-cycle, for the mitral (left) and tricuspid (right) — rvPhase's anchors
  const cross = (arr, closing, near) => {
    let best = near, bd = 1;
    for (let i = 0; i < OPEN_N; i++) {
      const a = arr[(i + OPEN_N - 1) % OPEN_N], b = arr[i];
      if (closing ? a > 0.1 && b <= 0.1 : a <= 0.1 && b > 0.1) {
        const ph = i / OPEN_N; let d = Math.abs(ph - near); d = Math.min(d, 1 - d);
        if (d < bd) { bd = d; best = ph; }
      }
    }
    return best;
  };
  tr.warp = [[cross(tr.tv, true, 0), cross(tr.av, true, 0)], [cross(tr.tv, false, 0.5), cross(tr.av, false, 0.5)]];
  _openCache.set(sum, tr);
  return tr;
}
function valveOpening(name, phase, tr, path) {
  const arr = name === 'tricuspid' ? tr.tv : name === 'pulmonic' ? tr.pv : VALVES[name].opensInDiastole ? tr.av : tr.sl;
  const x = (((phase % 1) + 1) % 1) * OPEN_N, i0 = Math.floor(x) % OPEN_N, f = x - Math.floor(x);
  let o = arr[i0] + (arr[(i0 + 1) % OPEN_N] - arr[i0]) * f;
  // pathology tweaks
  if (path) {
    if (name === 'aortic' && path.aorticStenosis) o *= AS_OPEN;   // restricted opening
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
// LVOT cross-section (cm2) for the LVOT PW velocity: the drawn LVOT measures
// 2.15 cm in mid-systolic PLAX (verify-anatomy 'LVOT diameter'), pi (D/2)^2
export const LVOT_CSA = Math.PI * 1.075 * 1.075;
// mitral effective orifice area for the INSTANTANEOUS inflow velocity v = Qmv / A
// (MS: hemodynamics.js MS_AREA, Bernoulli): the E and A peaks, their ratio and the
// E-wave deceleration all come straight from the modelled transmitral flow
const MVA_NORM = 4.0;
// right-heart effective orifice areas (cm^2) for the instantaneous velocity
// v = Q / A, at a normal RV size (the tricuspid one grows with the dilated RV): the
// tricuspid inflow at the leaflet tips (normal E ~0.5 m/s, A ~0.3) and the RVOT
// just below the pulmonary valve (2.4 cm: peak ~0.95 m/s, VTI ~16 cm)
const A_TV = 6.0, A_RVOT = 4.5;
const ET_EJECT = 0.33, T_FILL = 0.37; // systolic ejection / diastolic filling time (s)
const K_AO = 1.56;                    // peak/mean profile factor (peak ≈ k·mean)

// Tricuspid regurgitation: present with path.tr, and as the functional TR of the
// secondary (post-capillary) pulmonary hypertension of severe MS. The peak (m/s)
// is the Bernoulli velocity of the peak systolic RV-RA pressure difference of
// the solved right heart (hemodynamics.js): ~4.3 in pulmonary hypertension (high
// PVR), ~3.4 in severe MS (PA pressure riding on the LA pressure), and LOWER in
// isolated severe TR (~2.5), its gradient blunted by the large RA v-wave.
export function trPeak(path = {}) {
  const RH = rightHeartLoad(path);
  return path.tr || RH.msPH ? hemoSummary(path).trVmax : 0;
}
// TR jet envelope (0..1) on the right heart's clock: holosystolic (TV closure at
// the QRS to TV opening ~0.50). With pressure overload it is the rounded,
// mid-systolic-peaking parabola of a high RV-RA gradient; in severe (torrential)
// TR it peaks EARLY (~30 % into systole) and decays almost linearly as the RA
// v-wave equalises the pressures (the triangular "V-wave cutoff" envelope).
export function trEnvelope(phaseRV, path = {}) {
  const round = pulse(0.02, 0.48, phaseRV);
  const w = path.tr && !path.rvpo ? rightHeartLoad(path).trVol : 0;
  if (!w) return round;
  const u = (((phaseRV % 1) + 1) % 1 - 0.02) / 0.46;
  const tri = u <= 0 || u >= 1 ? 0 : smoothstep(0, 0.27, u) * Math.pow(clamp((1 - u) / 0.73, 0, 1), 0.85);
  return round + (tri - round) * w;
}

// Build the per-frame haemodynamic snapshot the flow field reads. Returns a
// FRESH object (never the hemodynamics() singleton) so two geometryAt() results
// — e.g. echo's end-diastolic and end-systolic geometries — never alias.
function hemoFrame(hd, path, rvd0, phaseRV) {
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
    : s.QmvMax / (MVA_NORM * 100);
  // MR: the jet velocity IS the instantaneous LV-LA gradient (Bernoulli), so its
  // envelope rises through isovolumic contraction and is blunted in late systole
  // as the LA v-wave climbs (asymmetric in severe MR)
  const vMRPeak = Math.sqrt(Math.max(0, s.mrGradMax) / 4);
  // TR jet peak (m/s) encodes RV systolic pressure: high with pressure overload
  // (pulmonary hypertension), moderate for isolated TR. Feeds the PASP estimate.
  const vTRPeak = trPeak(path);
  // instantaneous 0..1 envelopes from the modelled flow waveform (carry E/A,
  // ejection and regurgitant timing) so velocityAt only multiplies scalars.
  const ejAo = s.QaoMax > 0 ? clamp(hd.Qao / s.QaoMax, 0, 1) : 0;
  // AS jet speed in time, vAS(t) (m/s): Bernoulli on the instantaneous LV to
  // ascending-aortic pressure drop while the valve is open, so the CW envelope
  // traced by the learner gives the circulation's mean gradient (mean/peak
  // gradient ~0.6) as well as its peak. The jet's spatial profile multiplies it.
  const vAS = path.aorticStenosis && hd.Qao > 0.5 ? Math.min(vAoPeak, Math.sqrt(Math.max(0, hd.Plv - hd.Pao) / 4)) : 0;
  const ejAS = vAoPeak > 0 ? vAS / vAoPeak : 0;
  // LVOT PW velocity (m/s): the solved aortic flow over the LVOT area, so its
  // VTI x LVOT area is the stroke volume
  const vLvot = hd.Qao > 0 ? hd.Qao / (LVOT_CSA * 100) : 0;
  const ejMv = s.QmvMax > 0 ? clamp(hd.Qmv / s.QmvMax, 0, 1) : 0;
  const ejMr = hd.Qmr > 0.5 && !(hd.Qmv > 0) && vMRPeak > 0 ? clamp(Math.sqrt(Math.max(0, hd.Plv - hd.Pla) / 4) / vMRPeak, 0, 1) : 0;
  // right-heart velocities from the SOLVED right heart (review 12, MOT-4/PAT-10):
  // tricuspid inflow v = Qtv / A_TV and RV outflow v = Qpa / A_RVOT, so the
  // tricuspid E/A, the RVOT acceleration and ejection times, the PH mid-systolic
  // notch and the ASD's raised Qp all come from the circulation. (Each respiratory
  // phase of the tamponade case is its own solved circulation, so its tricuspid
  // E rises in inspiration by itself.) ejTv / ejPv are the same flows over the
  // case's peak, for the colour-flow and 3D streams that scale an envelope.
  const RVD = rvDrive(path);
  // (the tricuspid annulus dilates with the RV — the case's RV, not one breath's:
  // rvd0, the apnoeic drive of a tamponade breath, from geometryAt; the RVOT just
  // below the valve is held at the pulmonary annulus, as drawn — so a raised Qp
  // shows as velocity)
  const aTv = A_TV * (rvd0 || RVD).rvSize, aRvot = A_RVOT;
  const vTv = Math.max(0, hd.Qtv) / (aTv * 100), vPv = Math.max(0, hd.Qpa) / (aRvot * 100);
  const vTvPeak = (s.QtvMax || 0) / (aTv * 100);
  const ejTv = vTvPeak > 0 ? clamp(vTv / vTvPeak, 0, 1) : 0;
  const ejPv = RVD.QpaMax > 0 ? clamp(hd.Qpa / RVD.QpaMax, 0, 1) : 0;
  // INSTANTANEOUS regurgitant jet speeds (m/s). MR: Bernoulli on the LV-LA
  // pressure drop while the circulation carries regurgitant flow, so the CW
  // envelope is holosystolic, its upstroke follows dP/dt and the compliant LA's
  // v-wave blunts it late in systole. TR: the case's peak velocity on the right
  // heart's clock, shaped by trEnvelope (rounded with pressure overload, early-
  // peaking and triangular when torrential TR equalises RV and RA pressures).
  // (second guard, review 12 PHY-1: never while forward flow crosses the valve)
  const vMR = hd.Qmr > 0 && !(hd.Qmv > 0) ? Math.sqrt(Math.max(0, hd.Plv - hd.Pla) / 4) : 0;
  const ejTr = vTRPeak > 0 ? trEnvelope(phaseRV != null ? phaseRV : hd.phase, path) : 0;
  // (no TR jet while forward flow still crosses the open tricuspid valve)
  const vTR = hd.Qtv > 0.5 ? 0 : vTRPeak * ejTr;
  return {
    Plv: hd.Plv, Pao: hd.Pao, Vlv: hd.Vlv, Pla: hd.Pla,
    Qao: hd.Qao, Qmv: hd.Qmv, Qmr: hd.Qmr, rho: hd.rho, kv: hd.kv,
    EDV: s.EDV, ESV: s.ESV, SV: s.SV, EF: s.EF, forwardSV: s.forwardSV,
    regurgFraction: s.regurgFraction, gradient: s.gradient,
    PlvSys: s.PlvSys, PaoSys: s.PaoSys, PaoDia: s.PaoDia,
    vAoPeak, vMitPeak, vMRPeak, vTRPeak, vAS, ejAS, vLvot, ejAo, ejMv, ejMr, ejTv, ejPv, ejTr, vMR, vTR,
    vTv, vPv, Qtv: hd.Qtv, Qpa: hd.Qpa, aRvot, aTv,
    // shunt jets (m/s, + = left-to-right): Bernoulli on the modelled LA-RA / LV-RV
    // pressure difference, so the ASD signal is phasic and the VSD jet follows the
    // interventricular gradient
    vAsd: hd.vAsd || 0, vVsd: hd.vVsd || 0, vVsdPeak: s.vVsdPeak || 0, QpQs: s.QpQs || 1,
    // solved right heart: instantaneous RV / RA / PA and pericardial pressures,
    // and the case's RV systolic and mean RA pressure
    Prv: hd.Prv, Pra: hd.Pra, Ppa: hd.Ppa, Ppc: hd.Ppc || 0, RVSP: s.RVSP, RAP: s.RAP,
    // solved venous inflow (mL/s, + toward the atrium): the caval return into the
    // RA and the pulmonary-venous inflow into the LA. (The right heart is solved on
    // its own activation clock, so the caval flow at the cycle phase IS on the
    // right heart's clock.) The hepatic- and pulmonary-vein PW traces are these
    // flows over a vessel area.
    Qvr: hd.Qvr || 0, Qpv: hd.Qpv || 0,
    // solved regurgitant flows (mL/s) over their cycle peaks, the TR volume per
    // beat (mL) and the effective TR orifice (cm2): the colour jets are sized by
    // the instantaneous regurgitant flow, not by a grade label (review 12, PAT-4/PAT-11)
    AVA, QmrMax: s.QmrMax || 0, Qtr: hd.Qtr || 0, QtrMax: s.QtrMax || 0, trRegurgVol: s.trRegurgVol || 0, Atr: s.Atr || 0,
  };
}

// ---------------------------------------------------------------------------
// Live geometry for a given phase + pathology. Returns per-chamber ellipsoids
// (centre + radii) plus wall thicknesses, so both renderers stay consistent.
// ---------------------------------------------------------------------------
// Reference (normal) end-diastolic volume, for scaling remodelled ventricles.
let _edvRef = null;
const edvRef = () => (_edvRef != null ? _edvRef : (_edvRef = hemoSummary({}).EDV));

// The drawn right ventricle from the SOLVED right heart (review 12, MOT-4):
//  - size: the RV is scaled so its drawn end-diastolic volume holds the solved
//    RVEDV (cube root of the volume ratio to normal, the RV's dilatation);
//  - contractility: its short-axis (FAC) share follows the solved RV ejection
//    fraction relative to normal (rvC);
//  - long-axis share (TAPSE, rvL): afterload hits the longitudinal fibres
//    hardest, so above a normal RV systolic pressure it falls further (TAPSE
//    ~1.5 cm in severe PAH with an FAC of ~30-35 %);
//  - timing: kRV = (RVEDV - Vrv) / RV stroke volume, the solved RV volume.
// RV_DRAW holds the calibration so the drawn RV voxel volume and EF track the
// solved ones as the LV's track rho (verify-anatomy, 'Right-heart drive').
const RV_DRAW = { longGain: 1, gamma: 1.55, lvExp: 0.8, lvhK: 0.36, septK: 0.55, wallK: 0.37,
  afterloadRef: 35, afterloadK: 0.0115, lonMin: 0.3, cMin: 0.2, cMax: 1.3 };
let _rvRef = null;
const _rvDrive = new Map();
function rvDrive(path) {
  const s = hemoSummary(path);
  let r = _rvDrive.get(s);
  if (r) return r;
  if (!_rvRef) _rvRef = hemoSummary({});
  const ef = clamp(s.RVEF / _rvRef.RVEF, RV_DRAW.cMin, RV_DRAW.cMax);
  const lon = clamp(1 - RV_DRAW.afterloadK * Math.max(0, s.RVSP - RV_DRAW.afterloadRef), RV_DRAW.lonMin, 1);
  // (the volume-overload septum, flattened in diastole and recovering in systole,
  // and a hypertrophied septum thickening into the RV both add to the drawn RV's
  // ejection: the free wall's share is lowered by as much)
  const sep = rightHeartLoad(path).septum, swing = sep ? Math.max(0, sep(0) - sep(1)) : 0;
  const lvh = path.lvh || path.aorticStenosis ? 0.55 : 0;
  r = {
    rvSize: Math.cbrt(s.RVEDV / _rvRef.RVEDV),
    rvC: Math.max(0.05, Math.pow(ef, RV_DRAW.gamma) - RV_DRAW.septK * swing - RV_DRAW.wallK * lvh),
    rvL: RV_DRAW.longGain * Math.min(1, ef) * lon,
    RVEDV: s.RVEDV, RVSV: Math.max(1e-6, s.RVEDV - s.RVESV), QpaMax: 0,
  };
  for (let i = 0; i < 240; i++) r.QpaMax = Math.max(r.QpaMax, hemodynamics(i / 240, path).Qpa);
  _rvDrive.set(s, r);
  return r;
}

// ---- RWMA: compensatory hyperkinesis of the remote walls -------------------
// The circulation's LVEF is that of the whole ventricle; the anatomy then
// freezes the infarcted territory at its diastolic position (rwmaBlend), so on
// its own the drawn cavity would eject less than the PV loop says. The remote
// (non-infarcted) myocardium is hyperkinetic after an infarct: here it shortens
// by kappa times the global volume change, rho' = 1 - kappa (1 - rho), with
// kappa solved once per case so that the DRAWN cavity (voxel integration of the
// rendered LV lumen, clipped at the mitral and aortic annular planes as the
// audit does) ejects exactly the circulation's EF. Only the anatomy sees rho';
// the coarse proxies (lv.r: EF/LVID readouts, the Simpson truth) stay on rho,
// so the panel, the report and the image all carry one EF. Note: most of that
// agreement comes from the per-territory Emax cut in the circulation (the
// solved EF fell to what the drawn LV already ejected), not from visible remote
// hyperkinesis: the A4C end-systolic LV area moved by only ~2 % (review 12).
const _kappaCache = new Map();
function drawnLvVolume(A, h) {
  const M = A.valves.mitral.c, Av = A.valves.aortic.c, U = LM.U_AO;
  let n = 0;
  for (let x = -4; x <= 4; x += h) for (let z = -4; z <= 4; z += h) for (let y = -9; y <= 1.5; y += h) {
    if (lumenDist(x, y, z, A, 'LV') >= 0) continue;
    const ua = (x - Av[0]) * U[0] + (y - Av[1]) * U[1] + (z - Av[2]) * U[2];
    if (ua > 0 && (y > M[1] + 0.02 || y - M[1] > -0.3)) continue;
    n++;
  }
  return n * h * h * h;
}
export function rwmaRemoteBoost(path = {}) {
  if (!path.rwma) return 1;
  const key = JSON.stringify(path);
  let kappa = _kappaCache.get(key);
  if (kappa != null) return kappa;
  const s = hemoSummary(path), h = 0.14;
  const target = drawnLvVolume(geometryCore(0, path, 1).A, h) * s.ESV / s.EDV;
  const vAt = (kp) => drawnLvVolume(geometryCore(s.tMinVol, path, kp).A, h);
  // the drawn ES volume falls almost linearly with kappa: a few secant steps
  let k0 = 1, v0 = vAt(k0), k1 = 1.4, v1 = vAt(k1);
  for (let i = 0; i < 6 && Math.abs(v1 - target) > 0.3 && v1 !== v0; i++) {
    const k2 = clamp(k1 + (target - v1) * (k1 - k0) / (v1 - v0), 1, 2.5);
    k0 = k1; v0 = v1; k1 = k2; v1 = vAt(k1);
  }
  kappa = k1;
  _kappaCache.set(key, kappa);
  return kappa;
}

export function geometryAt(phase, path = {}) {
  return geometryCore(phase, path, path.rwma ? rwmaRemoteBoost(path) : 1);
}

function geometryCore(phase, path, kappa) {
  const kick = atrialKick(phase);

  // PHYSICS COUPLING: the lumped-parameter circulation solves the LV pressure–
  // volume loop. hd.rho = Vlv/EDV is the ABSOLUTE cavity-volume ratio (1 at
  // end-diastole, ESV/EDV at end-systole), so it carries the true ejection
  // fraction (low in DCM/ischemia). The cavity semi-axes are scaled so their
  // product tracks rho EXACTLY (sShort²·sLong = rho), with the long axis
  // shortening less than the short axis; hence the measured EF/LVIDs fall out of
  // the modelled volume rather than a prescribed phase curve. hd.kv (normalised)
  // still times the wall-thickening. The right ventricle follows its own solved
  // volume (rvDrive: kRV from Vrv, size from RVEDV, contractility from RVEF).
  const phaseRV = rvPhase(phase, path);
  const otr = openTrace(path);        // (first: it samples the hemodynamics() singleton)
  const RVD = rvDrive(path);
  // (tamponade: each respiratory phase is its own solved circulation, so the RV
  // fills more in inspiration; the tricuspid orifice the inflow velocity is read
  // through is the apnoeic case's — a breath does not resize the annulus. Before
  // hemodynamics(), as rvDrive samples its singleton.)
  const RVD0 = path.resp ? rvDrive({ ...path, resp: undefined }) : RVD;
  const hd = hemodynamics(phase, path);
  const kvR = (RVD.RVEDV - hd.Vrv) / RVD.RVSV, kR = clamp(kvR, 0, 1);
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
  // (severe AS: the stiff, hypertrophied LV's diastolic dysfunction enlarges the
  // LA too — LA AP ~4-4.5 cm — although its mean pressure is modelled near normal)
  // (the chronic DISTENDING pressure: in tamponade the pericardial pressure is subtracted)
  const laScale = (1 + 0.045 * Math.max(0, (hd.sum.Pla0 || 8) - (hd.sum.Ppc || 0) - 8)) * (path.aorticStenosis ? 1.2 : 1);
  let lvWallMul = 1.0;
  if (path.dilated) { lvWallMul = DCM_SEVERITIES[dcmSeverityOf(path)].wallMul; }
  if (path.lvh || path.aorticStenosis) { lvWallMul = 1.55; } // concentric hypertrophy (IVSd ~1.4)

  // Long axis shortens less than the short axis. aLong = 0.24 gives, for the
  // normal EF, MAPSE ~1.55 cm, GLS ~-19 % and fractional shortening ~28 % on
  // the 8.3 cm ventricle — all mid-normal. Concentric LVH (severe AS, the
  // hypertrophied heart) loses longitudinal (subendocardial) function first while
  // EF is kept by radial thickening: aLong 0.17 gives MAPSE ~1.1 cm and GLS
  // ~-15 % at a normal EF (reduced GLS with preserved EF in AS — a prognostic
  // marker; e.g. Ng et al., Eur Heart J Cardiovasc Imaging 2011; Magne et al.,
  // JACC Cardiovasc Imaging 2019).
  const aLong = path.aorticStenosis || path.lvh ? 0.17 : 0.24;
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

  // (the drawn RV shares its septum with the LV: a larger LV bows into it and a
  // hypertrophied septum takes up its cavity, so the RV is scaled up to hold the
  // solved RVEDV — the volume-exact coupling the LV has through rho)
  const rvFit = Math.cbrt(Math.pow(lvScale, -RV_DRAW.lvExp) / (1 - RV_DRAW.lvhK * Math.max(0, lvWallMul - 1)));
  const rvC = RVD.rvC;                    // RV contractility (solved RVEF, rvDrive)
  // (the coarse RV proxy, for flow gating, at the drawn RV's size and contraction)
  const rvR = BASE.rv.r.map((r, i) => {
    const shorten = (i === 1 ? 0.12 : 0.26) * rvC;
    return r * Math.min(1, RVD.rvSize * rvFit) * (1 - kR * shorten);
  });
  const rvWall = BASE.rv.wall * (1 + kR * 0.5);

  // anatomical SDF bundle (built once here so the chordae apparatus below can
  // hang off the same apex-anchored papillary tips + live valve planes).
  // (RWMA: the anatomy's remote walls shorten by kappa x the global change; see
  // rwmaRemoteBoost — the proxies above stay on rho)
  const rhoA = kappa === 1 ? rho : Math.max(0.05, 1 - kappa * (1 - rho));
  const sLongA = kappa === 1 ? sLong : Math.pow(rhoA, aLong);
  const sShortA = kappa === 1 ? sShort : Math.pow(rhoA, (1 - aLong) / 2);
  const lvWallA = kappa === 1 ? lvWall : BASE.lv.wall * lvWallMul * clamp(1 / sShortA, 1, 1.7);
  const A = anatomyParams(k, kick, path, { sShort: sShortA, sLong: sLongA, lvWall: lvWallA, lvScale, lvLong, laScale, kRV: kR, phaseRV, laAmp: laPhasic(path),
    rvSize: RVD.rvSize * rvFit, rvC: RVD.rvC, rvL: RVD.rvL, laTiming: { mvo: hd.sum.mvoPhase, aAt: hd.sum.aAt } }, phase);

  // atria: coarse axis-aligned proxies of the live anatomical atria (which sit on
  // their moving annuli and follow the reservoir/conduit/booster curve), so flow
  // gating, labels and measurements stay locked to the rendered atrial wall.
  const laR = [A.la.r1, A.la.rl, A.la.r2];
  const raR = [A.ra.r1, A.ra.rl, A.ra.r2];
  const valves = {
    mitral: valveOpening('mitral', phase, otr, path),
    tricuspid: valveOpening('tricuspid', phase, otr, path),
    aortic: valveOpening('aortic', phase, otr, path),
    pulmonic: valveOpening('pulmonic', phase, otr, path),
  };

  // rheumatic MS: the centre of the open funnel tip (the fish-mouth orifice), where
  // PSAX-MV planimetry is done — the view tracks it (views.js)
  if (path.mitralStenosis) {
    const Gt = { A, valves: { mitral: 0.35 }, path };
    const ta = avTipWorld('mitral', Gt, 1, 0, MS_TUBE / 2), tp = avTipWorld('mitral', Gt, -1, 0, MS_TUBE / 2);
    A.msTip = [(ta[0] + tp[0]) / 2, (ta[1] + tp[1]) / 2, (ta[2] + tp[2]) / 2];
  }

  return {
    phase, k, kick, kv, phaseRV, kRV: kR,
    // instantaneous lumped-parameter circulation state (PV loop, pressures,
    // stroke volume, EF, jet peak velocities) — read by the echo measurements,
    // the flow field and the UI PV-loop panel. A fresh snapshot (not the
    // hemodynamics() singleton) so distinct geometries never alias.
    hemo: hemoFrame(hd, path, RVD0, phaseRV),
    // Coarse per-chamber ellipsoid proxies — used by the velocity field, the
    // on-image labels and the LVIDd/EF measurements. Tissue classification and
    // the 3D surface use the anatomical SDF (G.A) below, not these.
    lv: { c: [BASE.lv.c[0], lvCy, BASE.lv.c[2]], r: lvR, wall: lvWall },
    // the same proxy at end-diastolic size: the filling stream is gated to it, so a
    // fixed Doppler gate in the inflow path sees the whole E wave from mitral
    // opening (when the cavity is still near end-systolic size), not only its tail
    lvEd: { c: [BASE.lv.c[0], lvApexY + lvR1d, BASE.lv.c[2]], r: [BASE.lv.r[0] * lvScale, lvR1d, BASE.lv.r[2] * lvScale] },
    rv: { c: BASE.rv.c, r: rvR, wall: rvWall },
    rvEd: { c: BASE.rv.c, r: BASE.rv.r },   // (end-diastolic size: tricuspid inflow gate, as lvEd)
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
    for (const sgn of sides) for (const q of [0.3, 0.72]) out.push(avTipWorld(name, G, sgn, s * q, name === 'mitral' && G.path && G.path.mitralStenosis ? MS_TUBE : 0));
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
  // One GLS per case: the territory pattern sets the segments' SHAPE, and a
  // global factor scales them so their mean is the long-axis shortening of the
  // same solved cavity (longitudinalStrain, the panel's GLS). DCM, LVH/AS and
  // RWMA then show the same GLS in the panel, the bullseye and the report.
  const rw = buildRwma(path);
  const w = AHA_SEGMENTS.map((s) => segAffected(s, rw));
  const shape = w.reduce((a, x) => a + NORMAL_SEG_STRAIN * (1 - x), 0) / w.length;
  const ls = longitudinalStrain(path);
  const gf = ls != null && shape < 0 ? ls / shape : 1;   // global contractility factor
  const segments = AHA_SEGMENTS.map((s, i) => (
    { name: s.name, angle: s.angle, level: s.level, strain: NORMAL_SEG_STRAIN * gf * (1 - w[i]) }));
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

// Tamponade: the pericardial pressure exceeds the transmural (distending) pressure
// of the thin-walled right heart when its chambers are least full — the RA free
// wall inverts from late diastole into systole (RA transmural < 0) and the RV free
// wall is pushed in in early diastole (RV transmural < ~1 mmHg). Each free wall is
// indented by a depth that follows that phase's solved pressures, over a Gaussian
// patch centred on the free wall (away from the septum and the AV groove).
const RA_INV_K = 0.8, RA_INV_MAX = 1.0;      // cm per mmHg below zero transmural, cap
const RV_COL_K = 1.2, RV_COL_MAX = 0.55;     // cm per mmHg below the buckling pressure
const RV_BUCKLE_P = 1.0;                     // mmHg (hemodynamics.js RV_BUCKLE)
const WALL_IN = 0.3;                         // drawn atrial / RV free-wall thickness inside the lumen surface (cm)
const _tampCache = new WeakMap();
function tamponadeDent(G) {
  let t = _tampCache.get(G);
  if (t) return t;
  const A = G.A, H = G.hemo || {};
  const raT = (H.Pra || 0) - (H.Ppc || 0), rvT = (H.Prv || 0) - (H.Ppc || 0);
  const dRa = clamp(-raT * RA_INV_K, 0, RA_INV_MAX), dRv = clamp((RV_BUCKLE_P - rvT) * RV_COL_K, 0, RV_COL_MAX);
  // free-wall points: march out of the lumen from the chamber centre
  // (into the lumen first, then out through the free wall)
  const wallOf = (c, u, which) => {
    let q = c, inside = false;
    for (let k = 0; k < 240; k++) {
      q = [q[0] + u[0] * 0.05, q[1] + u[1] * 0.05, q[2] + u[2] * 0.05];
      const d = lumenDist(q[0], q[1], q[2], A, which);
      if (d < 0) inside = true; else if (inside) break;
    }
    return q;
  };
  const uRa = LM.IAS_N;
  let uRv = [A.rv.c[0] - A.lv.c[0], 0, A.rv.c[2] - A.lv.c[2]];
  const L = Math.hypot(uRv[0], uRv[2]) || 1; uRv = [uRv[0] / L, 0, uRv[2] / L];
  // (the RA free wall itself: the crista terminalis is a fold carved from the
  // lumen, and the march must not stop on it)
  t = { dRa, dRv, uRa, uRv, wRa: dRa > 0 ? wallOf(A.ra.c, uRa, 'RA_NOFOLD') : null, wRv: dRv > 0 ? wallOf([A.lv.c[0], A.rv.c[1], A.lv.c[2]], uRv, 'RV') : null };
  _tampCache.set(G, t);
  return t;
}
// the free-wall indentation of a tamponade geometry (cm): { ra, rv, wRa, wRv, uRa, uRv }
export function tamponadeCollapse(G) {
  const t = tamponadeDent(G);
  return { ra: t.dRa, rv: t.dRv, wRa: t.wRa, wRv: t.wRv, uRa: t.uRa, uRv: t.uRv };
}
// classify override inside a dented free wall: null (no change), or a tissue
function tamponadeWall(px, py, pz, G, code) {
  const t = tamponadeDent(G), A = G.A;
  for (const [w, d, which, lum, sig] of [[t.wRa, t.dRa, 'RA', BODY.RA, 1.8], [t.wRv, t.dRv, 'RV', BODY.RV, 1.6]]) {
    if (!w) continue;
    const dx = px - w[0], dy = py - w[1], dz = pz - w[2], r2 = dx * dx + dy * dy + dz * dz;
    if (r2 > 9 * sig * sig) continue;
    const dent = d * Math.exp(-r2 / (sig * sig));
    if (dent < 0.05) continue;
    if (code !== lum && code !== BODY.MYO && code !== BODY.PAP && code !== BODY.OUT && code !== BODY.PERI && code !== BODY.FAT) continue;
    const ld = lumenDist(px, py, pz, A, which);
    if (ld > 0.45) continue;
    // (the drawn blood-pool edge sits ~WALL_IN inside the lumen surface, the
    // thin wall between it and the surface: both move in by the dent)
    const sh = ld + dent;                          // the lumen surface pushed in
    if (sh < -WALL_IN) { if (code === lum) return null; continue; }
    if (sh < -WALL_IN + 0.22) return TISSUE.MYO;                      // the inverted wall
    if (ld < 0.25) return TISSUE.PERICARDIUM;                         // fluid fills the dent
  }
  return null;
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
  if (path.tamponade) {
    const tw = tamponadeWall(px, py, pz, G, body.code);
    if (tw === TISSUE.MYO) return cls(TISSUE.MYO, 0.55);
    if (tw === TISSUE.PERICARDIUM) return cls(TISSUE.PERICARDIUM, 0.02);
  }
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
  // rheumatic MS: both thickened leaflets dome, so the fused orifice sits anterior
  // to the normal (posterior) coaptation line — central in the LV short axis
  const uc = uP + (mod && mod.ms ? 0.75 : P.coapt) * w;
  // functional (tethered) regurgitation: the papillary muscles pull the
  // coaptation point further into the ventricle (DCM mitral tenting)
  const tent = P.tent + (mod && mod.tent ? mod.tent : 0);
  const ac = ah + tent;
  const uh = sgn > 0 ? uA : uP;
  const reach = sgn > 0 ? uA - uc : uc - uP;               // horizontal reach to the coaptation line
  const Lb = Math.sqrt(reach * reach + tent * tent);
  const ovl = P.ovl * half + 0.08;
  let oo = o;
  if (mod && mod.ms && q2 > mod.msQ2) oo *= 0.12;           // rheumatic commissural fusion
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
    // (commissures fused beyond |q| = MS_Q), opening with the transmitral flow —
    // the central, medio-laterally elongated "fish-mouth" of the short axis
    const on = o / 0.35 > 1 ? 1 : o / 0.35;
    const e = 1 - q2 / mod.msQ2;
    // (the gap is solved at the live annular size: the orifice area stays the MVA
    // while the annulus breathes with the LV)
    const hg = e > 0 ? 0.5 * msGapAt(mod.msA, sc) * Math.sqrt(e) * on : 0;
    // The fused tips form one PLANAR ring across both leaflets and the fused
    // commissures (the funnel tip, perpendicular to the inflow axis), so a short
    // axis cut at the tips shows a closed fish-mouth rather than a broken C.
    // (sheared along c so the slit runs closer to medial-lateral in the short axis)
    _lp.ut = uc + sgn * hg + MS_SHEAR * q * P.Ric * sc;
    _lp.at = MS_DEPTH * sc * (0.35 + 0.65 * on);
    // Doming: the thin, still mobile body bellies into the LV (the anterior leaflet
    // toward the LVOT/septum) while the fused tip is held back at the orifice, so
    // the leaflet bends at a knee ~0.6 of the way down: the diastolic "hockey
    // stick" of the long-axis views. The knee stands MS_BOW off the hinge-tip
    // chord (anterior; the stiffer posterior leaflet domes less).
    const du = _lp.ut - uh, da = _lp.at - ah, Lc = Math.sqrt(du * du + da * da) || 1e-9;
    const bow = (sgn > 0 ? MS_BOW_A : MS_BOW_P) * on;
    _lp.ub = uh + du * 0.6 + sgn * bow * da / Lc;
    _lp.ab = ah + da * 0.6 - sgn * bow * du / Lc;
    // subvalvular fusion: below the tip ring the fused tips and matted chordae
    // continue as a short tunnel of the same lumen (avLeafletHit)
    _lp.tube = MS_TUBE * (0.4 + 0.6 * on);
  } else _lp.tube = 0;
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
    if (sgn < 0 && P.scallops && !(mod && mod.ms)) edgeCut = 0.9 + 0.1 * Math.abs(Math.cos(P.scallops * Math.PI * 0.5 * (q + 1)));
    const d1 = seg2d(u, a, L.uh, L.ah, L.ub, L.ab);
    const t1 = _segT;
    const d2 = seg2d(u, a, L.ub, L.ab, L.ut, L.at);
    const t2 = _segT;
    if (L.tube > 0 && seg2d(u, a, L.ut, L.at, L.ut, L.at + L.tube) <= (sgn > 0 ? P.thickA : P.thickP) * thickMul * (1.1 + MS_TIP_BULB)) return 1;
    const bodyFrac = 1 - (P.ovl + 0.08) / L.len;
    const f = d1 <= d2 ? t1 * bodyFrac : bodyFrac + t2 * (1 - bodyFrac);
    if (f > edgeCut) continue;
    const th0 = (sgn > 0 ? P.thickA : P.thickP) * thickMul;
    // thin membrane tapering to the free edge; rheumatic tips are bulbous
    const th = mod && mod.ms ? th0 * (1.1 + MS_TIP_BULB * f * f * f) : th0 * (1.05 - 0.4 * f);
    if ((d1 < d2 ? d1 : d2) <= th) return f;
  }
  return -1;
}
// Mitral pathology modifiers: rheumatic stenosis (commissural fusion, doming),
// and a posterior (P2) prolapse whose extent grows with the MR grade.
const PROLAPSE = { mild: 0, moderate: 0.55, severe: 1 };
// tricuspid malcoaptation gap (cm) in TR, by grade
const TV_MOD = { mild: { gap: 0.1 }, moderate: { gap: 0.25 }, severe: { gap: 0.45 } };
// Rheumatic MS leaflet shape (review 12, PAT-1 / PAT-5).
// - MS_Q: open commissural span (|q| < MS_Q; fused beyond).
// - Thickness: the body is ~1.2x a normal leaflet and still pliable, the
//   thickening is at the free edge: half-thickness th0 * (1.1 + MS_TIP_BULB f^3)
//   along the leaflet (f = 0 hinge .. 1 tip), i.e. ~2 mm bodies and ~4 mm tips
//   (they were 1.9x thick, ~6 mm bulbs, which bloomed over the orifice).
// - MS_DEPTH: the tip ring below the annulus centre. Below it the fused tips and
//   matted chordae continue for MS_TUBE as a short tunnel of the same lumen
//   (subvalvular fusion), so a short-axis cut a few degrees off the ring plane
//   still shows the limiting orifice as one closed fish-mouth.
// - MS_BOW_A / MS_BOW_P: diastolic doming, the knee's offset from the hinge-tip
//   chord (anterior leaflet: the "hockey stick"; the posterior leaflet is stiffer).
const MS_Q = 0.66, MS_TIP_BULB = 1.0, MS_DEPTH = 1.4, MS_THICK = 1.2, MS_TUBE = 0.45;
const MS_BOW_A = 0.45, MS_BOW_P = 0.2;
const MS_CALC_F = 0.85;   // leaflet fraction (0 hinge .. 1 tip) beyond which the leaflet is calcified
const MS_SHEAR = 0.25;                                    // tip-ring shear along the commissural axis
const MS_TIP_T = 0.5 * (AVL.mitral.thickA + AVL.mitral.thickP) * MS_THICK * (1.1 + MS_TIP_BULB);   // tip half-thickness
// Planimetry is done on the IMAGE: the sonographer traces the inner border of the
// bright leaflet ring, and that border lies inside the anatomical lumen by about
// the half-width of the point-spread function (the leaflets bloom into the
// orifice; high gain makes it worse: Martin RP et al., Am J Cardiol 1979;43:560).
// MS_AREA is the planimetered area (the Hakki and PHT areas are tied to it), so
// the lumen is built MS_RIM wider all round than the MS_AREA ellipse. MS_RIM is
// calibrated on the rendered PSAX-MV image (tools/verify-anatomy.mjs, 'B-mode
// planimetry'): the near-black orifice (grey <= 45) is ~0.85 cm2 and the orifice
// traced at the mid-grey edge ~1.1 cm2 for the severe valve's 0.95 cm2, in an
// opening of ~1.6 x 0.75 cm. The tissue-label lumen is then ~1.4 cm2.
const MS_RIM = 0.095;
// AP gap (cm) between the tip centrelines at annular scale sc for a planimetered
// orifice msA: an ellipse whose semi-axes are those of the lumen less MS_RIM. The
// lumen's long semi-axis is the open commissural span less the fused end (~0.7 of
// a tip half-thickness); its short semi-axis is half the gap less the tip
// half-thickness.
function msGapAt(msA, sc) {
  const halfL = MS_Q * AVL.mitral.Ric * sc - 0.7 * MS_TIP_T;
  return 2 * (msA / (Math.PI * (halfL - MS_RIM)) + MS_RIM + MS_TIP_T);
}
// The severe-MS anterior/posterior leaflet in its own long-axis section (q = 0),
// for the audit: hinge, knee and tip as [u (toward the aortic root), a (into the
// LV)] in cm, and the body / tip half-thickness.
export function msLeafletProfile(G, sgn = 1, q = 0) {
  const path = G.path || { mitralStenosis: true }, P = AVL.mitral;
  const L = avLeafletPoly(P, sgn, q, G.valves.mitral, G.A.lv.mvR / LM.MV_R, avMod(path));
  const th0 = (sgn > 0 ? P.thickA : P.thickP) * MS_THICK;
  const bodyFrac = 1 - (P.ovl + 0.08) / L.len, fm = 0.5 * bodyFrac;
  return { hinge: [L.uh, L.ah], knee: [L.ub, L.ab], tip: [L.ut, L.at], tube: L.tube, flow: LEAFLETS.mitral.flow,
    thBody: th0 * (1.1 + MS_TIP_BULB * fm * fm * fm), thTip: th0 * (1.1 + MS_TIP_BULB) };
}
const _modCache = new WeakMap();
function avMod(path) {
  let m = _modCache.get(path);
  if (!m) {
    // MS orifice: fused commissures leave a fish-mouth ~2 cm wide (|q| < MS_Q).
    // Its AP opening (msGapAt) is set so the orifice INSIDE the thickened tips (the
    // blood area a sonographer planimeters) matches the grade's MVA.
    const msA = MS_AREA[gradeOf(path)] || MS_AREA.severe;
    m = { ms: !!path.mitralStenosis, msQ2: MS_Q * MS_Q, msA,
      prolapse: path.mr ? PROLAPSE[gradeOf(path)] : 0,
      tent: path.dilated && !path.mr ? 0.55 : 0 };
    _modCache.set(path, m);
  }
  return m;
}
// World-space free-edge (tip) point of a leaflet at intercommissural fraction q.
function avTipWorld(name, G, sgn, q, da = 0) {
  const P = AVL[name], lf = LEAFLETS[name], vc = G.A.valves[name].c;
  const o = G.valves[name], sc = name === 'mitral' ? G.A.lv.mvR / LM.MV_R : 1;
  const mod = name === 'mitral' ? avMod(G.path || {}) : null;
  const L = avLeafletPoly(P, sgn, q, o, sc, mod);
  const c = q * P.Ric * sc;
  const t = [
    vc[0] + lf.ap[0] * L.ut + lf.flow[0] * (L.at + da) + lf.com[0] * c,
    vc[1] + lf.ap[1] * L.ut + lf.flow[1] * (L.at + da) + lf.com[1] * c,
    vc[2] + lf.ap[2] * L.ut + lf.flow[2] * (L.at + da) + lf.com[2] * c,
  ];
  // the mitral apparatus rides the angle-dependent annular excursion
  if (name === 'mitral') t[1] += mitralLift(t[0], t[1], t[2], G.A);
  else if (name === 'tricuspid') t[1] += tricuspidLift(t[0], t[1], t[2], G.A);
  return t;
}

// Semilunar cusp shape (valveTissue): shut-cusp coaptation point downstream of
// the hinge line (cm), belly sag toward the ventricle (cm), height of the
// coaptation zone above the coaptation point at the centre (cm), samples of the
// blended (opening / closing) profile.
const SL_COAPT = 0.28, SL_SAG = 0.13, SL_COAPT_H = 0.55, SL_N = 16;
const AS_OPEN = 0.28;                                     // AS: restricted opening (valveOpening)
// aortic cusps: nodulus of Arantius (cm) and the half-width of the shut cusps'
// coaptation seams (cm: two apposed lunulae, ~0.5 mm — finer than the renderer's
// coarsest sampling block, so its coarse sampler looks for valve tissue inside
// each block, see valveAt); calcific AS: the nodule threshold of the calcium
// lattice, and the half-width of the open cusps' commissural slits (cm)
const AO_NODULE = 0.05, AO_COAPT_W = 0.025, AS_NODULE = 0.06, AS_SLIT = 0.06, AS_SLIT_T = 0.3;
const _slA = new Float64Array(SL_N + 1), _slR = new Float64Array(SL_N + 1);
// Valve tissue alone at a point (classify() tests the valves before any other
// tissue, so a hit here is what classify returns there): the renderer's coarse
// sampler uses it to look for sub-block cusps cheaply. Returns the shared result
// object (tissue / echo / calc) or null.
export function valveAt(px, py, pz, G, path = {}) { return valveTissue(px, py, pz, G, path); }
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
    const dy = py - vc[1] - (name === 'mitral' ? mitralLift(px, py, pz, G.A) : name === 'tricuspid' ? tricuspidLift(px, py, pz, G.A) : 0);
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
      const f = avLeafletHit(AVL[name], u, a, c, open, sc, mod, ms ? MS_THICK : 1);
      if (f < 0) continue;
      _valveResult.tissue = TISSUE.VALVE;
      // rheumatic MS: the calcium sits in the fused free edges and the sub-valvular
      // tunnel (f >= MS_CALC_F); the leaflet bodies are fibrotic and pliable, bright
      // but not shadowing (Wilkins GT et al., Br Heart J 1988;60:299: thickening and
      // calcification graded from the leaflet margins toward the base). With the
      // whole leaflet flagged calcific, the apical views had the entire LA in an
      // acoustic shadow.
      _valveResult.echo = ms ? 0.98 : 0.8 + 0.06 * (1 - f);
      _valveResult.calc = ms && f >= MS_CALC_F ? 1 : 0;
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
    if (a < -0.25 || a > cd + 0.05) continue;             // the cusp span (+ the shut cusp's sag)
    let t = a <= 0 ? 0 : a / cd;                          // 0 at hinge, 1 at free edge
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
    // Cusp profile in the (axial a, radius) half-plane, hinge -> free edge.
    //  OPEN: a curtain along the sinus — an S-profile (smoothstep) that leaves the
    //    annulus nearly parallel to the outflow axis, bellies, then curves in to
    //    the free edge. Calcific AS: a straight, rigid cusp of its shut length
    //    swung about the hinge (thick, bright, barely mobile — not a dome).
    //  SHUT (semilunar): each cusp runs from its hinge inward to a central
    //    coaptation just downstream of the hinge line, its belly sagging toward
    //    the ventricle under the diastolic back-pressure (convex toward the
    //    outflow tract, the sinus pocket on the arterial side); the free edges
    //    then meet along the commissural lines (seams below).
    // In between, the two curves blend (sampled, nearest point).
    const asv = name === 'aortic' && !!path.aorticStenosis;
    let dist, bl = 0;
    if (!isCusp) {
      if (a < -0.05 || a > cd) continue;                  // only the leaflet span, hinge→edge
      const sc = t * t * (3 - 2 * t);
      dist = Math.abs(rad - (v.r + (edgeR - v.r) * sc));
    } else {
      const oN = asv ? open / AS_OPEN : open;              // opening relative to this valve's full excursion
      bl = smoothstep(0.02, 0.45, oN);                     // 0 shut .. 1 open shape
      const rho = rad / v.r;
      if (bl === 0 && rho <= 1) {
        // shut: the cusp height is a function of radius (shallow: vertical distance)
        const tb = 1 - rho;
        const ac = SL_COAPT * tb - SL_SAG * Math.sin(Math.PI * tb);
        const sl = (SL_COAPT - SL_SAG * Math.PI * Math.cos(Math.PI * tb)) / v.r;
        dist = Math.abs(a - ac) / Math.sqrt(1 + sl * sl);
        t = tb;
      } else {
        const eR = edgeR / v.r;
        // AS: rigid cusp of the shut length, tilted about the hinge to its orifice
        const aAS = Math.sqrt(Math.max(0, v.r * v.r + SL_COAPT * SL_COAPT - (v.r * (1 - eR)) ** 2));
        let best = 1e9, bi = 0;
        for (let i = 0; i <= SL_N; i++) {
          const ti = i / SL_N;
          const so = asv ? ti : ti * ti * (3 - 2 * ti);
          const ao = (asv ? aAS : cd) * ti, ro = (1 + (eR - 1) * so) * v.r;
          const ac = SL_COAPT * ti - SL_SAG * Math.sin(Math.PI * ti), rc = (1 - ti) * v.r;
          const ai = ac + (ao - ac) * bl, ri = rc + (ro - rc) * bl;
          _slA[i] = ai; _slR[i] = ri;
          const d2 = (a - ai) * (a - ai) + (rad - ri) * (rad - ri);
          if (d2 < best) { best = d2; bi = i; }
        }
        // refine on the two segments either side of the nearest sample
        t = bi / SL_N;
        for (let j = Math.max(0, bi - 1); j < Math.min(SL_N, bi + 1); j++) {
          const ea = _slA[j + 1] - _slA[j], er = _slR[j + 1] - _slR[j];
          const l2 = ea * ea + er * er || 1e-12;
          const u = clamp(((a - _slA[j]) * ea + (rad - _slR[j]) * er) / l2, 0, 1);
          const qa = _slA[j] + ea * u - a, qr = _slR[j] + er * u - rad;
          const d2 = qa * qa + qr * qr;
          if (d2 < best) { best = d2; t = (j + u) / SL_N; }
        }
        dist = Math.sqrt(best);
      }
    }
    // Thin membrane tapering to a fine free edge (thickest at the annular base).
    let th = thick * (0.45 + 0.55 * (1 - t) * (1 - t));
    // Nodulus of Arantius: a fibrous thickening at the CENTRE of each semilunar
    // cusp's free edge, where the three cusps meet — bulk up the free-edge tip at
    // the cusp centre so the coaptation shows the three little nodules. (Aortic:
    // ~1 mm, so the shut cusps read on TEE as a fine line, not a bar with a bulb.)
    const ao = name === 'aortic';
    let nod = 0;
    if (isCusp && t > 0.7) nod = (ao ? AO_NODULE : 0.11) * ((t - 0.7) / 0.3) * (1 - cptr);
    th += nod;
    // calcific AS: fibrous cusps 3-5 mm thick carrying discrete calcium nodules
    // (~3-4 mm masses on a 3-D lattice fixed to the valve, on the cusp bodies and
    // bases more than the free edges) — near-white with a distal acoustic shadow
    // (calc 2), the rest of the cusp fibrous (calc 0). Open, the cusps part along
    // their commissural lines from mid-cusp to the free edge (partial fusion at the
    // bases only): three thick cusp masses round a small stellate orifice.
    let nodule = false;
    if (asv) {
      const lump = (0.5 + 0.5 * Math.sin(dx * 4.3 + 1.3)) * (0.5 + 0.5 * Math.sin(dy * 3.8 + 0.7)) * (0.5 + 0.5 * Math.sin(dz * 3.5 + 2.1));
      nodule = lump > AS_NODULE && t < 0.85;
      th = th * 2.2 + 0.045 + (nodule ? 0.12 * (1 - 0.5 * t) : 0);
      if (bl > 0 && t > AS_SLIT_T && (Math.abs(kh) / 3) * rad < AS_SLIT) continue;
    }
    let hit = dist <= th;
    // Coaptation zone of the shut semilunar valve: neighbouring cusps' free
    // edges press together along the three commissural lines, from the cusp
    // surface up to the free margin, which rises from the central nodules toward
    // the commissures at the sinus wall — the short-axis "Mercedes" Y, and the
    // thin central closure line in long axis.
    if (!hit && isCusp && open < 0.25 * (asv ? AS_OPEN : 1)) {
      const angDist = Math.abs(kh) / 3;                   // angular distance to a commissure
      const rho = rad / v.r;
      if (rho < 0.92 && angDist * rad < (asv ? 0.126 : ao ? AO_COAPT_W : 0.04)) {
        const tb = 1 - rho;
        const base = SL_COAPT * tb - SL_SAG * Math.sin(Math.PI * tb);
        const top = SL_COAPT + SL_COAPT_H + rho * (cd - SL_COAPT - SL_COAPT_H);
        if (a >= base && a <= top) { hit = true; t = 1; }
      }
    }
    if (!hit) continue;
    let echo = 0.8 + 0.05 * (1 - t) + (nod > 0.03 ? 0.04 : 0); // nodule reads a touch brighter
    const calc = asv ? nodule : (name === 'mitral' && path.mitralStenosis);
    if (calc) echo = 0.98; // calcified / restricted → bright, thickened
    else if (asv) echo = 1.0; // the sclerotic, thickened AS cusp between the nodules (renderer: a few dB above a normal cusp)
    _valveResult.tissue = TISSUE.VALVE;
    _valveResult.echo = echo;
    _valveResult.calc = asv ? (nodule ? 2 : 0) : calc ? 1 : 0;
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
// left atrium; negative = retrograde, back up the vein. It is the solved
// pulmonary-venous inflow of the circulation (H.Qpv, mL/s) over the vein-mouth
// area, so the S wave (atrial relaxation and annular descent), the D wave (mitral
// opening), any reversal (severe MR in systole; an atrial kick that lifts the LA
// above the venous pressure) and their ratios follow the case's LA pressure.
// PV_AREA (cm2, the effective total area the sampled vein carries) puts the
// normal systolic peak at ~0.5 m/s. In the normal model the S and D peaks are
// about equal (S/D ~0.95), as in young adults with an E-dominant mitral inflow
// (E/A ~1.7); a raised LAP with restrictive filling (DCM) makes it D-dominant.
const PV_AREA = 2.9;
function pvFlowSpeed(H) {
  return H && H.Qpv ? H.Qpv / (PV_AREA * 100) : 0;
}

// Left-atrial-appendage flow (m/s, + = out of the appendage into the LA body):
// emptying with atrial contraction (normal peak ~0.5-0.6 m/s, > 0.4 in sinus
// rhythm), filling back into it in early systole as the LA relaxes (~0.4-0.5 m/s),
// and a small early-diastolic emptying wave as the mitral valve opens (Agmon Y et
// al., JACC 1999;34:1867-77; ASE/EACVI cardiac sources of embolism).
function laaFlowSpeed(phase) {
  const a = pulse(0.86, 1.0, phase), f = pulse(0.02, 0.20, phase), e = pulse(0.50, 0.68, phase);
  return 0.66 * a - 0.46 * f + 0.2 * e;
}

// Hepatic-vein flow speed (m/s, + toward the RA): the solved caval inflow of the
// circulation (H.Qvr, mL/s, read at the cycle phase like every other solved
// signal: hemoFrame keeps no separate right-heart clock) over an effective vessel
// area. The S wave is the RA's fall in pressure as the tricuspid annulus descends
// (x descent), the D wave its emptying once the tricuspid opens (y descent), and
// a reversal wherever the RA pressure rises above the systemic venous pressure:
// the atrial kick (AR, pushed back into the soft central veins of a normal RAP:
// hemodynamics.js Csv), the regurgitant volume of severe TR (systolic reversal),
// a stiff RA's tall a wave (pulmonary hypertension), and in tamponade the loss of
// the y descent (blunted D) with expiratory diastolic reversal, each respiratory
// beat read from its own solved circulation. HV_AREA (cm2) puts the normal S
// peak at ~0.37 m/s. The veins are distensible: with a raised RA pressure the
// IVC is plethoric and the hepatic veins dilated (ASE 2010 RA-pressure
// estimation), so the area scales with the square of the IVC radius the anatomy
// draws for the case's RAP (ivcOfRap: 0.9 cm at RAP <= 4.6, 1.25 cm at >= 9):
// the same backward volume of a tall RA wave reads ~half as fast in the
// plethoric veins of TR, PH, DCM or tamponade.
const HV_AREA = 4.0;
function hvFlowSpeed(H) {
  if (!H || !H.Qvr) return 0;
  const r = ivcOfRap(H.RAP || 0).r / 0.9;
  return H.Qvr / (HV_AREA * r * r * 100);
}
const _hvCache = new WeakMap();

// Off-axis Gaussian falloff: 1 on the centreline, ->0 as perpendicular distance
// d grows past the stream half-width w.
// Inflow-stream profile (0..1) across (perp, half-width w) and along the
// stream: flat over the central half and the first 2 cm below the annulus,
// Gaussian shoulders, then decaying toward the apex.
function inflowShape(perp, w, along) {
  const e = perp - 0.5 * w, h = 0.5 * w;
  const across = e <= 0 ? 1 : Math.exp(-(e / h) * (e / h));
  return across * Math.exp(-clamp(along - 2, 0, 4) / 3.5);
}
function gauss(d, w) {
  const t = d / w;
  return Math.exp(-t * t);
}

// Perpendicular distance from point (px,py,pz) to the line through anchor `a`
// with UNIT direction `u`. The signed distance travelled along the axis is
// stashed in _axis.along (module scope) to avoid allocating a return tuple.
const _axis = { along: 0 };
// The upper infundibulum's axis (A.rvot[1] -> A.rvot[2], just below the pulmonary
// valve) of the anatomy frame last asked for: its unit direction in _infU, its
// length returned. (velocityAt asks once per sample; the frame changes per phase.)
const _infU = [0, 0, 0];
let _infA = null, _infLen = 1;
function infundibulumAxis(A) {
  if (A !== _infA) {
    const a = A.rvot[1], b = A.rvot[2];
    const ux = b[0] - a[0], uy = b[1] - a[1], uz = b[2] - a[2];
    _infLen = Math.hypot(ux, uy, uz) || 1;
    _infU[0] = ux / _infLen; _infU[1] = uy / _infLen; _infU[2] = uz / _infLen;
    _infA = A;
  }
  return _infLen;
}
// Is a point at the edge of the RVOT -> PA stream in that lumen? Below the
// pulmonary valve (along < 0, cm along the trunk axis from the valve) the RV's
// outflow tract, above it the PA trunk — but not where the left PA runs over the
// left atrium (beyond ~2 cm, near the bifurcation). Exact lumen SDFs, tried in
// the order the point is most likely to be in.
function inRvotLumen(px, py, pz, A, along) {
  if (along < 0) return lumenDist(px, py, pz, A, 'RV') < 0 || lumenDist(px, py, pz, A, 'PA') < 0;
  if (lumenDist(px, py, pz, A, 'PA') < 0) return along < 2 || lumenDist(px, py, pz, A, 'LA') > 0;
  return lumenDist(px, py, pz, A, 'RV') < 0;
}
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
// RV outflow: just below the pulmonary valve the flow runs along the valve /
// trunk axis U_PA (review 12, TTE-5/PHY-11: the infundibulum-to-valve chord,
// PV - RVOT[1], lay ~32 deg off the trunk, so the PSAX-AV beam that is aligned with
// the drawn main PA met it at ~24 deg and no default RVOT trace was recorded)
const AX_RVOT = LM.U_PA;
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
// tricuspid -> into the RA: a central TR jet runs back up the RV inflow axis,
// bent a little toward the RA body
const AX_TR = (() => {
  const toRA = unit(_sub3(RA_MID, LM.T));
  return unit([0.35 * toRA[0] - 0.65 * AX_TRICUSPID_IN[0], 0.35 * toRA[1] - 0.65 * AX_TRICUSPID_IN[1], 0.35 * toRA[2] - 0.65 * AX_TRICUSPID_IN[2]]);
})();
// (perpendicular to the LV long axis, as the septum runs along it: the jet
// then lies in the short-axis plane as well as in the four-chamber plane)
const AX_VSD = unit([-1, -0.175, 0.15]);        // LV -> RV across the septum
// LA -> RA across the septum, the stream bending toward the tricuspid valve as
// it enters the RA (so the subcostal and ME views align with it)
const AX_ASD = unit([LM.IAS_N[0] + (LM.T[0] - LM.IAS.asdC[0]) / Math.hypot(...LM.T.map((v, i) => v - LM.IAS.asdC[i])),
  LM.IAS_N[1] + (LM.T[1] - LM.IAS.asdC[1]) / Math.hypot(...LM.T.map((v, i) => v - LM.IAS.asdC[i])),
  LM.IAS_N[2] + (LM.T[2] - LM.IAS.asdC[2]) / Math.hypot(...LM.T.map((v, i) => v - LM.IAS.asdC[i]))]);
// mid muscular septum, where the mid-papillary PSAX and the A4C planes cross, so the
// defect is seen in both standard views
const VSD_CORE = [-2.4, -3.6, LM.A4C_Z(-2.4, -3.6)];
const ASD_CORE = LM.IAS.asdC;        // the secundum ASD in the fossa ovalis
// Read-only flow axes and jet sites, for the spectral-Doppler gate placement
// (echo.js) so the gates sit on the same streams the velocity field draws.
export const FLOW_AXES = { mitralIn: AX_MITRAL_IN, tricuspidIn: AX_TRICUSPID_IN, lvot: AX_LVOT, rvot: AX_RVOT,
  mr: AX_MR, mrEcc: AX_MR_ECC, tr: AX_TR, vsd: AX_VSD, asd: AX_ASD };
export const FLOW_SITES = { vsd: VSD_CORE, asd: ASD_CORE, laMid: LA_MID, raMid: RA_MID };

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
// `variance` (0..1) is the spread of velocities within a resolution cell (the
// colour map's green axis and the spectral broadening); `turbulent` is kept as
// its thresholded form (variance >= 0.5) for older callers.
const _velResult = { vx: 0, vy: 0, vz: 0, speed: 0, flow: FLOW.NONE, turbulent: false, variance: 0 };
let _velHasBest = false;

// Lifted out of velocityAt so no closure is allocated per sample; state (the
// current best) lives in the module-scoped singleton above. Keeps the fastest
// contribution. dir components are expected pre-normalised (unit axes).
// `turb` is a variance 0..1, or a boolean (true = a fully turbulent cell).
function considerFlow(inside, dirx, diry, dirz, speed, flow, turb) {
  if (!inside || speed <= 0.01) return;
  if (_velHasBest && speed <= _velResult.speed) return;
  const vr = typeof turb === 'number' ? (turb < 0 ? 0 : turb > 1 ? 1 : turb) : (turb ? 0.85 : 0);
  _velResult.vx = dirx * speed;
  _velResult.vy = diry * speed;
  _velResult.vz = dirz * speed;
  _velResult.speed = speed;
  _velResult.flow = flow;
  _velResult.variance = vr;
  _velResult.turbulent = vr >= 0.5;
  _velHasBest = true;
}

// Shared regurgitant-jet contribution (MR / TR are geometrically identical). A
// regurgitant jet leaves the orifice as a NARROW vena contracta and widens
// downstream as a turbulent free jet (a flame of half-angle ~15 deg), its
// centreline speed decaying with travel; upstream, on the ventricular side, flow
// converges on the orifice through hemispheric isovelocity shells (PISA, v ~ 1/r^2).
// `O` is the (moving) orifice centre, `vInst` the instantaneous jet speed from
// the pressure drop across the valve, `len` the jet length (cm, from the
// regurgitant severity). With `wall` (an eccentric prolapse jet) the jet
// attaches to the atrial wall it is aimed at (Coanda effect) and runs along it
// as a thin sheet, instead of spreading into the middle of the atrium.
let _jetPhase = 0;                 // cycle phase for the jet's slow meander (set by velocityAt)
const _jetBasis = new WeakMap();   // per jet axis: two unit vectors across it
// a free jet keeps its vena-contracta speed along a short potential core
// (a few orifice widths) before it starts to decay
const JET_CORE = 0.4;
function regurgJet(px, py, pz, O, axis, downstreamGate, upstreamInside, p, vInst, len, wall, vc) {
  if (!(vInst > 0.05)) return;
  const dx = px - O[0], dy = py - O[1], dz = pz - O[2];
  const along = dx * axis[0] + dy * axis[1] + dz * axis[2];
  // turbulence potential: a jet well above ~1 m/s breaks up into a mosaic
  const vt = clamp((vInst - 1.0) / 1.5, 0, 1);
  if (downstreamGate && along >= 0) {
    const lenFree = wall ? Math.min(len, p.wallAttach) : len;
    // A free jet is a plume, not a beam-aligned cone: the vena contracta, then a
    // flame widening at ~0.3 x the distance travelled, its centreline meandering
    // slowly (shed eddies, drifting through the cycle), its edge scalloped, and
    // a rounded, fading head instead of a straight-edged wedge.
    let rx = dx - along * axis[0], ry = dy - along * axis[1], rz = dz - along * axis[2];
    let B = _jetBasis.get(axis);
    if (!B) { const e1 = Math.abs(axis[1]) < 0.9 ? _cross(axis, [0, 1, 0]) : _cross(axis, [1, 0, 0]); B = [e1, _cross(axis, e1)]; _jetBasis.set(axis, B); }
    const e1 = B[0], e2 = B[1];
    const ph = 6.2832 * _jetPhase;
    const m1 = 0.13 * along * Math.sin(1.9 * along - 2 * ph), m2 = 0.09 * along * Math.sin(1.3 * along + 3 * ph + 1.1);
    rx -= m1 * e1[0] + m2 * e2[0]; ry -= m1 * e1[1] + m2 * e2[1]; rz -= m1 * e1[2] + m2 * e2[2];
    const perp = Math.hypot(rx, ry, rz);
    let w = ((vc || p.vc) + along * p.spread) * (1 + 0.12 * Math.sin(3.3 * along - 4 * ph));
    const tip = along / (lenFree + 0.6);
    if (tip > 0.6) w *= Math.sqrt(Math.max(0, 1 - ((tip - 0.6) / 0.4) * ((tip - 0.6) / 0.4)));   // rounded head
    const q = perp / Math.max(w, 1e-3);
    if (q < 1.25 && along < lenFree + 1) {
      const prof = q < 0.8 ? 1 : (1.25 - q) / 0.45;   // flat turbulent core, short soft edge
      const fade = along < lenFree ? 1 : 1 - (along - lenFree);
      // variance: low at the (laminar) vena contracta, a mosaic core building
      // downstream, a laminar shell at the plume's edge
      const vr = vt * (0.15 + 0.85 * smoothstep(0.1, 0.8, along)) * (1 - smoothstep(0.55, 1.05, q)) * (0.4 + 0.6 * fade);
      considerFlow(true, axis[0], axis[1], axis[2], vInst * prof * fade / (1 + Math.max(0, along - JET_CORE) / p.L), p.flow, vr);
    }
    if (wall && along > 0.3) {
      const sd = Math.hypot(dx, dy, dz);
      if (along / sd > 0.55 && sd < len + 1) {
        const d = -lumenDist(px, py, pz, wall, 'LA_BODY');   // depth inside the LA lumen
        if (d > -0.1 && d < 0.9) {
          const t = (d - 0.35) / 0.5, g = 1 - t * t;
          if (g > 0) {
            const fade = sd < len ? 1 : 1 - (sd - len);
            const inv = 1 / sd;
            considerFlow(true, dx * inv, dy * inv, dz * inv, vInst * g * fade / (1 + sd / (1.6 * p.L)), p.flow, vt * (0.3 + 0.5 * g));
          }
        }
      }
    }
  }
  // proximal isovelocity surface area (PISA): radial convergence on the orifice
  // (laminar: clean aliasing shells, no variance)
  if (along < 0 && upstreamInside) {
    const dOrif = Math.hypot(dx, dy, dz);
    if (dOrif < p.pisaR) {
      const rr = p.pisa / (dOrif + 0.05);
      // (capped below the jet speed: the shells nearest the orifice are tiny,
      // so they never outshine the vena contracta on CW)
      const conv = vInst * Math.min(0.45, rr * rr);
      const inv = 1 / (dOrif + 1e-4);
      considerFlow(true, -dx * inv, -dy * inv, -dz * inv, conv, p.flow, conv > 1.2 ? 0.2 : 0);
    }
  }
}

// Per-lesion tuning for regurgJet: vena-contracta half-width vc (cm), flame
// spread (tan of the half-angle), centreline decay length L (cm), PISA radius
// scale and reach, and where an eccentric jet meets the wall.
const MR_JET_PARAMS = { vc: 0.24, spread: 0.27, L: 2.2, pisa: 0.3, pisaR: 1.6, wallAttach: 1.4, flow: FLOW.MR_JET };
const TR_JET_PARAMS = { vc: 0.24, spread: 0.27, L: 2.0, pisa: 0.28, pisaR: 1.5, wallAttach: 1.4, flow: FLOW.TR_JET };
const JET_LEN = { mild: 1.6, moderate: 3.0, severe: 4.6 };
// restrictive VSD: a ~0.5 cm vena contracta, a flame of ~12 deg half-angle,
// reaching ~3.5 cm from the defect (2-4 cm into the RV)
const VSD_JET_PARAMS = { vc: 0.25, spread: 0.22, L: 1.8, pisa: 0.3, pisaR: 1.5, wallAttach: 1.4, flow: FLOW.VSD_JET };
const VSD_JET_LEN = 2.8;
const _negAxes = new WeakMap();
function negAxis(a) { let n = _negAxes.get(a); if (!n) { n = [-a[0], -a[1], -a[2]]; _negAxes.set(a, n); } return n; }

// Radius (cm) of the modelled VSD channel. A ~1 cm defect with a ~4.5-5 m/s jet is
// a RESTRICTIVE defect: the circulation carries a modest L-R shunt (Qp:Qs ~1.4).
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
  if (flow === FLOW.AS_JET && path.aorticStenosis) { _cw.peak = H.vAoPeak; _cw.env = H.ejAS; }
  else if (flow === FLOW.MITRAL_IN && path.mitralStenosis) { _cw.peak = H.vMitPeak; _cw.env = H.ejMv; }
  // a regurgitant jet is read the same way: CW through it gives the modelled
  // ventricle-atrium peak (MR; TR Vmax, behind the PASP estimate) along its
  // envelope, whatever the frame
  else if (flow === FLOW.MR_JET && H.vMRPeak > 0) { _cw.peak = H.vMRPeak; _cw.env = H.ejMr; }
  else if (flow === FLOW.TR_JET && H.vTRPeak > 0) { _cw.peak = H.vTRPeak; _cw.env = H.ejTr; }
  return _cw;
}

// Secundum ASD shunt: a broad, laminar left-to-right column. On the LA side the
// flow converges on the defect only within ~1 cm of it; through the defect it
// runs across the septum (IAS_N) with a blunt profile as wide as the hole; on
// the RA side it streams on as a column that bends toward the tricuspid valve
// over 3-4 cm and slows as it goes. Only the core of the column (about 1 m/s)
// aliases at a standard colour scale. Variance stays low (laminar).
const _asdCache = new WeakMap();
function asdPath(G) {
  let c = _asdCache.get(G.A);
  if (c) return c;
  const C = ASD_CORE, N = LM.IAS_N, T = G.A.valves.tricuspid.c;
  const K = [C[0] + N[0] * 2.6, C[1] + N[1] * 2.6, C[2] + N[2] * 2.6];
  const pts = [];
  for (let i = 0; i <= 12; i++) {
    const t = i / 12, a = (1 - t) * (1 - t), b = 2 * (1 - t) * t, d = t * t;
    pts.push([a * C[0] + b * K[0] + d * T[0], a * C[1] + b * K[1] + d * T[1], a * C[2] + b * K[2] + d * T[2]]);
  }
  const segs = [];
  let L = 0;
  for (let i = 0; i < 12; i++) {
    const a = pts[i], e = pts[i + 1], u = [e[0] - a[0], e[1] - a[1], e[2] - a[2]], l = Math.hypot(u[0], u[1], u[2]);
    segs.push({ a, u: [u[0] / l, u[1] / l, u[2] / l], l, s0: L });
    L += l;
  }
  c = { segs, L };
  _asdCache.set(G.A, c);
  return c;
}
function asdFlow(px, py, pz, G, H) {
  const v = H.vAsd || 0;
  if (Math.abs(v) < 0.02) return;
  const C = ASD_CORE, N = LM.IAS_N, U = LM.IAS.asdU, Bk = LM.IAS.asdB;
  const Ru = LM.IAS.asdR, Rb = LM.IAS.asdRb;
  const dx = px - C[0], dy = py - C[1], dz = pz - C[2];
  if (dx * dx + dy * dy + dz * dz > 30) return;
  const s = iasSide(px, py, pz);
  const du = dx * U[0] + dy * U[1] + dz * U[2], db = dx * Bk[0] + dy * Bk[1] + dz * Bk[2];
  const sg = v >= 0 ? 1 : -1, av = Math.abs(v);
  // the upstream side is the LA for left-to-right flow (the RA during the brief reversal)
  const up = s * sg;
  if (up <= 0.15) {
    // within and just upstream of the defect: across the septum, blunt profile
    const r2 = (du / Ru) * (du / Ru) + (db / Rb) * (db / Rb);
    if (up > -0.5 && r2 < 1.6) {
      const conv = up > -0.1 ? 1 : Math.exp((up + 0.1) / 0.15);       // converging only within the last few mm
      const prof = r2 < 1 ? 1 - 0.75 * r2 : 0.25 * (1.6 - r2) / 0.6;
      considerFlow(true, N[0] * sg, N[1] * sg, N[2] * sg, av * prof * conv, FLOW.ASD_JET, 0.05);
    }
    if (up > -0.1) return;
    return;
  }
  if (sg < 0) return;                                   // the brief reversal stays at the defect
  // RA side: nearest point on the column's centreline
  const P = asdPath(G);
  let best = 1e9, bs = 0, bu = null, bx = 0, by = 0, bz = 0;
  for (const g of P.segs) {
    const qx = px - g.a[0], qy = py - g.a[1], qz = pz - g.a[2];
    const t = clamp(qx * g.u[0] + qy * g.u[1] + qz * g.u[2], 0, g.l);
    const rx = qx - g.u[0] * t, ry = qy - g.u[1] * t, rz = qz - g.u[2] * t;
    const d2 = rx * rx + ry * ry + rz * rz;
    if (d2 < best) { best = d2; bs = g.s0 + t; bu = g.u; bx = rx; by = ry; bz = rz; }
  }
  const f = bs / P.L;                                    // 0 at the defect -> 1 at the valve
  if (f > 0.85) return;
  // the column fans out as it goes (more across the long axis of the oval
  // defect), slowing as a free jet does (~1/width): only its core near the
  // defect is above a standard colour scale
  const gu = 1 + f, gb = 1 + 1.5 * f;
  const ru = (bx * U[0] + by * U[1] + bz * U[2]) / (Ru * 1.2 * gu), rb = (bx * Bk[0] + by * Bk[1] + bz * Bk[2]) / (Rb * 0.8 * gb);
  const r2 = ru * ru + rb * rb;
  if (r2 >= 1.3) return;
  const prof = r2 < 1 ? 1 - 0.75 * r2 : 0.25 * (1.3 - r2) / 0.3;
  const sp = av * prof / Math.pow(gu * gb, 0.6) * (f > 0.7 ? (0.85 - f) / 0.15 : 1);
  considerFlow(true, bu[0], bu[1], bu[2], sp, FLOW.ASD_JET, f < 0.2 && r2 < 0.4 ? 0.25 : 0.08);
}

// Ascending aorta, arch and descending thoracic aorta: a blunt (plug) systolic
// profile along the vessel, the stroke volume's forward flow (~1 m/s), with a
// small early-diastolic reversal after the aortic valve closes. Read by the
// suprasternal and TEE descending-aorta views (red up the ascending limb, blue
// down the descending one from the notch).
const _aoCache = new WeakMap();
function aortaPath(A) {
  let c = _aoCache.get(A);
  if (c) return c;
  const segs = [];
  const add = (a, b, r1, r2) => {
    const u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]], l = Math.hypot(u[0], u[1], u[2]);
    if (l > 1e-6) segs.push({ a, u: [u[0] / l, u[1] / l, u[2] / l], l, r1, r2 });
  };
  const ao = A.ao;
  for (const g of ao.ascSegs || []) add(g.a, g.b, g.r1, g.r2);
  const arch = ao.arch || [];
  for (let i = 0; i < arch.length; i++) {
    if (i > 0 && Math.hypot(arch[i][0][0] - arch[i - 1][1][0], arch[i][0][1] - arch[i - 1][1][1], arch[i][0][2] - arch[i - 1][1][2]) > 1e-3) break;   // head/neck branches follow
    add(arch[i][0], arch[i][1], arch[i][2], arch[i][2]);
  }
  if (ao.dta) { add(ao.dta.top, ao.dta.a, ao.dta.r, ao.dta.r); add(ao.dta.a, ao.dta.b, ao.dta.r, ao.dta.r); }
  let lo = [1e9, 1e9, 1e9], hi = [-1e9, -1e9, -1e9];
  for (const g of segs) for (const q of [g.a, [g.a[0] + g.u[0] * g.l, g.a[1] + g.u[1] * g.l, g.a[2] + g.u[2] * g.l]]) for (let k = 0; k < 3; k++) { lo[k] = Math.min(lo[k], q[k] - 2); hi[k] = Math.max(hi[k], q[k] + 2); }
  c = { segs, lo, hi };
  _aoCache.set(A, c);
  return c;
}
function aorticFlow(px, py, pz, G, H, phase) {
  if (!G.A || !G.A.ao) return;
  let v = 0;
  if (H.ejAo > 0.001) v = clamp(1.05 * (H.forwardSV || 70) / 70, 0.6, 1.3) * H.ejAo;
  else v = -0.14 * pulse(0.39, 0.56, phase);            // early-diastolic reversal
  if (Math.abs(v) < 0.02) return;
  const asV = G.path && G.path.aorticStenosis ? G.A.valves.aortic.c : null;
  const P = aortaPath(G.A);
  if (px < P.lo[0] || py < P.lo[1] || pz < P.lo[2] || px > P.hi[0] || py > P.hi[1] || pz > P.hi[2]) return;
  for (const g of P.segs) {
    const qx = px - g.a[0], qy = py - g.a[1], qz = pz - g.a[2];
    const t = qx * g.u[0] + qy * g.u[1] + qz * g.u[2];
    if (t < -0.2 || t > g.l + 0.2) continue;
    const rx = qx - g.u[0] * t, ry = qy - g.u[1] * t, rz = qz - g.u[2] * t;
    const r = g.r1 + (g.r2 - g.r1) * clamp(t / g.l, 0, 1);
    const q = Math.hypot(rx, ry, rz) / r;
    if (q >= 1) continue;
    const sg = v >= 0 ? 1 : -1;
    if (asV && v > 0) {
      // severe AS: beyond the valve the lumen holds ONE jet and, around it, slow
      // disturbed recirculating flow, not a laminar plug at 0.6-0.8 m/s. The jet
      // leaves the root along the LVOT axis (velocityAt) and is then turned by
      // the aortic wall to run up the ascending aorta: here it follows the
      // vessel, widening ~0.25 cm per cm from the orifice radius and slowing.
      // Around it, within 4 cm of the valve, the plug is replaced by eddies of
      // <= 0.3 m/s (forward and back), blending into the reattached stream by 5.5 cm.
      const dV = Math.hypot(px - asV[0], py - asV[1], pz - asV[2]);
      if (dV < 5.5) {
        const rr = q * r, w = Math.sqrt(H.AVA / Math.PI) + 0.25 * dV, qj = rr / w;
        if (H.vAS > 0.05 && qj < 1.25 && dV < 5) {
          const prof = (qj < 0.8 ? 1 : (1.25 - qj) / 0.45) * (dV > 4 ? 5 - dV : 1);
          considerFlow(true, g.u[0], g.u[1], g.u[2], H.vAS * prof / (1 + Math.max(0, dV - 0.8) / 3.5), FLOW.AS_JET,
            clamp((H.vAS - 1.0) / 1.0, 0, 1) * 0.9);
        }
        const k = smoothstep(4, 5.5, dV), fw = v * (1 - q * q * q) * k;
        const e = Math.sin(2.3 * px + 1.7 * pz + 12.57 * phase) * Math.sin(1.9 * py - 2.6 * pz + 6.28 * phase);
        const sw = Math.min(0.3, 0.3 * v) * Math.abs(e) * (1 - k);
        if (fw >= sw) considerFlow(true, g.u[0], g.u[1], g.u[2], fw, FLOW.AORTA, 0.3 * (1 - k));
        else { const se = e >= 0 ? 1 : -1; considerFlow(true, g.u[0] * se, g.u[1] * se, g.u[2] * se, sw, FLOW.AORTA, 0.3); }
        continue;
      }
    }
    considerFlow(true, g.u[0] * sg, g.u[1] * sg, g.u[2] * sg, Math.abs(v) * (1 - q * q * q), FLOW.AORTA, 0);
  }
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
  _jetPhase = phase;

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
    // anchored on the MOVING annulus, so the inflow starts at the valve in every frame
    const cM = G.A.valves.mitral.c;
    const mPerp = offAxisDist(px, py, pz, cM, AX_MITRAL_IN);
    const mAlong = _axis.along;
    // (MS: the stream is as wide as the stenotic orifice, sqrt(MVA / normal area):
    // a ~1 cm candle-flame jet for a 1 cm2 valve, not the normal 2.5 cm inflow)
    const msW = path.mitralStenosis ? Math.sqrt((MS_AREA[gradeOf(path)] || MS_AREA.severe) / MVA_NORM) : 1;
    const mW = 1.25 * msW * (1 - 0.5 * clamp(mAlong / 4.5, 0, 1)); // narrows toward apex
    // FLAT core (review 12, PHY-7): over the first 2 cm and the central half of
    // the stream the speed is the solved Q/A, so a PW gate at the leaflet tips
    // reads vMitPeak x ejMv in every loading state (the small inspiratory LV of
    // tamponade no longer slides a sloping profile past the gate), and the trace
    // E equals the reported E
    const mShape = inflowShape(mPerp, mW, mAlong);
    const mSpeed = H.vMitPeak * H.ejMv * mShape;
    // the LV proxy ellipsoid stops short of the annulus, so the basal inflow
    // tract just below the valve (within the inflow cone) counts as LV too
    const inMitTract = inLV || (mAlong < 1.6 && mPerp < 1.3 && !inAo);
    if (inMitTract && mAlong > -0.05 && mSpeed > 0.01) {
      // superimpose the divergence-free filling vortex: a pure swirl about
      // VORTEX_AXIS (tangential = axis × r_perp) whose magnitude rises then
      // decays with radius (Lamb-Oseen-like), so the inflow curls behind the
      // anterior leaflet instead of reading as a straight jet.
      const rx = px - VORTEX_CORE[0] - (cM[0] - VALVES.mitral.c[0]);
      const ry = py - VORTEX_CORE[1] - (cM[1] - VALVES.mitral.c[1]);
      const rz = pz - VORTEX_CORE[2] - (cM[2] - VALVES.mitral.c[2]);
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
      // (MS: a turbulent mosaic jet below the orifice)
      considerFlow(true, vx / sp, vy / sp, vz / sp, sp, FLOW.MITRAL_IN,
        path.mitralStenosis && mSpeed > 1.2 && mAlong > 0.15);
    }
    // atrial-side flow convergence: blood in the LA accelerates toward the
    // orifice (speed ~ 1/(distance + 0.6 cm), equal to the annular speed at the
    // orifice), so the inflow colour starts in the atrium and narrows into the
    // valve instead of being cut flat at the annulus. In MS this is the PISA dome.
    const vAnnM = H.vMitPeak * H.ejMv;
    if (mAlong < 0 && vAnnM > 0.02 && !inRV && !inRA && !inAo && iasSide(px, py, pz) < 0) {
      const dO = Math.hypot(px - cM[0], py - cM[1], pz - cM[2]);
      if (dO < 2.2 && dO > 1e-3) {
        const k = 0.6 / dO;
        const dx = (cM[0] - px) * k + 0.4 * AX_MITRAL_IN[0], dy = (cM[1] - py) * k + 0.4 * AX_MITRAL_IN[1], dz = (cM[2] - pz) * k + 0.4 * AX_MITRAL_IN[2];
        const dl = Math.hypot(dx, dy, dz) || 1;
        if (path.mitralStenosis) {
          // MS: hemispheric proximal isovelocity shells converging on the stenotic
          // orifice, v = v_orifice (r0 / r)^2 with r0 the orifice radius (capped at
          // the orifice speed): ONE laminar red-to-blue aliasing shell ~1 cm from
          // the orifice, the innermost few mm at jet speed (mosaic)
          const r0 = Math.sqrt((MS_AREA[gradeOf(path)] || MS_AREA.severe) / Math.PI);
          const sp = vAnnM * Math.min(1, (r0 / dO) * (r0 / dO)) * (dO > 1.8 ? (2.2 - dO) / 0.4 : 1);
          const inv = 1 / dO;
          considerFlow(true, (cM[0] - px) * inv, (cM[1] - py) * inv, (cM[2] - pz) * inv, sp, FLOW.MITRAL_IN, sp > 1.5);
        } else {
        const sp = vAnnM * 0.6 / (dO + 0.6) * (dO > 1.8 ? (2.2 - dO) / 0.4 : 1);
        considerFlow(true, dx / dl, dy / dl, dz / dl, sp, FLOW.MITRAL_IN, false);
        }
      }
    }

    // tricuspid inflow -> RV apex (right heart not lumped-modelled; the filling
    // envelope on the right heart's clock, at a lower normal velocity)
    const cT = G.A.valves.tricuspid.c;
    const tPerp = offAxisDist(px, py, pz, cT, AX_TRICUSPID_IN);
    const tAlong = _axis.along;
    const tW = 1.35 * (1 - 0.5 * clamp(tAlong / 4.5, 0, 1));
    const tShape = inflowShape(tPerp, tW, tAlong);   // flat core, as the mitral stream
    const tSpeed = (H.vTv != null ? H.vTv : 0.55 * H.ejTv) * tShape;
    considerFlow((inRV || (tAlong < 1.6 && tPerp < 1.3 && !inAo)) && tAlong > -0.05, AX_TRICUSPID_IN[0], AX_TRICUSPID_IN[1], AX_TRICUSPID_IN[2],
      tSpeed, FLOW.TRICUSPID_IN, false);
    const vAnnT = H.vTv != null ? H.vTv : 0.55 * H.ejTv;   // RA-side convergence, as for the mitral
    if (tAlong < 0 && vAnnT > 0.02 && !inLV && !inLA && !inAo && iasSide(px, py, pz) > 0) {
      const dO = Math.hypot(px - cT[0], py - cT[1], pz - cT[2]);
      if (dO < 2.2 && dO > 1e-3) {
        const k = 0.6 / dO;
        const dx = (cT[0] - px) * k + 0.4 * AX_TRICUSPID_IN[0], dy = (cT[1] - py) * k + 0.4 * AX_TRICUSPID_IN[1], dz = (cT[2] - pz) * k + 0.4 * AX_TRICUSPID_IN[2];
        const dl = Math.hypot(dx, dy, dz) || 1;
        const sp = vAnnT * 0.6 / (dO + 0.6) * (dO > 1.8 ? (2.2 - dO) / 0.4 : 1);
        considerFlow(true, dx / dl, dy / dl, dz / dl, sp, FLOW.TRICUSPID_IN, false);
      }
    }
  }

  // --- systolic ejection: a THIN column hugging the outflow centreline only,
  //     NOT a whole-ventricle fill. Present only during ejection (0.04-0.34).
  if (H.ejAo > 0.001) {
    // LVOT -> aorta (centreline anchored at the aortic valve). Clip the stream
    // to the aortic (downstream, +along) side of the valve plane so it reads as
    // ONE directed jet aimed into the aorta rather than two symmetric lobes
    // flanking the annulus (the off-axis Gaussian alone is symmetric in `along`).
    const lPerp = offAxisDist(px, py, pz, G.A.valves.aortic.c, AX_LVOT);   // the moving valve
    const lAlong = _axis.along;
    // LVOT/aortic-valve jet peak from continuity (v = SV / (AVA·ET)): with the
    // SAME stroke volume forced through a reduced AVA, AS auto-scales to ~4 m/s
    // and normal sits ~1 m/s — severity now follows physiology, not a constant.
    // (AS: the jet's speed in time, vAS(t), times its spatial profile. On the
    // ventricular side of the orifice the flow converges on it (proximal
    // acceleration, v ~ v_orifice (r0 / (r0 + d))^2 with r0 the radius of the
    // ~0.6 cm2 orifice), so the jet speed is reached only at the valve and an
    // LVOT PW gate the annulus swings toward late in systole reads the LVOT flow,
    // not the jet)
    const AS_R0 = 0.45;
    if (path.aorticStenosis) {
      // ONE turbulent jet (review 12, PAT-8 / TEE-12): a flat vena-contracta core
      // of the orifice radius sqrt(AVA / pi), widening ~0.25 cm per cm into the
      // ascending aorta behind a short shear layer, its centreline speed kept
      // over ~2 orifice diameters and then decaying. Turbulent (a mosaic) from
      // the orifice on while it runs above ~2 m/s: a laminar 4 m/s core would
      // wrap from beam to beam into regular stripes, which no real jet shows.
      const r0 = Math.sqrt(H.AVA / Math.PI);
      const turb = clamp((H.vAS - 1.0) / 1.0, 0, 1);
      if (lAlong >= 0) {
        const w = r0 + 0.25 * lAlong, q = lPerp / w;
        if (q < 1.25 && inAo) {
          const prof = q < 0.8 ? 1 : (1.25 - q) / 0.45;
          considerFlow(true, AX_LVOT[0], AX_LVOT[1], AX_LVOT[2],
            H.vAS * prof / (1 + Math.max(0, lAlong - 0.8) / 3.5), FLOW.AS_JET,
            turb * (0.6 + 0.3 * smoothstep(0, 0.6, lAlong)));
        }
      } else if (lAlong > -1.5 && ((inLV && py > G.lv.c[1]) || inAo || lPerp < 0.9)) {
        // flow acceleration below the valve: hemispheric shells converging on the
        // orifice (v = v_orifice (r0 / (r0 + d))^2), a small aliased dome in the
        // LVOT. An LVOT PW gate 1 cm below the valve still reads the LVOT flow.
        const ac = G.A.valves.aortic.c;
        const ex = ac[0] - px, ey = ac[1] - py, ez = ac[2] - pz, dO = Math.hypot(ex, ey, ez);
        if (dO < 0.6 && dO > 1e-3) {
          // (confined to the last ~5 mm: the LVOT PW gate, which the annulus
          // swings toward late in systole, stays outside it; the streamlines
          // turn axial as they reach the orifice)
          const sp = H.vAS * (AS_R0 / (AS_R0 + dO)) ** 2 * (1 - smoothstep(0.3, 0.6, dO));
          const kr = Math.min(1, dO / 0.8) / dO;
          const ux = ex * kr + AX_LVOT[0], uy = ey * kr + AX_LVOT[1], uz = ez * kr + AX_LVOT[2], ul = Math.hypot(ux, uy, uz) || 1;
          considerFlow(true, ux / ul, uy / ul, uz / ul, sp, FLOW.AS_JET, sp > 1.3 ? turb * 0.6 : 0);
        }
      }
    } else {
    const lSpeed = H.vAoPeak * H.ejAo * gauss(lPerp, 0.7);
    considerFlow(((inLV && py > G.lv.c[1]) || inAo) && lAlong > -0.25,
      AX_LVOT[0], AX_LVOT[1], AX_LVOT[2], lSpeed, FLOW.LVOT, 0);
    }
    // subvalvular LVOT: flow accelerating toward the valve over the last ~2 cm,
    // at the LVOT velocity by continuity (forward SV through a ~2 cm LVOT): ~1 m/s
    // whether or not the valve beyond is stenotic (the PW-in-LVOT / CW-across-AV
    // pair of the continuity equation)
    if (lAlong <= -0.25 && lAlong > -2.2 && (inLV || (lPerp < 0.9 && !inRA)) && py > G.lv.c[1]) {
      considerFlow(true, AX_LVOT[0], AX_LVOT[1], AX_LVOT[2],
        H.vLvot * gauss(lPerp, 0.8) * clamp((lAlong + 2.2) / 0.9, 0, 1), FLOW.LVOT, false);
    }
  }
  // RV ejection: on the right heart's clock it starts ~10 ms before and ends
  // ~30 ms after the aortic, so it has its own gate
  if (H.ejPv > 0.001 || H.vPv > 0.005) {
    // RVOT -> pulmonary artery (centreline anchored at the pulmonic valve, on the
    // valve / trunk axis), reaching 2.4 cm down the subpulmonary tract and 3.5 cm
    // up the trunk, so no stray systolic colour appears mid-RV near the tricuspid valve.
    const rPerp = offAxisDist(px, py, pz, G.A.valves.pulmonic.c, AX_RVOT);
    const rAlong = _axis.along;
    // (the solved pulmonary flow over the RVOT area: timing, AT, notch and peak)
    const vRvot = H.vPv != null ? H.vPv : H.ejPv * 0.9;
    // The outflow profile is blunt (plug flow across the ~2.2 cm subpulmonary
    // tract and trunk): flat over the central 0.45 cm, tapering toward the walls,
    // so the stream fills the lumen and still crosses a plane that misses the
    // valve centre by ~1 cm (the PSAX-AV; review 12 repair). Never into the
    // aortic root beside it; at the stream's edge only inside the RV / PA lumen
    // (and not where the left PA runs over the left atrium).
    const rSpeed = vRvot * gauss(Math.max(0, rPerp - 0.45), 0.55);
    if (rAlong > -2.4 && rAlong < 3.5 && rPerp < 1.3 && !inAo && rSpeed > 0.01 &&
        (!_velHasBest || rSpeed > _velResult.speed) && (rPerp < 0.7 || inRvotLumen(px, py, pz, G.A, rAlong))) {
      considerFlow(true, AX_RVOT[0], AX_RVOT[1], AX_RVOT[2], rSpeed, FLOW.RVOT, false);
    }
    // The infundibulum: blood converges up the subpulmonary tract onto the valve
    // along the tract's own axis (upper infundibulum -> just below the valve),
    // a little slower than at the valve and broader. considerFlow keeps the
    // faster stream, so the PW gate below the valve still reads the trunk-axis
    // stream above.
    const iLen = infundibulumAxis(G.A);
    const iPerp = offAxisDist(px, py, pz, G.A.rvot[2], _infU);
    const iAlong = _axis.along;
    if (iAlong > -iLen - 0.5 && iAlong < 0.3 && iPerp < 1.4) {
      const iSpeed = 0.85 * vRvot * gauss(iPerp, 0.9);
      if (iSpeed > 0.01 && (!_velHasBest || iSpeed > _velResult.speed) && lumenDist(px, py, pz, G.A, 'RV') < 0) {
        considerFlow(true, _infU[0], _infU[1], _infU[2], iSpeed, FLOW.RVOT, false);
      }
    }
  }

  // --- pathological jets: a fast vena-contracta core tapering laterally and
  //     distally (Gaussian off-axis + travel decay). ---
  // (regurgitant jets only through a SHUT valve: review 12, PHY-1 / PAT-4)
  if ((path.mr || path.dilated) && !(G.valves && G.valves.mitral >= 0.1)) {
    // mitral regurgitation (primary, or functional in DCM): systolic jet LV -> LA. The region reaches from the
    // mitral annulus DEEP into the LA body, with a long travel-decay so the
    // brightest turbulent mosaic sits well within the atrium (mid-LA). The
    // downstream colour is confined STRICTLY to the LA side of the crux: the
    // sample must be inside the LA lumen (or, in the annulus->LA transition,
    // above the annulus AND on the LV/LA side of the septum) and must never fall
    // inside the RV/RA proxies — so no blue/mosaic leaks into the RV or septum.
    const cMV = G.A.valves.mitral.c;
    const aboveMV = (px - cMV[0]) * VALVES.mitral.n[0] + (py - cMV[1]) * VALVES.mitral.n[1] + (pz - cMV[2]) * VALVES.mitral.n[2] > 0.1;
    const downstream = !inRV && !inAo && iasSide(px, py, pz) < 0 && (inLA || aboveMV);
    // speed from the instantaneous LV-LA pressure drop (holosystolic, rounded);
    // jet length from the modelled regurgitant fraction. A posterior-leaflet
    // prolapse jet is eccentric and wall-hugging; functional MR (DCM) is central.
    // (the plume grows and shrinks with the instantaneous regurgitant flow)
    const qMR = H.QmrMax > 0 ? Math.sqrt(clamp(H.Qmr / H.QmrMax, 0, 1)) : 1;
    const lenMR = clamp(1.4 + 7 * (H.regurgFraction || 0) * qMR, 1.4, 5.2);
    regurgJet(px, py, pz, cMV, path.mr ? AX_MR_ECC : AX_MR, downstream, inLV, MR_JET_PARAMS,
      H.vMR != null ? H.vMR : H.ejMr * H.vMRPeak, lenMR, path.mr ? G.A : null,
      clamp(0.1 + 0.55 * (H.regurgFraction || 0), 0.15, 0.4));   // vena-contracta half-width: wider with severity
  }
  if (H.vTRPeak > 0 && !(G.valves && G.valves.tricuspid >= 0.1)) {
    // tricuspid regurgitation: systolic jet RV -> RA, firing into the RA cavity.
    // The TR peak velocity encodes RV systolic pressure: with pressure overload
    // (pulmonary hypertension) the jet runs fast (~4.3 m/s → PASP ~80 mmHg via
    // 4v²+RAP); isolated TR at normal PA pressure runs ~3.0 m/s.
    const cTV = G.A.valves.tricuspid.c;
    const aboveTV = (px - cTV[0]) * VALVES.tricuspid.n[0] + (py - cTV[1]) * VALVES.tricuspid.n[1] + (pz - cTV[2]) * VALVES.tricuspid.n[2] > 0.1;
    const downstream = !inLV && !inAo && iasSide(px, py, pz) > 0 && (inRA || aboveTV);
    // Jet size from the SOLVED circulation, not from a grade label (review 12,
    // PAT-11): the length follows the regurgitant volume per beat (~55 mL fills
    // the atrium: severe) and the instantaneous regurgitant flow; the vena
    // contracta is the radius of the effective regurgitant orifice, sqrt(Atr/pi).
    // So the functional TR of pulmonary hypertension (~15 mL through 0.1 cm2) is
    // a small fast jet and torrential TR (~60 mL) fills the atrium.
    const qTR = H.QtrMax > 0 ? Math.sqrt(clamp(H.Qtr / H.QtrMax, 0, 1)) : 1;
    const lenTR = H.QtrMax > 0 ? 1.4 + 3.2 * clamp(H.trRegurgVol / 55, 0, 1) * qTR : (JET_LEN[gradeOf(path)] || JET_LEN.severe);
    const vcTR = H.Atr > 0 ? clamp(Math.sqrt(H.Atr / Math.PI), 0.12, 0.45) : 0.38;
    regurgJet(px, py, pz, cTV, AX_TR, downstream, inRV, TR_JET_PARAMS,
      H.vTR != null ? H.vTR : H.vTRPeak * pulse(0.02, 0.48, G.phaseRV != null ? G.phaseRV : phase),
      lenTR, null, vcTR);
  }
  if (path.vsd) {
    // ventricular septal defect: LV -> RV across the septum. Flow converges on the
    // defect from the LV side (laminar), leaves it as a narrow fast jet and
    // spreads into the RV as a turbulent plume.
    let vJet = H.vVsd != null ? H.vVsd : 4.5 * pulse(0.02, 0.44, phase);
    // restrictive VSD: the systolic jet speed IS the LV-RV gradient (Bernoulli),
    // ~4.5-5 m/s. In diastole the few-mmHg gradient drives a slow to-and-fro
    // signal, well under 1 m/s (the simplified Bernoulli relation overstates such
    // low, viscous-loss-dominated flows), so it is damped below ~2 m/s.
    const avj = Math.abs(vJet);
    vJet *= 0.5 + 0.5 * smoothstep(1.5, 3.0, avj);
    const peri = path.vsd === 'perimembranous';           // else the mid-muscular defect
    const ms = peri ? membSite(G) : null;
    const core = peri ? ms.c : VSD_CORE, ax = peri ? ms.rad : AX_VSD;
    // (review 12, PHY-5) gated on the jet axis, not on a ~1 cm box round the
    // defect: laminar convergence in the LV lumen (PISA), the narrow jet through
    // the septal channel, then a turbulent plume that spreads and slows over
    // ~3 cm in the RV lumen (regurgJet: vena contracta, flame, rounded head).
    const av = Math.abs(vJet);
    if (av > 0.05) {
      const jax = vJet >= 0 ? ax : negAxis(ax);
      const dx = px - core[0], dy = py - core[1], dz = pz - core[2];
      const dn = dx * jax[0] + dy * jax[1] + dz * jax[2];
      if (dn > -1.6 && dn < 4.6 && dx * dx + dy * dy + dz * dz < 25) {
        const perp2 = dx * dx + dy * dy + dz * dz - dn * dn;
        // the channel through the septum (classified as LV blood) on either side of the core
        const channel = Math.abs(dn) < 0.9 && perp2 < VSD_RADIUS * VSD_RADIUS;
        const up = vJet >= 0 ? 'LV' : 'RV', down = vJet >= 0 ? 'RV' : 'LV';
        const downstream = dn >= 0 && (channel || lumenDist(px, py, pz, G.A, down) < 0);
        const upstream = dn < 0 && (channel || lumenDist(px, py, pz, G.A, up) < 0);
        regurgJet(px, py, pz, core, jax, downstream, upstream, VSD_JET_PARAMS, av, VSD_JET_LEN, null);
      }
    }
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
    // (within ~4 cm of the LA centre, further in a dilated atrium, whose vein
    // mouths sit further out: the solved D wave of DCM was cut off otherwise)
    const reach = Math.max(4, Math.max(G.la.r[0], G.la.r[1], G.la.r[2]) + 2);
    if (dcx * dcx + dcy * dcy + dcz * dcz < reach * reach) {
      const sp = pvFlowSpeed(H);
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
  // --- left-atrial-appendage flow -----------------------------------------
  // Along the appendage axis (tip -> ostium for emptying) inside its segments,
  // strongest on the axis; gated to within ~5 cm of the ostium.
  if (A && A.la && A.la.aa && A.la.aa.length) {
    const aa = A.la.aa, o = aa[0].a, dx = px - o[0], dy = py - o[1], dz = pz - o[2];
    if (dx * dx + dy * dy + dz * dz < 25) {
      const sp = laaFlowSpeed(phase);
      if (sp !== 0) {
        for (let i = 0; i < aa.length; i++) {
          const g = aa[i], ux = g.a[0] - g.b[0], uy = g.a[1] - g.b[1], uz = g.a[2] - g.b[2], L = Math.hypot(ux, uy, uz) || 1;
          const u = [ux / L, uy / L, uz / L];
          const perp = offAxisDist(px, py, pz, g.b, u), al = _axis.along;
          if (al < -0.2 || al > L + 0.4) continue;
          const t = al < 0 ? 0 : al > L ? 1 : al / L, r = g.r2 + (g.r1 - g.r2) * t;  // (b end: r2; a end: r1)
          if (perp > r) continue;
          const v = sp * (1 - 0.6 * (perp / r) * (perp / r));
          const dir = v >= 0 ? 1 : -1;
          considerFlow(true, u[0] * dir, u[1] * dir, u[2] * dir, v < 0 ? -v : v, FLOW.LAA_FLOW, false);
        }
      }
    }
  }
  // --- hepatic-vein flow (subcostal) ---------------------------------------
  if (A && A.ivc) {
    let hv = _hvCache.get(A);
    if (!hv) { hv = hepaticVeinSegs(A); _hvCache.set(A, hv); }
    const j = hv[0].a, dx = px - j[0], dy = py - j[1], dz = pz - j[2];
    if (dx * dx + dy * dy + dz * dz < 42) {               // within ~6.5 cm of the confluence
      const sp = hvFlowSpeed(H);
      if (sp !== 0) {
        for (let i = 0; i < hv.length; i++) {
          const v = hv[i];
          const perp = offAxisDist(px, py, pz, v.a, v.u), t = _axis.along;
          if (t < 0 || t > v.len) continue;
          const r = v.r1 + (v.r2 - v.r1) * t / v.len;
          if (perp > r) continue;
          const s = sp * gauss(perp, r * 0.9);
          const dir = s >= 0 ? -1 : 1;                     // + = toward the confluence (RA)
          considerFlow(true, v.u[0] * dir, v.u[1] * dir, v.u[2] * dir, s < 0 ? -s : s, FLOW.HV_FLOW, false);
        }
      }
    }
  }
  if (path.asd) asdFlow(px, py, pz, G, H);
  aorticFlow(px, py, pz, G, H, phase);

  return _velHasBest ? _velResult : null; // may be null (no meaningful flow)
}
