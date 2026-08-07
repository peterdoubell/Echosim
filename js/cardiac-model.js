// cardiac-model.js
// A dependency-free analytic model of the human heart used as the single source
// of truth for BOTH the 3D anatomy view and the 2D echo sector renderer.
//
// Everything lives in "heart space" (centimetres):
//   +y = towards the base (top),  -y = towards the apex (bottom)
//   +x = patient's left  (the LV / systemic side)
//   +z = anterior (towards the chest wall / transducer)
//
// The model is intentionally *schematic* rather than anatomically exact: the
// chambers are ellipsoids, the valves are thin annular disks with hinged
// leaflets, and blood flow is a piece-wise bulk-velocity field. This keeps the
// classification maths cheap enough to sample tens of thousands of points per
// frame while still reproducing the pattern-recognition cues echocardiography
// students actually rely on.

import { clamp, smoothstep, pulse, unit } from './mathutils.js';
import { anatomyParams, bodyClassify, epiDist, BODY, CFG } from './anatomy.js';
import { hemodynamics, hemoSummary, pvLoop, gradeOf, eaRatio } from './hemodynamics.js';

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
// Derive the coarse RV / aorta flow-gating proxies from CFG rather than hand-
// tuned magic numbers: the RV single-ellipsoid is the midpoint + bounding span
// of the two CFG crescent lobes, and the aortic-root ellipsoid is grown from the
// aortic-valve centre + sinus radius. Both then track the anatomy automatically.
const midOf = (a, b) => [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2];
const spanOf = (cA, rA, cB, rB) =>
  [0, 1, 2].map((i) => Math.max(rA[i], rB[i]) + Math.abs(cA[i] - cB[i]) / 2);


const BASE = {
  lv: { c: CFG.lv.c, r: CFG.lv.r, wall: CFG.lv.wall }, // LVIDd = 2*r[0] = 4.6 cm (mid-normal)
  // single ellipsoid bounding the two-lobe RV crescent (flow-gating proxy only),
  // derived from the CFG rvA/rvB lobes (midpoint centre + bounding span)
  rv: { c: midOf(CFG.rvA.c, CFG.rvB.c), r: spanOf(CFG.rvA.c, CFG.rvA.r, CFG.rvB.c, CFG.rvB.r), wall: CFG.rvWall },
  la: { c: CFG.la.c, r: [CFG.la.r, CFG.la.r, CFG.la.r] },
  ra: { c: CFG.ra.c, r: [CFG.ra.r, CFG.ra.r, CFG.ra.r] },
  // aortic-root proxy (sinuses + ascending column) for LVOT / aorta flow gating,
  // grown from the aortic-valve centre + sinus radius (CFG.aoValve / aoSinusR)
  aorta: {
    c: [CFG.aoValve[0], CFG.aoValve[1] + 0.6, CFG.aoValve[2] + 0.06],
    r: [CFG.aoSinusR + 0.35, CFG.aoSinusR + 0.65, CFG.aoSinusR + 0.35],
  },
};

// Valve definitions: a plane point + normal, an annulus radius, the two
// chambers it separates, and when it is open. Leaflets swing about the annulus.
// The tricuspid annulus is seated ~0.95 cm apical to the mitral (lower y) — the
// normal septal-leaflet offset, exaggerated here so the annular STEP at the crux
// of the heart reads clearly in A4C (the tricuspid inserts more apically on the
// septum than the mitral). Used for the A4C crux cue.
const VALVES = {
  mitral:    { c: [1.25, 1.35, -0.25], n: unit([0.18, 1, -0.28]), r: 1.45, opensInDiastole: true },
  tricuspid: { c: [-1.5, 0.4, 0.35],   n: unit([-0.12, 1, 0.2]),  r: 1.55, opensInDiastole: true },
  aortic:    { c: [0.42, 1.5, 0.14],   n: unit([0.05, 1, 0.05]),  r: 1.05, opensInDiastole: false },
  pulmonic:  { c: [-0.95, 1.95, 1.0],  n: unit([-0.2, 1, 0.55]),  r: 1.0,  opensInDiastole: false },
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
const _mFlow = unit([0.05, -1, 0.12]);        // mitral: into the LV
const _tvFlow = unit([-0.05, -1, 0.05]);      // tricuspid: into the RV
const _mAP = _perp([VALVES.aortic.c[0] - VALVES.mitral.c[0],
                    VALVES.aortic.c[1] - VALVES.mitral.c[1],
                    VALVES.aortic.c[2] - VALVES.mitral.c[2]], _mFlow);
const _aoFlow = unit([0.05, 1, 0.05]);        // aortic: up the root
const _puFlow = unit([-0.2, 1, 0.55]);        // pulmonic: up the PA
const _aoB = _basis(_aoFlow), _puB = _basis(_puFlow);
const LEAFLETS = {
  // AV valves (mitral/tricuspid): flow = downstream leaflet-extension axis; ap =
  // in-plane axis toward the longer (anterior) leaflet; com = intercommissural axis
  // (posterior scallops); cd = coaptation depth; asym = anterior/posterior ratio.
  mitral:    { flow: _mFlow, cd: 1.5,  thick: 0.075, thickP: 0.09,
               ap: _mAP, com: _cross(_mFlow, _mAP), scallops: 3, asym: 1.55 },
  // tricuspid: the septal leaflet (toward the IVS, +ap) is the shortest and least
  // mobile; the anterior/lateral leaflet is longer — so asym < 1 shortens +ap.
  tricuspid: { flow: _tvFlow, cd: 1.35, thick: 0.08,
               ap: _perp([VALVES.mitral.c[0] - VALVES.tricuspid.c[0],
                          VALVES.mitral.c[1] - VALVES.tricuspid.c[1],
                          VALVES.mitral.c[2] - VALVES.tricuspid.c[2]], _tvFlow), asym: 0.7 },
  // Semilunar valves (aortic/pulmonic): three cusps. u/w give the angular frame;
  // cusps=3 triggers the trilobed (triangular) orifice + commissural coaptation
  // seams that read as the short-axis "Mercedes" Y when the cusps shut.
  aortic:    { flow: _aoFlow, cd: 0.95, thick: 0.065, ap: null, asym: 1,
               cusps: 3, u: _aoB.u, w: _aoB.w, commOff: 0.52 },
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
// Contraction fraction k: 0 at end-diastole (full), 1 at end-systole (smallest).
function contraction(phase) {
  // ejection roughly 0.02 -> 0.34, relaxation to 0.5
  const up = smoothstep(0.0, 0.14, phase);       // rapid contraction
  const down = 1 - smoothstep(0.34, 0.5, phase); // relaxation
  return Math.min(up, down);
}
// Atrial contraction (atrial kick) late in diastole.
function atrialKick(phase) {
  return pulse(0.86, 1.0, phase);
}

// Valve opening 0..1 for a given valve at a phase.
function valveOpening(name, phase, path) {
  const v = VALVES[name];
  let o;
  if (v.opensInDiastole) {
    // AV valves: open during filling (E wave then A wave), shut in systole
    const eWave = pulse(0.5, 0.66, phase);
    const aWave = pulse(0.84, 1.0, phase);
    o = clamp(Math.max(eWave, 0.75 * aWave), 0, 1);
    o = phase < 0.48 ? 0.02 : Math.max(0.15, o); // held open through diastasis
    if (phase > 0.48 && phase < 0.86) o = Math.max(o, 0.35);
  } else {
    // semilunar valves: open during ejection
    o = smoothstep(0.02, 0.08, phase) * (1 - smoothstep(0.32, 0.37, phase));
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
const MVA_NORM = 5.0, MVA_MS = 1.8;   // mitral effective orifice area
const ET_EJECT = 0.33, T_FILL = 0.37; // systolic ejection / diastolic filling time (s)
const K_AO = 1.56, K_MV = 2.2;        // peak/mean profile factors (peak ≈ k·mean)

// Build the per-frame haemodynamic snapshot the flow field reads. Returns a
// FRESH object (never the hemodynamics() singleton) so two geometryAt() results
// — e.g. echo's end-diastolic and end-systolic geometries — never alias.
function hemoFrame(hd, path) {
  const s = hd.sum;
  // effective orifice area follows the severity grade (continuity: same SV through
  // a smaller area → higher velocity), so the Doppler peak and BSE band both track it.
  const g = gradeOf(path);
  const AVA = path.aorticStenosis ? ({ mild: 0.88, moderate: 0.71, severe: 0.53 })[g] : AVA_NORM;
  const MVA = path.mitralStenosis ? ({ mild: 2.6, moderate: 2.0, severe: MVA_MS })[g] : MVA_NORM;
  // peak jet speeds (m/s): continuity for forward flows, Bernoulli (ΔP=4v²) for MR
  const vAoPeak = K_AO * s.forwardSV / (AVA * ET_EJECT) / 100;
  const vMitPeak = K_MV * s.forwardSV / (MVA * T_FILL) / 100;
  const vMRPeak = Math.sqrt(Math.max(0, s.PlvSys - s.PlaAtSys) / 4);
  // TR jet peak (m/s) encodes RV systolic pressure: high with pressure overload
  // (pulmonary hypertension), moderate for isolated TR. Feeds the PASP estimate.
  const vTRPeak = path.tr ? (path.rvpo ? 4.3 : 3.0) : 0;
  // instantaneous 0..1 envelopes from the modelled flow waveform (carry E/A,
  // ejection and regurgitant timing) so velocityAt only multiplies scalars.
  const ejAo = s.QaoMax > 0 ? clamp(hd.Qao / s.QaoMax, 0, 1) : 0;
  const ejMv = s.QmvMax > 0 ? clamp(hd.Qmv / s.QmvMax, 0, 1) : 0;
  const ejMr = s.QmrMax > 0 ? clamp(hd.Qmr / s.QmrMax, 0, 1) : 0;
  return {
    Plv: hd.Plv, Pao: hd.Pao, Vlv: hd.Vlv, Pla: hd.Pla,
    Qao: hd.Qao, Qmv: hd.Qmv, Qmr: hd.Qmr, rho: hd.rho, kv: hd.kv,
    EDV: s.EDV, ESV: s.ESV, SV: s.SV, EF: s.EF, forwardSV: s.forwardSV,
    regurgFraction: s.regurgFraction, gradient: s.gradient,
    PlvSys: s.PlvSys, PaoSys: s.PaoSys, PaoDia: s.PaoDia,
    vAoPeak, vMitPeak, vMRPeak, vTRPeak, ejAo, ejMv, ejMr,
  };
}

// ---------------------------------------------------------------------------
// Live geometry for a given phase + pathology. Returns per-chamber ellipsoids
// (centre + radii) plus wall thicknesses, so both renderers stay consistent.
// ---------------------------------------------------------------------------
export function geometryAt(phase, path = {}) {
  const k = contraction(phase);      // kinematic timing (wall-thickening phase)
  const kick = atrialKick(phase);

  // PHYSICS COUPLING: the lumped-parameter circulation solves the LV pressure–
  // volume loop. hd.rho = Vlv/EDV is the ABSOLUTE cavity-volume ratio (1 at
  // end-diastole, ESV/EDV at end-systole), so it carries the true ejection
  // fraction (low in DCM/ischemia). The cavity semi-axes are scaled so their
  // product tracks rho EXACTLY (sShort²·sLong = rho), with the long axis
  // shortening less than the short axis; hence the measured EF/LVIDs fall out of
  // the modelled volume rather than a prescribed phase curve. hd.kv (normalised)
  // still times the wall-thickening / RV pulsation.
  const hd = hemodynamics(phase, path);
  const rho = hd.rho;                 // Vlv / EDV (absolute cavity-volume ratio)
  const kv = hd.kv;                   // normalised contraction fraction (timing)

  let lvScale = 1.0, lvWallMul = 1.0;
  if (path.dilated) { lvScale = CFG.dilatedScale; lvWallMul = 0.82; } // dilated cavity
  if (path.lvh || path.aorticStenosis) { lvWallMul = 1.7; } // hypertrophy

  const aLong = 0.30;                       // long axis shortens less than short axis
  const sLong = Math.pow(rho, aLong);
  const sShort = Math.pow(rho, (1 - aLong) / 2); // sShort²·sLong = rho (volume-exact)
  const lvR = [
    BASE.lv.r[0] * lvScale * sShort,
    BASE.lv.r[1] * lvScale * sLong,
    BASE.lv.r[2] * lvScale * sShort,
  ];
  // apex-anchored longitudinal contraction: the LV centre shifts apically as the
  // long axis shortens (apex fixed), so the coarse proxy tracks the SDF and the
  // mitral annulus descends the full shortening (MAPSE). Volume-preserving.
  const lvR1d = BASE.lv.r[1] * lvScale;         // end-diastolic long semi-axis
  const lvApexY = BASE.lv.c[1] - lvR1d;         // fixed apex
  const lvCy = lvApexY + lvR1d * sLong;         // anchored centre y
  // wall thickens as the short axis shrinks (approx. muscle-volume conservation)
  const lvWall = BASE.lv.wall * lvWallMul * clamp(1 / sShort, 1, 1.7);

  const rvR = BASE.rv.r.map((r, i) => {
    const shorten = i === 1 ? 0.12 : 0.26;
    return r * (1 - kv * shorten);
  });
  const rvWall = BASE.rv.wall * (1 + k * 0.5);

  // atria: fill through systole (reservoir), empty during diastole, kick at end
  // (reuse the already-computed contraction fraction k rather than recomputing)
  const atrialFill = 0.85 + 0.35 * (1 - k) - 0.35 * kick;
  const laR = BASE.la.r.map((r) => r * (path.dilated ? CFG.laDilation : 1) * atrialFill);
  const raR = BASE.ra.r.map((r) => r * atrialFill);

  // anatomical SDF bundle (built once here so the chordae apparatus below can
  // hang off the same apex-anchored papillary tips + valve plane).
  const A = anatomyParams(kv, kick, path, { sShort, sLong, lvWall });

  return {
    phase, k, kick, kv,
    // instantaneous lumped-parameter circulation state (PV loop, pressures,
    // stroke volume, EF, jet peak velocities) — read by the echo measurements,
    // the flow field and the UI PV-loop panel. A fresh snapshot (not the
    // hemodynamics() singleton) so distinct geometries never alias.
    hemo: hemoFrame(hd, path),
    // Coarse per-chamber ellipsoid proxies — used by the velocity field, the
    // on-image labels and the LVIDd/EF measurements. Tissue classification and
    // the 3D surface use the anatomical SDF (G.A) below, not these.
    lv: { c: [BASE.lv.c[0], lvCy, BASE.lv.c[2]], r: lvR, wall: lvWall },
    rv: { c: BASE.rv.c, r: rvR, wall: rvWall },
    la: { c: BASE.la.c, r: laR },
    ra: { c: BASE.ra.c, r: raR },
    aorta: BASE.aorta,
    // anatomical SDF parameter bundle — classify() and the 3D marching-cubes
    // surface both sample this for a topologically correct heart. The same
    // volume-exact LV scaling (sShort/sLong) is threaded through so the SDF
    // cavity, the coarse proxy above and the measured EF all track the PV loop.
    A,
    // subvalvular apparatus: chordae tendineae fanning from the papillary tips
    // to the mitral leaflet free edges (built from A so they track the annulus).
    chordae: buildChordae(A, path),
    valves: {
      mitral: valveOpening('mitral', phase, path),
      tricuspid: valveOpening('tricuspid', phase, path),
      aortic: valveOpening('aortic', phase, path),
      pulmonic: valveOpening('pulmonic', phase, path),
    },
  };
}

// ---------------------------------------------------------------------------
// Chordae tendineae. Thin fibrous strands from each papillary-muscle tip to the
// mitral leaflet free edges — echo shows them as fine linear echoes tethering
// the leaflets into the LV cavity (and, tellingly, they go SLACK in diastole and
// pull TAUT in systole). Modelled as a small fan of capsule segments so the
// tissue sampler can paint them over the LV blood pool. Built once per frame off
// the apex-anchored papillary tips (A.pap[i].b) and the mitral coaptation zone.
function buildChordae(A, path) {
  const segs = [];
  // free-edge anchor points for a valve: the coaptation zone + two points spread
  // along the leaflet (ap) axis, so the fan spans the leaflet free edge.
  const anchorsFor = (name) => {
    const v = VALVES[name], lf = LEAFLETS[name];
    const cy = (A.axial && (name === 'mitral' || name === 'tricuspid'))
      ? A.axial.apexY + (v.c[1] - A.axial.apexY) * A.axial.lsy : v.c[1];
    const f = lf.flow, ap = lf.ap, d = lf.cd * 0.88, sp = 0.4;
    const b = [v.c[0] + f[0] * d, cy + f[1] * d, v.c[2] + f[2] * d];
    return [b,
      [b[0] + ap[0] * sp, b[1] + ap[1] * sp, b[2] + ap[2] * sp],
      [b[0] - ap[0] * sp, b[1] - ap[1] * sp, b[2] - ap[2] * sp]];
  };
  // mitral: both LV papillary tips -> anterior/posterior leaflet free edges
  const rMit = path.mitralStenosis ? 0.07 : 0.035;  // fine strands; fused in MS
  const mAnch = anchorsFor('mitral');
  for (const p of A.pap) for (const q of mAnch) segs.push({ a: p.b, b: q, r: rMit });
  // tricuspid: the RV anterior papillary (fused to the moderator band) -> the
  // tricuspid anterior leaflet free edge, completing the septum→band→PM→leaflet chain
  if (A.rvPap) {
    const tAnch = anchorsFor('tricuspid');
    for (const q of tAnch) segs.push({ a: A.rvPap.b, b: q, r: 0.033, tv: true });
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
// region-name → central angle (model frame); 'apical' handled separately.
const RWMA_ANGLE = {
  septal: Math.PI, anteroseptal: Math.PI * 0.75, anterior: Math.PI / 2,
  lateral: 0, inferior: -Math.PI / 2, posterior: -Math.PI / 2,
};

// how "affected" (0..1) a segment is by the RWMA spec {region, sev}
function segAffected(seg, rw) {
  if (!rw) return 0;
  if (rw.region === 'apical') {
    return rw.sev * (seg.level >= 2 ? 1 : seg.level === 1 ? 0.35 : 0.05);
  }
  const ang = RWMA_ANGLE[rw.region] != null ? RWMA_ANGLE[rw.region] : Math.PI;
  if (seg.level === 3) return rw.sev * 0.3;           // apex: partial
  let d = seg.angle - ang;
  while (d > Math.PI) d -= 2 * Math.PI;
  while (d < -Math.PI) d += 2 * Math.PI;
  const t = d / 1.15;                                  // matches RWMA half-width
  const levelF = seg.level === 2 ? 0.8 : 1;            // apical territory slightly less
  return rw.sev * Math.exp(-t * t) * levelF;
}

// Returns { segments: [{name, angle, level, strain}], gls } for a pathology.
let _rsPath = null, _rsVal = null;
export function regionalStrain(path = {}) {
  if (path === _rsPath && _rsVal) return _rsVal;
  let gf = 1;                                           // global contractility factor
  if (path.dilated) gf = 0.42;
  else if (path.aorticStenosis && path.lvh) gf = 0.85;
  else if (path.lvh) gf = 0.9;
  const rw = path.rwma ? {
    region: (typeof path.rwma === 'object' ? path.rwma.region : (typeof path.rwma === 'string' ? path.rwma : 'septal')) || 'septal',
    sev: (typeof path.rwma === 'object' && path.rwma.severity != null) ? path.rwma.severity : 0.85,
  } : null;
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
};

export function classify(px, py, pz, G, path = {}) {
  const A = G.A;

  // --- pericardial effusion: a DEPENDENT echo-free crescent, not a rind ---
  // Fluid fills the gap just outside the epicardial surface, widest posteriorly
  // (-z) / inferiorly (low y) behind the LV and near-zero anteriorly/superiorly.
  if (path.effusion) {
    const de = epiDist(px, py, pz, A); // signed distance to the whole-heart epicardium
    if (de > 0) {
      const post = clamp(0.5 - (pz - A.lv.c[2]) * 0.55, 0, 1);
      const infer = clamp(0.5 - (py - A.lv.c[1]) * 0.32, 0, 1);
      const gap = 0.15 + 1.7 * clamp(Math.max(post, infer), 0, 1);
      if (de < gap && py < A.lv.c[1] + 1.6) return cls(TISSUE.PERICARDIUM, 0.02);
    }
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
  const body = bodyClassify(px, py, pz, A);
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

function aortaLumen(px, py, pz, G) {
  const a = G.aorta;
  // root as ellipsoid + vertical column up to the arch
  if (ellip(px, py, pz, a.c, a.r) <= 1) return true;
  // ascending column
  const dx = px - a.c[0], dz = pz - a.c[2];
  const rad = Math.hypot(dx, dz);
  if (py > a.c[1] && py < a.c[1] + 2.4 && rad < a.r[0] * 0.82) return true;
  return false;
}

// Separate singleton from _clsResult so a valve hit and a subsequent non-valve
// classify() return can never share (alias) the same object.
const _valveResult = { tissue: TISSUE.VALVE, echo: 0 };

function valveTissue(px, py, pz, G, path) {
  // apex-anchored contraction descends the AV valve planes with the annulus: a
  // valve at end-diastolic height y maps to apexY + (y−apexY)·lsy (same material
  // map as the LV wall). Mitral follows the LV; aortic sits at the LV base too.
  const ax = G.A && G.A.axial;
  for (const name in VALVES) {
    const v = VALVES[name];
    const lf = LEAFLETS[name];
    const cy = (ax && (name === 'mitral' || name === 'aortic'))
      ? ax.apexY + (v.c[1] - ax.apexY) * ax.lsy : v.c[1];
    // vector from annulus centre to the sample, and its split into the leaflet
    // axial (downstream) component `a` and the in-plane radial vector.
    const dx = px - v.c[0], dy = py - cy, dz = pz - v.c[2];
    const f = lf.flow;
    const a = dx * f[0] + dy * f[1] + dz * f[2];          // axial: 0 at hinge, +downstream
    const rx = dx - a * f[0], ry = dy - a * f[1], rz = dz - a * f[2];
    const rad = Math.hypot(rx, ry, rz);                   // in-plane radius from axis
    if (rad > v.r + 0.2) continue;                        // outside the annulus footprint
    const open = G.valves[name];                          // 0 shut … 1 fully open
    // Anterior/posterior asymmetry + posterior scalloping/thickness.
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
    const s = t * t * (3 - 2 * t);
    const shellR = v.r + (edgeR - v.r) * s;               // curved radius at this depth
    // Thin membrane tapering to a fine free edge (thickest at the annular base).
    let th = thick * (0.45 + 0.55 * (1 - t) * (1 - t));
    // Nodulus of Arantius: a fibrous thickening at the CENTRE of each semilunar
    // cusp's free edge, where the three cusps meet — bulk up the free-edge tip at
    // the cusp centre so the coaptation shows the three little nodules.
    let nod = 0;
    if (isCusp && t > 0.7) nod = 0.11 * ((t - 0.7) / 0.3) * (1 - cptr);
    th += nod;
    let hit = Math.abs(rad - shellR) <= th;
    // Semilunar commissural coaptation seams: three radial lines meeting centrally
    // as the cusps shut — the short-axis "Mercedes" Y, and the closure line in LAX.
    if (!hit && isCusp) {
      const closed = 1 - open;
      const angDist = Math.abs(kh) / 3;                   // angular distance to a commissure
      if (closed > 0.25 && a > 0.4 * cd && rad < v.r * 0.9 && angDist * rad < 0.07) hit = true;
    }
    if (!hit) continue;
    let echo = 0.8 + 0.05 * (1 - t) + (nod > 0.03 ? 0.12 : 0); // nodule reads brighter
    if ((name === 'aortic' && path.aorticStenosis) || (name === 'mitral' && path.mitralStenosis)) {
      echo = 0.98; // calcified / restricted → bright, thickened
    }
    _valveResult.tissue = TISSUE.VALVE;
    _valveResult.echo = echo;
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
const AX_MITRAL_IN = unit([0.05, -1, 0.12]);    // mitral annulus -> LV apex
const AX_TRICUSPID_IN = unit([-0.05, -1, 0.05]);// tricuspid annulus -> RV apex
const AX_LVOT = unit([-0.12, 1, 0.05]);         // LVOT -> aorta
const AX_RVOT = unit([0.0, 1, 0.35]);           // RVOT -> pulmonary artery
const AX_MR = unit([0.1, 1, -0.55]);            // mitral -> into the LA
const AX_TR = unit([-0.1, 1, 0.15]);            // tricuspid -> into the RA
const AX_VSD = unit([-1, 0.1, 0.15]);           // LV -> RV across the septum
const AX_ASD = unit([-1, 0.05, 0.3]);           // LA -> RA across the septum
const VSD_CORE = [-1.1, -0.4, 0.3];  // on the interventricular septum (LV septal endocardium ~x -1.1)
const ASD_CORE = [-0.4, 2.6, -0.35]; // on the interatrial septum, between LA and RA
// Regurgitant-jet spread regions (ellipsoid centres): MR reaches deep into the
// LA body, TR into the RA body. Used to bound each jet's colour to its receiving
// atrium so no mosaic leaks across the septum into the other side of the heart.
const MR_REGION = [1.1, 2.75, -0.95];  // LA body (mitral regurgitation target)
const TR_REGION = [-1.85, 2.55, 0.3];  // RA body (tricuspid regurgitation target)

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
  const jet = env != null ? env : pulse(0.02, 0.36, phase);
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

    // tricuspid inflow -> RV apex (right heart not lumped-modelled; use the
    // diastolic filling envelope for correct timing at a lower normal velocity)
    const tPerp = offAxisDist(px, py, pz, VALVES.tricuspid.c, AX_TRICUSPID_IN);
    const tAlong = _axis.along;
    const tW = 1.35 * (1 - 0.5 * clamp(tAlong / 4.5, 0, 1));
    const tShape = gauss(tPerp, tW) * Math.exp(-clamp(tAlong, 0, 6) / 4.5);
    const tSpeed = 0.55 * H.ejMv * tShape;
    considerFlow(inRV && tAlong > -0.3, AX_TRICUSPID_IN[0], AX_TRICUSPID_IN[1], AX_TRICUSPID_IN[2],
      tSpeed, FLOW.TRICUSPID_IN, false);
  }

  // --- systolic ejection: a THIN column hugging the outflow centreline only,
  //     NOT a whole-ventricle fill. Present only during ejection (0.04-0.34).
  const ejec = pulse(0.04, 0.34, phase);
  if (ejec > 0.001) {
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

    // RVOT -> pulmonary artery (centreline anchored at the pulmonic valve).
    // Tightened to the true outflow centreline: a narrow Gaussian, clipped to
    // the downstream (+along) side and to the high subpulmonary tract (py > 0.8)
    // so no stray systolic colour blob appears mid-RV near the tricuspid valve.
    const rPerp = offAxisDist(px, py, pz, VALVES.pulmonic.c, AX_RVOT);
    const rAlong = _axis.along;
    const rSpeed = H.ejAo * 0.9 * gauss(rPerp, 0.5); // systolic timing from the ejection envelope
    considerFlow(inRV && py > 0.8 && rAlong > -0.7, AX_RVOT[0], AX_RVOT[1], AX_RVOT[2],
      rSpeed, FLOW.RVOT, false);
  }

  // --- pathological jets: a fast vena-contracta core tapering laterally and
  //     distally (Gaussian off-axis + travel decay). ---
  if (path.mr) {
    // mitral regurgitation: systolic jet LV -> LA. The region reaches from the
    // mitral annulus DEEP into the LA body, with a long travel-decay so the
    // brightest turbulent mosaic sits well within the atrium (mid-LA). The
    // downstream colour is confined STRICTLY to the LA side of the crux: the
    // sample must be inside the LA lumen (or, in the annulus->LA transition,
    // above the annulus AND on the LV/LA side of the septum) and must never fall
    // inside the RV/RA proxies — so no blue/mosaic leaks into the RV or septum.
    const region = ellip(px, py, pz, MR_REGION, [1.7, 3.0, 1.6]) <= 1;
    const downstream = region && !inRV && !inRA && px > 0.2 && (inLA || py > 1.4);
    // envelope + peak from the modelled regurgitant flow / LV-LA gradient: the
    // regurgitant fraction sets how much colour fills the LA; the systolic
    // LV-LA pressure difference sets the ~5 m/s vena-contracta velocity.
    regurgJet(px, py, pz, phase, VALVES.mitral.c, AX_MR, downstream, inLV, MR_JET_PARAMS,
      H.ejMr, H.vMRPeak);
  }
  if (path.tr) {
    // tricuspid regurgitation: systolic jet RV -> RA, firing into the RA cavity.
    // The TR peak velocity encodes RV systolic pressure: with pressure overload
    // (pulmonary hypertension) the jet runs fast (~4.3 m/s → PASP ~80 mmHg via
    // 4v²+RAP); isolated TR at normal PA pressure runs ~3.0 m/s.
    const trPeak = path.rvpo ? 4.3 : 3.0;
    const region = ellip(px, py, pz, TR_REGION, [1.7, 2.7, 1.6]) <= 1;
    const downstream = region && (inRA || py > 1.1);
    regurgJet(px, py, pz, phase, VALVES.tricuspid.c, AX_TR, downstream, inRV, TR_JET_PARAMS, null, trPeak);
  }
  if (path.vsd) {
    // ventricular septal defect: LV -> RV across septum, systole, turbulent core
    const jet = pulse(0.02, 0.4, phase);
    const region = ellip(px, py, pz, VSD_CORE, [1.0, 1.6, 1.2]) <= 1;
    const perp = offAxisDist(px, py, pz, VSD_CORE, AX_VSD);
    const along = _axis.along;
    const shape = gauss(perp, 0.55) * Math.exp(-clamp(Math.abs(along), 0, 3) / 2.2);
    const speed = jet * 4.5 * shape;                      // restrictive VSD peak ~4.5 m/s
    considerFlow(region, AX_VSD[0], AX_VSD[1], AX_VSD[2], speed, FLOW.VSD_JET, jet > 0.15);
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
