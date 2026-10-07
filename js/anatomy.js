// anatomy.js — an anatomically-structured signed-distance model of the heart
// and the structures an echocardiographer sees around it.
//
// Both the 2D echo cross-section and the 3D marching-cubes surface sample this
// one field, so the two views are always consistent.
//
// Heart space (cm) is aligned with the LEFT-VENTRICULAR LONG AXIS:
//   +y  apex -> base (the LV long axis; apex at the bottom)
//   +x  toward the LV lateral wall (patient's left, tipped posterior-superior)
//   +z  anterior (toward the chest wall / RV, tipped superior)
// The LV long axis is the line x = z = 0.
//
// Geometry is built from LANDMARKS placed at adult reference dimensions (ASE /
// BSE chamber quantification, CT anatomy) rather than tuned per view, and every
// standard window is then derived from those landmarks (views.js). The spatial
// relationships a sonographer relies on are explicit: aorto-mitral fibrous
// continuity, septal-aortic continuity, the apical offset of the tricuspid
// septal leaflet, an RVOT that wraps anterior to the aortic root, atria that sit
// on their annuli and share a flat septum, the coronary sinus in the posterior
// AV groove, the descending aorta behind the left atrium, the pericardium, and
// the liver / IVC / hepatic veins that form the subcostal window.

import { smin, smax, ssub, sdSphere, sdEllipsoid, sdCapsule, sdRoundCone, sdCylinder, sdFrustum } from './sdf.js';
import { smoothstep } from './mathutils.js';

// ---- BSE/ASE normal reference ranges -------------------------------------
// Adult echocardiographic normal ranges (ASE/EACVI 2015 chamber quantification,
// men where sex-specific). Exported so the Measurements UI can draw the "normal
// band" next to each value and the anatomy audit (tools/verify-anatomy.mjs) can
// check the model against them. Centimetres except EF (%), areas (cm^2).
export const REF = {
  lviddNormal:  [4.2, 5.8],  // LV internal diameter, end-diastole
  lvidsNormal:  [2.5, 4.0],  // LV internal diameter, end-systole
  lvWallNormal: [0.6, 1.0],  // IVSd / PWd (septal / posterior wall thickness)
  efNormal:     [52, 72],    // ejection fraction (%)
  lvLengthNormal: [7.2, 9.4], // LV length, apex -> mitral annular plane (A4C, ED)
  laNormal:     [3.0, 4.0],  // LA antero-posterior diameter (PLAX, end-systole)
  raAreaNormal: [10, 18],    // RA area (<= 18 cm^2 normal)
  aoAnnulusNormal: [2.0, 2.9], // aortic annulus
  aoRootNormal: [2.9, 3.7],  // aortic sinus-of-Valsalva diameter
  aoAscNormal:  [2.2, 3.6],  // proximal ascending aorta
  lvotNormal:   [1.8, 2.4],  // LVOT diameter
  paNormal:     [1.5, 2.5],  // main pulmonary-artery diameter
  rvd1Normal:   [2.5, 4.1],  // RV basal diameter (RVD1)
  rvd3Normal:   [5.9, 8.3],  // RV base-apex length (RVD3)
  rvWallNormal: [0.1, 0.5],  // RV free-wall thickness
  tvOffsetNormal: [0.3, 1.2], // TV septal-leaflet apical displacement vs mitral
  ivcNormal:    [1.2, 2.1],  // IVC diameter (expiration)
  csNormal:     [0.4, 1.0],  // coronary sinus diameter
  dtaNormal:    [1.8, 2.8],  // descending thoracic aorta diameter
};

// ---- tiny vector helpers (build-time only; the hot path stays scalar) -------
const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const mul = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
const mad = (a, b, s) => [a[0] + b[0] * s, a[1] + b[1] * s, a[2] + b[2] * s];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const unit = (a) => { const l = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0] / l, a[1] / l, a[2] / l]; };
const lerp3 = (a, b, t) => [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t];

// ---- the heart in the chest ------------------------------------------------
// The long axis runs from the base (right, posterior, superior) to the apex
// (left, anterior, inferior): ~45 deg off the sagittal plane and ~30 deg below
// the horizontal (apex direction 0.61 L, 0.61 A, 0.50 I). The rotation ABOUT the
// long axis is fixed clinically: heart +z is body-anterior projected
// perpendicular to the long axis — the parasternal short-axis "12 o'clock" (the
// anterior wall) — so +x is the lateral wall (3 o'clock), -x the septum
// (9 o'clock) and -z the inferior wall (6 o'clock). These are the patient's body
// axes in that heart space. They place the body-referenced structures — great
// vessels, pulmonary veins, descending aorta, diaphragm, liver — and set which
// way fluid is DEPENDENT.
const HEART_ROLL = 45 * Math.PI / 180;
const _roll = (v) => [v[0] * Math.cos(HEART_ROLL) - v[2] * Math.sin(HEART_ROLL), v[1], v[0] * Math.sin(HEART_ROLL) + v[2] * Math.cos(HEART_ROLL)];
// unrolled axes: for short structures attached to the heart whose ostia were
// laid out in heart space (the pulmonary veins) — they rotate WITH the heart
const BL0 = [0.632, -0.612, -0.476], BP0 = [0.0, 0.612, -0.791], BS0 = [0.775, 0.500, 0.386];
const BL = _roll([0.632, -0.612, -0.476]);  // patient left
const BP = _roll([0.0, 0.612, -0.791]);     // patient posterior
const BS = _roll([0.775, 0.500, 0.386]);    // patient superior
export const BODY_AX = {
  L: BL, R: mul(BL, -1), P: BP, A: mul(BP, -1), S: BS, I: mul(BS, -1),
};

// ---- landmarks (end-diastolic reference geometry) --------------------------
// LV: a truncated prolate profile about x = z = 0. ED cavity length 8.3 cm
// (apex endocardium -> annular plane), widest 2 cm below the annulus with an
// endocardial semi-axis of 2.45 (LVIDd 4.9 cm), closing in over the top 2 cm to
// the mitral inflow and LV outflow funnels. Integrated volume ~115 mL — the
// circulation model's EDV — so geometry and haemodynamics agree in absolute
// terms, not just as ratios. Length/width 1.7 (normal sphericity > 1.5).
const LVP = {
  apexY: -8.3,       // endocardial apex (ED)
  len: 8.3,          // apex -> annular plane
  sEq: 6.3,          // equator (widest) height above the apex
  sTop: 7.95,        // top of the body profile (funnels carry it to the valves)
  b: 2.45,           // equatorial endocardial radius (LVIDd 4.9)
  rTop: 2.0,         // body radius where the base closes in
  wall: 0.9,         // basal / mid wall (IVSd = PWd = 0.9)
  apexWallFrac: 0.7, // apical wall thins to ~0.63 cm
};

// Mitral annulus centre: posterior to the LV long axis, so the LV outflow tract
// occupies the anteroseptal base and the A4C plane (apex-mitral-tricuspid)
// passes BEHIND the LVOT. Annulus 3.25 cm in A4C.
const M0 = [0.15, 0.0, -0.9];
const MV_R = 1.62;
// Aortic valve: anterior and rightward of the mitral, sharing the aorto-mitral
// curtain (annuli ~0.1 cm apart). The root axis is tipped ~49 deg anterior of the
// LV long axis — an aorto-septal angle of ~131 deg (normal 126 +/- 6 deg; it
// narrows with age and in the "sigmoid septum") — while keeping the RVOT, the
// aortic root and the LA on ONE line perpendicular to the root: the classic
// Ao/LA M-mode line of the PLAX view.
const A0 = [-1.52, 0.45, 0.88];              // ~10 o'clock (anteroseptal) of the LV axis
// the root axis lies IN the long-axis plane through the LV axis and the aortic
// valve (so the PLAX/A3C plane shows the whole root and ascending aorta)
const AO_TILT = 49 * Math.PI / 180;
const D_AO0 = unit([A0[0], 0, A0[2]]);
const U_AO = unit(add(mul([0, 1, 0], Math.cos(AO_TILT)), mul(D_AO0, Math.sin(AO_TILT))));
// the Ao/LA M-mode line: perpendicular to the root, pointing away from the chest wall
const AO_LA_DIR = unit(add(mul([0, 1, 0], Math.sin(AO_TILT)), mul(D_AO0, -Math.cos(AO_TILT))));
const AO = { annR: 1.15, sinusR: 1.65, stjR: 1.4, ascR: 1.5 };  // d 2.3 / 3.3 / 2.8 / 3.0
// LV outflow tract: from the anteroseptal LV cavity up to the aortic annulus,
// d 2.2 cm, running close to the septum.
const LVOT0 = [-0.95, -2.1, 0.55];
const LVOT_R = 1.1;
// Angle-dependent mitral annular excursion. The anterior annulus is the fibrous
// aorto-mitral curtain: it is continuous with the aortic annulus and descends with
// the root (~55 % of the long-axis shortening), while the muscular posterolateral
// annulus moves the full amount — so lateral MAPSE exceeds septal, the annular
// saddle deepens in systole and the curtain keeps its length. Weight at angle phi
// from the aorto-mitral direction: w = 1 - (1 - MV_W_CURTAIN)·((1 + cos phi)/2)^4
// (~0.55 at the curtain, ~0.7 at the A4C septal hinge, ~1 posterolaterally).
const MV_AP = unit([A0[0] - M0[0], 0, A0[2] - M0[2]]);   // horizontal, toward the aorta
const MV_W_CURTAIN = 0.55;
// mean weight round the annulus (the centre's share): mean of ((1+cos)/2)^4 = C(8,4)/4^4
const MV_W_MEAN = 1 - (1 - MV_W_CURTAIN) * 70 / 256;
function mvWeight(cosPhi) {
  const g = (1 + cosPhi) / 2, g2 = g * g;
  return 1 - (1 - MV_W_CURTAIN) * g2 * g2;
}
// Tricuspid annulus: right of the crux and 0.85 cm APICAL to the mitral — the
// septal-leaflet offset that identifies the morphological RV on A4C (Ebstein
// anomaly exaggerates it). Annulus 3.4 cm; inflow directed at the RV apex.
const T0 = [-4.25, -0.85, -0.1];
const TV_R = 1.7;

// PSAX-at-aortic-valve frame: looking up the root axis from the apex, anterior
// (toward the probe) up, patient-left to the right. Used to place the pulmonary
// valve and RVOT in the positions the PSAX-AV view shows them.
// root frame (cusps, PV, coronary ostia) is rigid with the heart: unrolled anterior
const _ANT0 = mul(BP0, -1);
const E_ANT = unit(sub(_ANT0, mul(U_AO, dot(_ANT0, U_AO))));
const E_SCR = unit(cross(U_AO, E_ANT));
const avFrame = (deg, r, h) => {
  const a = deg * Math.PI / 180;
  return add(add(A0, mul(U_AO, h)), add(mul(E_SCR, Math.cos(a) * r), mul(E_ANT, Math.sin(a) * r)));
};
// Pulmonary valve: anterior, leftward and ~1.2 cm higher than the aortic valve
// (~1-2 o'clock in PSAX-AV). The RVOT arcs over the front of the aortic root from
// the tricuspid side (~10 o'clock) to it — the RV "wraps" the aorta.
const PV0 = avFrame(58, 3.05, 1.25);
const PV_R = 1.1;                                    // annulus d 2.2
// Main PA runs superiorly, posteriorly and to the left to bifurcate ~3.5 cm
// above the LA (the RPA then passes behind the ascending aorta, over the LA roof).
const U_PA = unit(add(add(mul(BS, 0.7), mul(BP, 0.9)), mul(BL, 0.25)));
const PA_BIF = mad(PV0, U_PA, 3.8);
// the RPA passes to the right BEHIND the ascending aorta and SVC, beneath the arch
const U_RPA = unit(add(add(mul(BODY_AX.R, 0.85), mul(BP, 0.45)), mul(BS, -0.35)));
const U_LPA = unit(add(add(mul(BL, 0.6), mul(BP, 0.7)), mul(BS, 0.25)));
// (the infundibulum is separated from the LVOT / aortic root by the muscular
// outlet septum, so on PLAX its centre lies ~3.4 cm in front of the aortic-root
// centre, clear of the basal LV)
const RVOT_PTS = [avFrame(132, 3.4, -1.0), avFrame(90, 4.1, 0.5), mad(PV0, U_PA, -0.25)];

// Atria. Each sits on its (moving) AV annulus with a roof fixed in the chest, so
// the systolic descent of the AV plane STRETCHES them: that is how the atrial
// reservoir fills, and it keeps the atria attached to the valves at every phase.
// Lengths/widths are the END-SYSTOLIC (maximum) values measured clinically.
const MAPSE_REF = 1.55;                              // cm, from the circulation (see anatomyParams)
const TAPSE_REF = 2.5;
// The LA rises from the mitral annulus roughly in line with the LV long axis
// (A4C) and lies directly behind the aortic root on the Ao/LA M-mode line (PLAX,
// PSAX-AV) — "posterior" in the body is +y as much as -z in heart space.
const A_LA = unit([0.07, 1, 0.13]);
const A_RA = unit([-0.12, 1, 0.08]);
const LA_FLOOR0 = add(M0, [-0.15, 0, -0.05]);
const RA_FLOOR0 = add(T0, [0.45, 0.1, -0.2]);
const LA_ES = { len: 5.0, w1: 1.55, w2: 1.75 };       // A4C major 5.0, (cut) minor ~4.0, PLAX AP ~3.4
// the LA body sits back along the Ao/LA line, behind the root (PLAX AP ~3.3 cm)
const LA_POST = 0.8;
const RA_ES = { len: 4.0, w1: 2.3, w2: 1.9 };        // A4C major ~4.6, (cut) minor ~3.7
// Roofs follow 35 % of their annulus's excursion (the atria are not pinned),
// placed so the end-systolic floor-to-roof length is exactly the ES value above.
const ROOF_FOLLOW = 0.35;
const LA_ROOF = mad(add(LA_FLOOR0, [0, -(1 - ROOF_FOLLOW) * MAPSE_REF, 0]), A_LA, LA_ES.len);
// tricuspid hinge: fraction of TAPSE at the annulus centre / septal rim, and the
// in-plane direction from the septal toward the lateral (free-wall) rim
const TV_SEPT_FRAC = 0.62, TV_CENTRE_FRAC = 0.609;   // septal rim moves with the mitral annulus (fibrous skeleton); RA floor follows the old centre
// the RA body lies right-posterior of the aortic root: the non-coronary sinus
// indents only its anteromedial wall, it does not sit in the chamber's middle
const RA_POST = 1.0;
const RA_LAT = 0.2;
const RA_AWAY = (() => {
  const mid = mad(RA_FLOOR0, A_RA, 2.0), root = mad(A0, U_AO, 1.0);
  const v = sub(mid, root);
  return unit(sub(v, mul(A_RA, dot(v, A_RA))));
})();
const RA_ROOF = mad(add(RA_FLOOR0, [0, -(1 - ROOF_FOLLOW) * TV_CENTRE_FRAC * TAPSE_REF, 0]), A_RA, RA_ES.len);
// Interatrial septum: a plane through the crux and the posterior aortic root
// (the non-coronary sinus abuts it), facing the RA. The atria overlap across it
// and are CUT by it, so they share one flat septal wall — thin over the fossa
// ovalis, thicker at the limbus.
const IAS_P = [-2.05, 1.3, -0.1];
const IAS_N = unit([-0.92, -0.05, 0.4]);
// fossa ovalis: the thin central membrane, low and posterior on the septum
// (posterior to the aortic root, where the A4C / subcostal planes cross the septum)
const FOSSA0 = [-2.33, 1.65, -0.8];
// horizontal direction within the septal plane, toward the posterior-right
// (septal frame: IAS_N across the septum, SEPT_V back along it)
const SEPT_V = unit(cross(IAS_N, unit(sub(BS, mul(IAS_N, dot(BS, IAS_N))))));
const IAS = {
  tLimbus: 0.16, tFossa: 0.05, rFossa: 0.6,
  pen: { la: 0.5, ra: 1.3 },                          // septal-face overlap across the plane (cm)
  fossaC: sub(FOSSA0, mul(IAS_N, dot(sub(FOSSA0, IAS_P), IAS_N))),
};
// Caval orifices in the septal frame (IAS_N across the septum, SEPT_V back along
// it, from the fossa centre), cm. svcSd/ivcSd: distance on the RA side of the
// septal plane; svcV/ivcV: offset along SEPT_V (the SVC then slides further back
// until it clears the aortic root by aoGap); svcUp: SVC orifice height above the
// reference RA roof; ivcTilt: the IVC runs down and slightly back (SEPT_V share).
// sv*: the sinus venarum capsule between the orifices (fractions svT of the way
// SVC -> IVC, svBulge off the septum, its upper end svFwd forward along the
// septum so the RA covers the septum up to the SVC); flare*: the SVC orifice
// funnel toward the tricuspid valve.
const CAVA = {
  svcSd: 1.0, svcV: 0.5, svcUp: 1.2, svcR: 0.8, aoGap: 0.15,
  ivcSd: 0.95, ivcV: -1.0, ivcTilt: 0.25,
  svBulge: 0.2, svR: 1.2, svT: [0.15, 0.7], svFwd: 0.5,
  flareUp: 1.0, flareLen: 1.5, flareR: 1.1,
};

// Aorta beyond the root: ascending aorta, arch and descending thoracic aorta.
// The DTA runs cranio-caudally behind the LA, slightly left of the spine; on
// PLAX it is the round echo-free structure just behind the atrioventricular
// groove — the landmark that separates a pericardial effusion (tracks ANTERIOR
// to it) from a left pleural effusion (posterior to it).
// beyond the sino-tubular junction the ascending aorta turns cranially (it lies
// just right of the sternum, not far out to the right)
// beyond the STJ the ascending aorta runs almost straight up, just right of the
// sternum (a body-fixed direction, independent of the root's tilt)
const U_ASC = unit(add(add(mul(BS, 0.85), mul(BODY_AX.R, 0.15)), mul(BODY_AX.A, 0.35)));
const AO_STJ = mad(A0, U_AO, 2.1);
const AO_ASC_TOP = mad(AO_STJ, U_ASC, 4.2);
// (posterolateral-left of the oesophagus, which lies directly against the LA)
const DTA_P = add(add(add(LA_FLOOR0, mul(BP, 4.7)), mul(BL, 0.22)), mul(BS, -0.6));
const DTA_TOP = mad(DTA_P, BS, 7.0);
const DTA_BOT = mad(DTA_P, BS, -9.0);
const AO_ARCH_MID = add(lerp3(AO_ASC_TOP, DTA_TOP, 0.5), mul(BS, 2.0));
const DTA_R = 1.12;                                   // d 2.25
// The arch as a smooth Catmull-Rom curve (ascending top -> arch -> descending),
// tapering from the ascending calibre to the descending, with its three head and
// neck branches (brachiocephalic, left common carotid, left subclavian) rising
// from the top — the suprasternal long-axis landmarks.
const ARCH_SEGS = (() => {
  const ctl = [mad(AO_ASC_TOP, U_ASC, -1.0), AO_ASC_TOP, AO_ARCH_MID, DTA_TOP, mad(DTA_TOP, BS, -1.5)];
  const cr = (p0, p1, p2, p3, t) => {
    const t2 = t * t, t3 = t2 * t;
    return [0, 1, 2].map((k) => 0.5 * (2 * p1[k] + (p2[k] - p0[k]) * t +
      (2 * p0[k] - 5 * p1[k] + 4 * p2[k] - p3[k]) * t2 + (3 * p1[k] - p0[k] - 3 * p2[k] + p3[k]) * t3));
  };
  const pts = [];
  for (let i = 1; i < ctl.length - 2; i++) for (let j = 0; j < 6; j++) pts.push(cr(ctl[i - 1], ctl[i], ctl[i + 1], ctl[i + 2], j / 6));
  pts.push(DTA_TOP);
  const segs = [];
  for (let i = 0; i < pts.length - 1; i++) {
    const f = (i + 0.5) / (pts.length - 1);
    segs.push([pts[i], pts[i + 1], AO.ascR * 0.95 + (DTA_R - AO.ascR * 0.95) * f]);
  }
  // head and neck vessels from the top of the arch
  const top = pts[Math.round(pts.length * 0.45)];
  const along = unit(sub(DTA_TOP, AO_ASC_TOP));
  for (const [off, r] of [[-1.3, 0.55], [0.0, 0.38], [1.2, 0.42]]) {
    const o = mad(mad(top, along, off), BS, 0.6);
    segs.push([o, mad(mad(o, BS, 3.5), along, off * 0.3), r]);
  }
  return segs;
})();

// Body-level references for the diaphragm / liver (subcostal window). The
// heart's diaphragmatic surface — the LV and RV inferior walls — is flat and
// runs ALONG the long axis, lying on the central tendon; so locally the
// diaphragm is a plane containing the long axis, facing inferoseptally (body
// inferior, projected perpendicular to the long axis). DIAPH_H is the height of
// the lowest point of the heart on that normal, plus the pericardium. The liver
// fills the right upper quadrant beneath it, and its left lobe reaches under the
// heart: the window the subcostal beam images through.
// The standard A4C is aimed A4C_TILT cm posterior (heart -z) of the crux, off
// the LVOT and aortic root (views.js uses the same tilt).
const A4C_TILT = 0.15;
const DIAPH_N = unit(add(BODY_AX.I, mul(BODY_AX.A, 0.0)));  // central tendon: near-horizontal
// Diaphragm height along DIAPH_N: just below the heart's most inferior
// epicardium (computed once from the ED geometry, so it follows the heart's
// orientation in the body rather than a hand-tuned number).
let DIAPH_H = null;
const DIAPH_BELOW = -0.18;                           // tendon plane sits 1.8 mm into the sac's lowest point
function diaphH() {
  if (DIAPH_H !== null) return DIAPH_H;
  const A = anatomyParams(0, 0, {}, null, 0);
  let best = -1e9;
  for (let x = -7; x <= 7; x += 0.2) for (let y = -10; y <= 6; y += 0.2) for (let z = -7; z <= 7; z += 0.2) {
    const h = x * DIAPH_N[0] + y * DIAPH_N[1] + z * DIAPH_N[2];
    if (h > best - 0.2 && epiDist(x, y, z, A) < 0 && h > best) best = h;
  }
  DIAPH_H = best + DIAPH_BELOW;
  return DIAPH_H;
}
// Gastric contact for the transgastric window: the probe sits on the inner surface
// of the stomach wall (GASTRIC.p, facing GASTRIC.n toward the heart), below the
// inferior wall. The stomach itself is ONE organ, LIV_LOBE.stomach (left upper
// quadrant, under the left dome, behind the left lobe), fitted so its lumen
// surface passes through p with its normal close to n; its wall is t thick.
const GASTRIC = (() => {
  const n = unit([0.1, 0, 1]);
  return { n, p: mad([0, -3.6, 0], n, -4.6), t: 0.45 };
})();
const DIAPH_WRAP = 6.5;                              // how far the diaphragm rises to meet the sac (cm)
const MIDLINE = -2.0;                                  // p . L of the body mid-sagittal plane

// Everything the rest of the app needs to place valves, flows, labels and views.
export const LM = {
  M: M0, MV_R, MV_AP, A: A0, U_AO, AO, LVOT0, LVOT_R, T: T0, TV_R, PV: PV0, PV_R, U_PA,
  PA_BIF, U_RPA, U_LPA, RVOT: RVOT_PTS, E_ANT, E_SCR,
  LV: LVP, apex: [0, LVP.apexY, 0],
  LA_ROOF, RA_ROOF, A_LA, A_RA, LA_FLOOR0, RA_FLOOR0, LA_ES, RA_ES,
  IAS_P, IAS_N, IAS, AO_STJ, AO_ASC_TOP, AO_ARCH_MID, DTA_P, DTA_TOP, DTA_BOT, DTA_R,
  GASTRIC, DIAPH_N, A4C_TILT, ANT_HEART: _ANT0, A4C_Z: (x, y) => onA4C(x, y), MIDLINE, MAPSE_REF, TAPSE_REF, TV_SEPT_FRAC,
};

// Coarse per-chamber sizes kept for the flow/measurement proxies in
// cardiac-model.js (single source of truth, no silent drift).
const CFG = {
  lv: { c: [0, LVP.apexY + LVP.len / 2, 0], r: [LVP.b, LVP.len / 2, LVP.b], wall: LVP.wall },
  // RV body: widest at ~9:30 o'clock (PSAX), so the crescent runs from the
  // anterior insertion (~12) round the septum to the inferior insertion (~7:30)
  rvBody: { c: [-2.22, -3.35, 0.6], r: [4.5, 4.05, 3.7] },
  rvWall: 0.38,      // RV free-wall thickness (normal <= 0.5 cm)
  dilatedScale: 1.34, // DCM: LVIDd 6.6 cm (severe), more spherical
  dilatedLong: 1.12,
  laDilation: 1.2,   // LA enlargement in DCM
};
// Four-chamber plane (through the apex, mitral and tricuspid centres), as z(x, y).
const A4C_N = (() => {
  const ap = [0, LVP.apexY - LVP.wall * LVP.apexWallFrac, 0];
  const crux = add(mul(add(M0, T0), 0.5), [0, 0, -A4C_TILT]);
  const n = unit(cross(sub(crux, ap), sub(M0, T0)));
  return { n, p: ap };
})();
function onA4C(x, y) {
  const { n, p } = A4C_N;
  return p[2] - (n[0] * (x - p[0]) + n[1] * (y - p[1])) / n[2];
}
// moderator band's free-wall insertion (anterior papillary muscle base)
const MOD_FW = [-6.2, -4.75, onA4C(-6.2, -4.75)];
const TV_LAT = (() => {
  const n0 = unit(sub(T0, [-2.3, CFG.rvBody.c[1] - CFG.rvBody.r[1] + 0.6, 1.9]));
  const v = sub(T0, M0);
  return unit(sub(v, mul(n0, dot(v, n0))));
})();
export { CFG };

// ---- atrial phase function ------------------------------------------------
// The atria are RESERVOIRS, so their volume curve runs roughly ANTI-phase to the
// ventricles': they are distending while the ventricles eject, and emptying
// while the ventricles fill. Three overlapping mechanical roles, with phase 0 at
// AV-valve closure (QRS / onset of ventricular systole):
//
//   reservoir (0.00 -> 0.51)  mitral + tricuspid shut, so pulmonary-venous and
//                             caval inflow distends the atrium; descent of the
//                             AV plane during ejection adds to it. Atrial volume
//                             therefore PEAKS at AV-valve opening — the instant
//                             the ventricle is at its smallest.
//   conduit   (0.51 -> 0.69)  AV valves open (the circulation opens the mitral at
//                             ~0.50, after an ~85 ms IVRT): rapid early-diastolic emptying
//                             (the E wave), then a near-flat diastasis.
//   booster   (0.86 -> 1.00)  atrial systole (the A wave) empties the atrium to
//                             its minimum, reached exactly at the next AV-valve
//                             closure.
//
// The amplitudes satisfy A_RES === A_COND + A_BOOST, so the curve is exactly
// periodic — f(1) === f(0) — with no step or slope break across the wrap.
//
// Calibration: these are LINEAR volume-equivalent scale factors, so volume goes
// as f^3. Peak 1.17 / minimum 0.89 gives a total LA emptying fraction of
// 1 - 0.89^3/1.17^3 ~= 56 %, split into a passive (conduit) emptying fraction of
// ~38 % and an active (booster) emptying fraction of ~30 % — normal adult values.
const A_MIN = 0.89, A_RES = 0.28, A_COND = 0.13, A_BOOST = 0.15;
const A_PEAK = A_MIN + A_RES;

export function atrialFill(phase) {
  const p = phase - Math.floor(phase);          // wrap into [0,1)
  const reservoir = smoothstep(0.00, 0.51, p);
  const conduit   = smoothstep(0.51, 0.69, p);
  const booster   = smoothstep(0.89, 1.00, p);
  return A_MIN + A_RES * reservoir - A_COND * conduit - A_BOOST * booster;
}

// An oriented ellipsoid anchored between a floor point on the annulus and a
// fixed roof. Its length is whatever the anchors dictate; the two widths are set
// so the volume follows the atrial curve (V ∝ fill^3) exactly:
//   len·w1·w2 = (fill/peak)^3 · lenES·w1ES·w2ES.
const FOSSA_PULL = 1;
const ATR_GAIN = 1.0;   // uniform (whole-chamber) share of the volume change
const ATR_REMODEL = 1.6;   // width exponent of the dilatation scale
const ATR_FREE = 0.5;   // free-wall share: compression along the septal normal
const ATR_LONG = 0.5;   // roof share: compression along the long axis toward the annulus
function atrium(floor, roof, es, fill, lateral, widen, post = 0, postDir = null) {
  const ax = sub(roof, floor);
  const len = Math.hypot(ax[0], ax[1], ax[2]);
  const a = mul(ax, 1 / len);
  const e1 = unit(cross(a, [0, 0, 1]));
  const e2 = cross(e1, a);
  // the septal anchor (fixed septum) and the vein / appendage mouths add a
  // near-constant volume, so the free walls must move more than the bare
  // ellipsoid for the DRAWN chamber to reach the modelled emptying fractions
  const vol = Math.pow(fill / A_PEAK, 3 * ATR_GAIN);
  // chronic atrial dilatation (pressure/volume overload) enlarges the chamber in all
  // three directions: the roof rises away from the annulus and the walls bow out
  const grow = widen > 1 ? Math.pow(widen, ATR_REMODEL) : widen;
  const w = Math.sqrt(vol * es.len / len) * grow;
  // a dilating atrium balloons laterally (e1 lies in the four-chamber plane) and
  // lengthens, rather than growing backwards out of the 4C plane
  const lenG = len * (widen > 1 ? widen : 1);
  const lat = widen > 1 ? Math.pow(widen, 0.6) : 1;
  const c = add(add(mad(floor, a, lenG / 2), mul(e1, lateral * grow)), mul(postDir || e2, postDir ? post : -post));
  // the septum stays put while the free wall moves: extra emptying compresses the
  // chamber toward the septal plane (applied in sdAtrium about that plane)
  const f = Math.pow(fill / A_PEAK, 3 * ATR_FREE);
  const g = Math.pow(fill / A_PEAK, 3 * ATR_LONG);
  return { c, a, e1, e2, rl: lenG / 2 + 0.15, r1: es.w1 * w * lat, r2: es.w2 * w, len: lenG, f, g, fl: floor };
}

// The interatrial septum is a shared wall that barely moves: the atrial volume
// change is carried by the free walls. So at every phase each atrium is
// stretched toward the septal plane (its free walls stay where they are) until
// its septal face crosses the plane by `pen` cm (sgn +1 for the LA, on the -n
// side; -1 for the RA); the septal cut then leaves one flat, broad septum
// between the chambers instead of an extracardiac gap.
function anchorToSeptum(E, sgn, pen) {
  const n = IAS_N;
  const an = dot(E.a, n), e1n = dot(E.e1, n), e2n = dot(E.e2, n);
  const h = Math.sqrt((E.r1 * e1n) ** 2 + (E.rl * an) ** 2 + (E.r2 * e2n) ** 2); // support along n
  const sd = dot(sub(E.c, IAS_P), n) * sgn;                 // centre depth, <0 on own side
  E.sn = mul(n, sgn);
  E.sw = Math.max(0, pen - (sd + h));
}

// Build the per-phase parameter bundle. k = contraction 0..1, kick = atrial
// kick 0..1, path = pathology flags. `mech` carries the volume-exact LV scaling
// from the lumped-parameter circulation (sShort/sLong/lvWall) so the SDF cavity
// tracks the modelled PV loop; when omitted a kinematic scaling is used (keeps
// anatomyParams callable standalone). `phase` drives the atrial reservoir curve;
// it cannot be recovered from k (two-valued in phase) so it is passed explicitly.
export function anatomyParams(k, kick, path = {}, mech = null, phase = 0) {
  // Remodelling: when coupled to the circulation the LV size follows the modelled
  // EDV (eccentric dilatation in DCM / chronic MR, a smaller cavity in concentric
  // LVH / AS / MS) and the LA follows the chronic LA pressure; standalone, the
  // pathology flags give representative values.
  let lvScale = 1.0, lvLong = 1.0, lvWallMul = 1.0, laScale = 1.0, contract = 1.0;
  if (path.dilated) { contract = 0.26; lvScale = CFG.dilatedScale; lvLong = CFG.dilatedLong; lvWallMul = 0.82; laScale = CFG.laDilation; }
  if (path.lvh || path.aorticStenosis) lvWallMul = 1.55;
  if (mech && mech.lvScale) { lvScale = mech.lvScale; lvLong = mech.lvLong; }
  if (mech && mech.laScale) laScale = mech.laScale;

  // RV pressure/volume overload (pulmonary hypertension): the RV dilates and the
  // interventricular septum flattens toward the LV, so the LV reads D-shaped in
  // PSAX (the "D-sign"). rvpoScale enlarges the RV; septFlat clips the LV septum.
  // Volume overload (TR, secundum ASD) dilates the RV and RA too, without the
  // pressure-overload septal flattening; the RA enlarges with any of the three.
  const rvpoScale = (path.rvpo ? 1.3 : 1) * (path.asd ? 1.22 : path.tr && !path.rvpo ? 1.18 : 1);
  const raScale = path.rvpo ? 1.3 : path.tr ? 1.35 : path.asd ? 1.25 : 1;
  const septFlat = path.rvpo ? 0.42 : 0;

  // LV cavity scale: volume-exact from the circulation when coupled.
  const sS = mech ? mech.sShort : 1 - k * 0.28 * contract;
  const sL = mech ? mech.sLong : 1 - k * 0.19 * contract;
  const lvWall = mech ? mech.lvWall : LVP.wall * lvWallMul * (1 + k * 0.5 * contract);

  // Apex-anchored longitudinal contraction: the apex stays put and the base
  // (mitral annulus) descends by the full long-axis shortening — MAPSE. A
  // remodelled ventricle keeps its base on the fibrous skeleton and grows (or
  // shrinks) toward the apex.
  const apexY = LVP.apexY - (lvLong - 1) * LVP.len;
  const lsy = sL;
  const axialMap = (y0) => apexY + (y0 - LVP.apexY) * lsy * lvLong;
  const dM = axialMap(M0[1]) - M0[1];                      // posterolateral annular displacement (<= 0 in systole)
  // the annulus narrows ~12 % in systole, almost all of it in the muscular
  // posterior annulus: the fibrous anterior rim stays on the aorto-mitral curtain,
  // so the centre drifts toward the aorta as the radius shrinks
  const mvR0 = MV_R * (path.dilated ? 1.15 : 1);           // annular dilatation in DCM
  const mvR = mvR0 * (1 - 0.12 * k);
  const mvShift = 0.85 * (mvR0 - mvR);
  // the centre descends by the mean of the angle-dependent excursion (mitralLift
  // adds each sector's share on top); the base level is the posterolateral one
  const Mlive = [M0[0] + MV_AP[0] * mvShift, M0[1] + MV_W_MEAN * dM, M0[2] + MV_AP[2] * mvShift];
  const Mbase = [M0[0], M0[1] + dM, M0[2]];
  // aortic root moves less than the mitral annulus (tethered to the arch), and
  // the aorto-mitral curtain moves with it
  const Alive = [A0[0], A0[1] + MV_W_CURTAIN * dM, A0[2]];
  // the LVOT's lower end is part of the contracting LV body, so its offset from
  // the long axis follows the short-axis scaling (the basal septum moves in)
  const LVOTlive = [LVOT0[0] * sS, axialMap(LVOT0[1]), LVOT0[2] * sS];

  // RV: apex-anchored too. The tricuspid annulus HINGES rather than translating:
  // its lateral rim descends toward the RV apex by TAPSE (~2.1 cm) while the
  // septal rim, tethered to the fibrous skeleton, moves ~45 % of that — so the
  // annulus centre drops ~0.73 TAPSE and the annular plane tilts toward the free
  // wall in systole. (TAPSE is, by definition, the lateral-annulus excursion.)
  // the right heart runs on its own, slightly offset clock (cardiac-model rvPhase:
  // T1 after M1, RV ejection starting earlier and ending later than the LV's)
  const kR = mech && mech.kRV != null ? mech.kRV : k;
  const rvS = 1 - kR * 0.24;                               // RV short-axis shortening
  const tapse = TAPSE_REF * kR;
  const tvDrop = TV_CENTRE_FRAC * tapse;
  const Tlive = [T0[0], T0[1] - tvDrop, T0[2]];
  const PVlive = mad(PV0, [0, 1, 0], -0.25 * tapse);
  const rvB = CFG.rvBody;
  // a dilated LV displaces the RV crescent outward rather than crushing it
  const push = (lvScale - 1) * LVP.b * 0.9;
  const rvDir = unit([rvB.c[0], 0, rvB.c[2]]);
  const rvApexY = rvB.c[1] - rvB.r[1] - (lvLong - 1) * LVP.len * 0.5;
  const rvRy = rvB.r[1] * rvpoScale - tapse * 0.5 + (lvLong - 1) * LVP.len * 0.25;
  const rvW = rvpoScale * (1 + (lvScale - 1) * 0.35);
  const rv = {
    c: [rvB.c[0] - (rvpoScale - 1) * 1.2 + rvDir[0] * push, rvApexY + rvRy, rvB.c[2] + (rvpoScale - 1) * 0.6 + rvDir[2] * push],
    r: [rvB.r[0] * rvS * rvW, rvRy, rvB.r[2] * rvS * rvW],
  };
  // annulus normal, facing the RA, tilted by the lateral-minus-septal excursion
  const tvN = unit(mad(unit(sub(Tlive, [-2.3, rvApexY + 0.6, 1.9])), TV_LAT, (1 - TV_SEPT_FRAC) * tapse / (2 * TV_R)));
  const tvLat = mad(Tlive, unit(sub(TV_LAT, mul(tvN, dot(TV_LAT, tvN)))), TV_R);   // lateral rim (TAPSE point)
  const inflow = { c: mad(Tlive, tvN, -1.6), r: [2.25 * rvS * rvpoScale, 1.8, 2.1 * rvS * rvpoScale] };
  const rvot = [
    mad(RVOT_PTS[0], [0, 1, 0], -0.25 * tapse),
    mad(RVOT_PTS[1], [0, 1, 0], -0.25 * tapse),
    mad(RVOT_PTS[2], [0, 1, 0], -0.25 * tapse),
  ];
  const rvotR = [1.3 * (1 - kR * 0.18) * rvpoScale, 1.1 * (1 - kR * 0.12) * rvpoScale, PV_R * 1.02];

  // Atria: sit on the live annuli, fixed roofs, volume from the reservoir curve.
  const aFill = atrialFill(phase);
  // the RA empties through the tricuspid, so its reservoir curve follows the
  // right heart's clock
  const aFillR = mech && mech.phaseRV != null ? atrialFill(mech.phaseRV) : aFill;
  // (the LA body bulges posteriorly, behind the aortic root, beyond its annulus)
  const la = atrium(add(LA_FLOOR0, [0, dM, 0]), add(LA_ROOF, [0, ROOF_FOLLOW * dM, 0]), LA_ES, aFill, -0.1, laScale, LA_POST, AO_LA_DIR);
  const ra = atrium(add(RA_FLOOR0, [0, -tvDrop, 0]), add(RA_ROOF, [0, -ROOF_FOLLOW * tvDrop, 0]), RA_ES, aFillR, RA_LAT, raScale, RA_POST, RA_AWAY);
  // As the atria empty they shrink about the free walls, not the septum: pull
  // each centre (within the septal plane) toward the fossa so the swept septal
  // face still covers it.
  for (const [E, f] of [[la, aFill], [ra, aFillR]]) {
    const pull = Math.min(1, FOSSA_PULL * (1 - f / A_PEAK));
    let t = sub(FOSSA0, E.c); t = sub(t, mul(IAS_N, dot(t, IAS_N)));
    E.c = mad(E.c, t, pull);
  }
  anchorToSeptum(la, 1, IAS.pen.la);
  anchorToSeptum(ra, -1, IAS.pen.ra);
  ra.f = Math.sqrt(ra.f);   // the thin-walled RA empties less by free-wall collapse (keeps a chamber at ED)
  ra.g = 1;                 // ...and its roof (the caval junctions) stays put
  const iasT = path.asd ? -0.12 : IAS.tFossa;

  const papY0 = -5.3, papY1 = -2.3;                         // papillary base/tip heights (ED)
  // Papillary muscles: ANTEROLATERAL (~4 o'clock in PSAX) and POSTEROMEDIAL
  // (~8 o'clock) — both in the posterior half of the LV, where the mitral
  // leaflets they tether sit. Each is a cone whose broad, trabeculated base is
  // embedded in the wall at the junction of the apical and mid thirds; only the
  // tip stands free, ~2 cm below the annulus under its commissure. The bases move
  // with the apex-anchored wall; the muscles CONTRACT in systole, so each tip
  // follows its commissural annulus down and the annulus-to-tip distance (the
  // chordal length) stays constant, and they thicken with the wall.
  // Systolic thickening (with the wall, capped at 1.4x) is carried by the free
  // head; the embedded, trabeculated base barely fattens and grows into the wall
  // (its axis moves out by the added radius), so the muscles stand out more in
  // systolic short axis without pinching the mid cavity into an hourglass.
  const papThick = Math.min(1.4, Math.max(1, lvWall / (LVP.wall * lvWallMul)));
  const papBase = 1 + (papThick - 1) * 0.15;
  const out = (papBase - 1) * 0.62;                        // added base radius
  const papAt = (ang, rb, rt) => {
    const c = Math.cos(ang), s = Math.sin(ang);
    // the tip descends with the annular sector above it (mvWeight at its angle)
    const tx = c * rt + M0[0] * 0.3 - M0[0], tz = s * rt + M0[2] * 0.3 - M0[2];
    const w = mvWeight((tx * MV_AP[0] + tz * MV_AP[2]) / Math.hypot(tx, tz));
    return {
      a: [c * (rb * sS * lvScale + out), axialMap(papY0 - (lvLong - 1) * 1.5), s * (rb * sS * lvScale + out)],
      b: [c * rt * sS * lvScale + M0[0] * 0.3, apexY + (papY1 - LVP.apexY) * lvLong + w * dM, s * rt * sS * lvScale + M0[2] * 0.3],
    };
  };
  const alpm = papAt(-0.52, 2.2, 1.8), pmpm = papAt(-2.62, 2.2, 1.75);

  const lvLive = {
      c: [0, apexY + LVP.len * lvLong * lsy * 0.5, 0],
      ax: 0, az: 0, apexY,
      sR: sS * lvScale, sL: sL * lvLong, sS, lvScale,
      sRDia: lvScale, sLDia: lvLong,
      wall: lvWall,
      wallDia: LVP.wall * lvWallMul,
      rwma: buildRwma(path),
      septFlat,
      M: Mlive, A: Alive, lvot: LVOTlive,
      mvK: dM, base: Mbase[1],             // annular excursion (see mitralLift) and posterolateral base level
      // systolic thickening goes inward at the apex: the endocardial apex creeps
      // basally by the apical thickening, so the epicardial apex stays put
      apexLift: Math.max(0, lvWall - LVP.wall * lvWallMul) * LVP.apexWallFrac,
      apexShift: apexY - LVP.apexY,        // remodelling: the apex (and the apical window) moves out
      mvR,                                 // the mitral annulus narrows ~25 % in area in systole
      lvotR: LVOT_R * sS,                  // lower LVOT contracts with the body; the top stays at the annulus
      // coarse live semi-axes (for the RWMA angular frame + legacy callers)
      r: [LVP.b * sS * lvScale, LVP.len * 0.5 * sL * lvLong, LVP.b * sS * lvScale],
  };
  const aoLive = buildAorta(Alive, path);
  const cav = (() => {
    // The caval junctions are fixed to the mediastinum (the IVC to the diaphragm),
    // so they are placed on a reference RA (the normal-sized end-diastolic
    // chamber) rather than the beating one, against the resting aortic root; a
    // dilated RA grows around them rather than lifting the SVC orifice away.
    const rr = atrium(RA_FLOOR0, RA_ROOF, RA_ES, atrialFill(0), RA_LAT, 1, RA_POST, RA_AWAY);
    const hS = Math.hypot(rr.r1 * dot(rr.e1, BS), rr.rl * dot(rr.a, BS), rr.r2 * dot(rr.e2, BS));
    // Both caval orifices sit ~1 cm on the RA side of the interatrial septum,
    // the IVC in line below the fossa ovalis (its Eustachian valve points at
    // it) and the SVC above and behind it, against the right-posterior wall of
    // the ascending aorta: the SVC -> fossa -> IVC axis of the ME bicaval view.
    // Heights come from the RA roof and floor.
    const seat = (q, sd, v) => mad(mad(q, IAS_N, sd - dot(sub(q, IAS_P), IAS_N)), SEPT_V, v - dot(sub(q, IAS.fossaC), SEPT_V));
    let sa = seat(mad(rr.c, BS, hS + CAVA.svcUp), CAVA.svcSd, CAVA.svcV);
    let away = sub(rr.c, T0); away = unit(sub(away, mul(BS, dot(away, BS))));  // from the TV, level
    // d 1.8 cm; plethoric (~2.5 cm) with raised RA pressure (pulmonary hypertension,
    // severe TR), when its wider wall keeps it a little further off the septum
    const ivcR = (path.rvpo || path.tr) ? 1.25 : 0.9;
    const ia = seat(mad(mad(mad(rr.c, BS, -(hS + 0.4)), away, 1.0), BP, 1.0), Math.max(CAVA.ivcSd, ivcR + 0.05), CAVA.ivcV);
    // the SVC abuts the aorta without indenting it: slide the orifice back along
    // the septum until its first 3 cm clear the aortic wall
    const aoA = { ao: buildAorta(A0, path) };
    const sDir = unit(add(mul(BS, 4.5), mul(BP, 0.3)));
    const gap = (q) => {
      let g = 1e9;
      for (let t = 0; t <= 3.0; t += 0.25) {
        const p = mad(q, sDir, t);
        g = Math.min(g, dAOroot(p[0], p[1], p[2], aoA) - aoA.ao.wall - CAVA.svcR - CAVA.aoGap);
      }
      return g;
    };
    for (let it = 0; it < 40 && gap(sa) < 0; it++) sa = mad(sa, SEPT_V, 0.05);
    // the junctions follow only the roof's share of the tricuspid excursion, and
    // only along the caval axis (they slide, they do not swing out of line)
    const su = unit(sub(BS, mul(IAS_N, dot(BS, IAS_N))));
    const tether = mul(su, -ROOF_FOLLOW * tvDrop * su[1]);
    sa = add(sa, tether);
    return {
      svc: { a: sa, b: mad(sa, sDir, 4.5), r: CAVA.svcR },     // d 1.6
      ivc: { a: add(ia, tether), b: mad(add(ia, tether), unit(mad(BODY_AX.I, SEPT_V, CAVA.ivcTilt)), 11.0), r: ivcR },
    };
  })();
  // the band's free-wall end moves in with the RV free wall (short-axis shortening)
  const modFW = [rvB.c[0] + (MOD_FW[0] - rvB.c[0]) * rvS, MOD_FW[1] + tapse * 0.4, MOD_FW[2]];
  const params = {
    k, kick, contract, phase,
    lv: lvLive,
    rv, inflow, rvot, rvotR, rvWall: CFG.rvWall * (path.rvpo ? 1.8 : 1) * (1 + kR * 0.45), rvCarve: 0.28,
    la: { ...la, aa: laaSegs(la) }, laScale,
    pv: buildPV(la),
    // sinus venarum: the smooth venous back of the RA between the caval orifices,
    // reaching forward to the septum below the SVC (CAVA.sv*)
    ra: { ...ra, aa: raaSegs(ra), sv: { a: mad(mad(lerp3(cav.svc.a, cav.ivc.a, CAVA.svT[0]), IAS_N, CAVA.svBulge), SEPT_V, -CAVA.svFwd), b: mad(lerp3(cav.svc.a, cav.ivc.a, CAVA.svT[1]), IAS_N, CAVA.svBulge), r: CAVA.svR },
      svcFlare: { a: mad(cav.svc.a, unit(sub(cav.svc.b, cav.svc.a)), CAVA.flareUp), b: mad(cav.svc.a, unit(sub(Tlive, cav.svc.a)), CAVA.flareLen) } },
    crista: null,                          // (seated on the live RA wall below)
    eustachian: buildEustachian(ra),
    // SVC: enters the POSTERO-superior RA roof, beside (not in line with) the
    // RA long axis, and runs cranially behind the ascending aorta
    // the venae cavae enter the RA's superior and posterior-inferior poles, in
    // line with each other (the ME bicaval axis); the IVC orifice lies ~3 cm
    // from the tricuspid annulus, across the cavotricuspid isthmus
    svc: cav.svc, ivc: cav.ivc,
    ias: { p: IAS_P, n: IAS_N, tLimbus: IAS.tLimbus, tFossa: iasT, rFossa: IAS.rFossa, fc: IAS.fossaC },
    cs: buildCS(Mbase, add(RA_FLOOR0, [0, -tvDrop, 0])),
    // aortic root: three sinuses of Valsalva around the valve, sino-tubular
    // junction, ascending aorta, arch and descending thoracic aorta (one vessel).
    ao: aoLive,
    // epicardial coronary arteries in their grooves (LM, LAD, LCx, RCA, PDA)
    cor: buildCoronaries(lvLive, aoLive, Tlive, tvN),
    pa: {
      pv: PVlive,
      main: { a: mad(PVlive, U_PA, -0.1), b: PA_BIF, r1: 1.1 * (1 + (path.rvpo ? 0.35 : 0)), r2: 1.05 * (1 + (path.rvpo ? 0.35 : 0)) },
      branch: [[PA_BIF, mad(PA_BIF, U_RPA, 4.8), 0.8], [PA_BIF, mad(PA_BIF, U_LPA, 2.8), 0.68]],
      wall: 0.16,
    },
    valves: {
      mitral: { c: Mlive }, aortic: { c: Alive }, tricuspid: { c: Tlive, n: tvN, lat: tvLat }, pulmonic: { c: PVlive },
    },
    pap: [
      { a: alpm.a, b: alpm.b, r: 0.5, rb: 0.62 * papBase, rt: 0.3 * papThick },
      { a: pmpm.a, b: pmpm.b, r: 0.5, rb: 0.62 * papBase, rt: 0.3 * papThick },
    ],
    // moderator band: crosses the lower third of the RV from the septum to the
    // base of the anterior papillary muscle on the free wall.
    // It lies in the four-chamber plane (the A4C / subcostal landmark for the
    // morphological RV and a classic pseudo-mass), so z is taken on that plane.
    mod: { a: [-2.9, axialMap(-5.3), onA4C(-2.9, -5.3)], b: modFW, r: 0.3 },
    // RV anterior papillary muscle — continuous with the moderator band, rising
    // toward the tricuspid, whose anterior-leaflet chordae it anchors.
    rvPap: { a: modFW, b: [-4.3, -2.7 - tapse * 0.6, 1.2], r: 0.33 },
    axial: { apexY, lsy: lsy * lvLong },
    trab: !path.noTrab,
  };
  params.crista = buildCrista(params);
  return params;
}

// Left atrial appendage: a narrow oval ostium on the anterolateral LA opening
// into a neck, then a multi-lobed body that hooks back on itself in the left AV
// groove over the circumflex — the "chicken-wing" morphology (~48%). The narrow
// neck is why it is THE site of thrombus in atrial fibrillation.
function laaSegs(la) {
  // ostium on the ANTERIOR LA wall near the annulus (so the A2C, rotated ~15 deg
  // lateral of 12 o'clock, cuts the ostium and neck on its anterior side); the
  // body runs forward and down in the left AV groove over the circumflex,
  // beneath the PA trunk, and its tip hooks back on itself
  const u = unit(add(add(mul(la.e1, -0.2), mul(la.e2, 0.95)), mul(la.a, -0.2)));
  const s = 1 / Math.hypot(dot(u, la.e1) / la.r1, dot(u, la.a) / la.rl, dot(u, la.e2) / la.r2);
  const o = mad(la.c, u, s * 0.9);
  const p1 = mad(o, unit([0.3, -0.6, 1]), 1.1);
  const p2 = mad(p1, unit([0.6, -0.7, 0.8]), 0.9);
  const p3 = mad(p2, unit([-0.2, -0.5, 0.7]), 0.8);
  return [
    { a: o, b: p1, r1: 0.5, r2: 0.4 },
    { a: p1, b: p2, r1: 0.4, r2: 0.33 },
    { a: p2, b: p3, r1: 0.33, r2: 0.2 },
  ];
}
// Right atrial appendage: BROAD-BASED and triangular, wide-mouthed — the key
// discriminator from the narrow-necked LAA — lying anterior over the aortic root.
function raaSegs(ra) {
  const o = add(ra.c, add(mul(ra.e2, ra.r2 * 0.6), mul(ra.a, ra.rl * 0.35)));
  // (it lies to the RIGHT of the aortic root, its tip over the right AV groove)
  return [{ a: o, b: add(o, [-0.95, -0.3, 1.3]), r1: 0.95, r2: 0.36 }];
}

// Pulmonary veins draining into the posterior LA, left and right, superior and
// inferior. Ostia are seated on the live atrial wall (direction `u` from the LA
// centre); each vein then runs toward its lung (`d`, from the body axes). The
// right superior vein enters beside the septum, which is why it is the one seen
// on an A4C for pulmonary-vein Doppler. Ostial calibre ~1.0-1.2 cm.
const PV_DEF = [
  { name: 'LSPV', u: [0.72, 0.42, -0.55], d: unit(add(add(mul(BL0, 0.7), mul(BS0, 0.4)), mul(BP0, 0.45))), r: 0.58, len: 2.0 },
  { name: 'LIPV', u: [0.7, -0.15, -0.7], d: unit(add(add(mul(BL0, 0.7), mul(BS0, -0.35)), mul(BP0, 0.45))), r: 0.52, len: 1.9 },
  { name: 'RSPV', u: [-0.55, 0.55, -0.62], d: unit(add(add(mul(BL0, -0.7), mul(BS0, 0.4)), mul(BP0, 0.45))), r: 0.6, len: 2.0 },
  { name: 'RIPV', u: [-0.55, -0.05, -0.83], d: unit(add(add(mul(BL0, -0.7), mul(BS0, -0.35)), mul(BP0, 0.45))), r: 0.55, len: 1.9 },
].map((v) => ({ ...v, u: unit(v.u) }));
function laPoint(la, u, f) {
  // point on the LA ellipsoid surface in (heart-space) direction u, scaled by f
  const lu = dot(u, la.e1) / la.r1, la_ = dot(u, la.a) / la.rl, lw = dot(u, la.e2) / la.r2;
  const t = 1 / Math.hypot(lu, la_, lw);
  return mad(la.c, u, t * f);
}
function buildPV(la) {
  const veins = PV_DEF.map((v) => {
    const a = laPoint(la, v.u, 0.82);
    const b = mad(mad(a, v.u, 0.6), v.d, v.len);
    const ax = unit(sub(a, b));                               // inflow axis (vein -> atrium)
    return { name: v.name, r: v.r, a, b, ax };
  });
  // the muscular ridge between the LSPV ostium and the LAA ostium — the
  // "warfarin" (coumadin) ridge, a normal structure regularly mistaken for a mass
  // It runs ALONG the lateral wall between the two ostia, so its ends are taken
  // on the atrial surface (directions interpolated between the ostia) rather
  // than on the straight chord, which would cross the cavity.
  const uL = unit(sub(veins[0].a, la.c)), uA = unit(sub(laaSegs(la)[0].a, la.c));
  const ridge = { a: laPoint(la, unit(lerp3(uL, uA, 0.3)), 0.97), b: laPoint(la, unit(lerp3(uL, uA, 0.6)), 0.97), r: 0.14 };
  return { veins, ridge };
}

// Crista terminalis: the C-shaped muscular ridge from the SVC orifice around the
// lateral RA wall to the IVC, dividing the smooth sinus venarum from the
// trabeculated pectinate atrium — a classic pseudo-mass on echo. Waypoints are
// unit directions in the RA frame (a = long axis, e1 lateral, e2 anterior). Each
// is seated on the live lumen surface (cavae, sinus venarum and septal stretch
// included), so the ridge is attached to the wall at every phase and stands
// CRISTA.proud into the cavity; never a free intracavitary bar.
const CRISTA_U = [[0.15, 0.95, 0.1], [-0.55, 0.6, 0.15], [-0.85, 0.0, 0.1], [-0.6, -0.6, 0.05], [0.0, -0.95, -0.1]];
const CRISTA = { r: 0.25, proud: 0.4, sub: 3, shell: 0.45 };   // capsule radius, ridge height above the wall (cm), points per span, deepest paint
function raFramePoint(ra, u, f) {
  const d = add(add(mul(ra.e1, u[0] * ra.r1), mul(ra.a, u[1] * ra.rl)), mul(ra.e2, u[2] * ra.r2));
  return mad(ra.c, d, f);
}
function buildCrista(A) {
  const ra = A.ra, c = ra.c;
  // the wall point along the ray from the RA centre through u (bisection on the lumen)
  const wall = (u) => {
    const dir = unit(sub(raFramePoint(ra, u, 1), c));
    let lo = 0, hi = 0.1;                   // march out to the first exit, then bisect
    for (; hi < 4.0; lo = hi, hi += 0.1) {
      const q = mad(c, dir, hi);
      if (dRAlumen(q[0], q[1], q[2], A) >= 0) break;
    }
    for (let i = 0; i < 12; i++) {
      const m = 0.5 * (lo + hi), q = mad(c, dir, m);
      if (dRAlumen(q[0], q[1], q[2], A) < 0) lo = m; else hi = m;
    }
    return mad(c, dir, lo - (CRISTA.proud - CRISTA.r));   // the capsule's inner face stands `proud` off the wall
  };
  const pts = [];
  for (let i = 0; i < CRISTA_U.length - 1; i++) {
    for (let k = 0; k < CRISTA.sub; k++) pts.push(wall(lerp3(CRISTA_U[i], CRISTA_U[i + 1], k / CRISTA.sub)));
  }
  pts.push(wall(CRISTA_U[CRISTA_U.length - 1]));
  const segs = [];
  for (let i = 0; i + 1 < pts.length; i++) segs.push({ a: pts[i], b: pts[i + 1], r: CRISTA.r });
  return segs;
}
// Eustachian valve: the crescentic flap at the IVC-RA junction directed toward
// the fossa ovalis (in fetal life it streams IVC blood across the foramen ovale).
function buildEustachian(ra) {
  return { a: raFramePoint(ra, [-0.2, -0.85, -0.35], 0.9), b: raFramePoint(ra, [0.55, -0.75, -0.3], 0.82), r: 0.09 };
}

// Coronary sinus: the main cardiac vein, running in the POSTERIOR atrioventricular
// groove (between the LA and the LV posterior wall) and opening into the RA just
// above the septal tricuspid leaflet, beside the IVC. d ~0.8 cm; on PLAX it is the
// small circle at the posterior AV groove, and it DILATES with a persistent left
// SVC or raised right-sided pressure.
function buildCS(M, raFloor) {
  const pts = [];
  const R = 3.05, y0 = M[1] + 0.15;
  const angs = [-0.15, -0.75, -1.35, -1.95, -2.45];
  for (const t of angs) pts.push([M[0] + Math.cos(t) * R, y0 + 0.08 * t, M[2] + Math.sin(t) * (R - 0.2)]);
  pts.push(add(raFloor, [0.55, 0.25, -1.05]));             // ostium, low on the posterior septal RA
  const segs = [];
  for (let i = 0; i < pts.length - 1; i++) segs.push({ a: pts[i], b: pts[i + 1], r: i < 2 ? 0.34 : 0.42 });
  return segs;
}

// The aorta as one vessel: root (three sinuses of Valsalva — so the root is
// trilobed, not a tube), sino-tubular junction, ascending aorta, arch and the
// descending thoracic aorta.
function buildAorta(A, path) {
  const cuspDirs = [0, 1, 2].map((i) => {
    const a = (90 + i * 120) * Math.PI / 180;               // RCC anterior, LCC left-posterior, NCC right-posterior
    return add(mul(E_SCR, Math.cos(a)), mul(E_ANT, Math.sin(a)));
  });
  const sinusH = 0.85, sinusOff = 0.62, sinusR = AO.sinusR - 0.62; // sinus spheres -> root d 3.3
  const dil = path.dilatedAorta ? 1.25 : 1;
  const stj = mad(A, U_AO, 2.1);
  return {
    valve: A,
    sinus: cuspDirs.map((d) => mad(mad(A, U_AO, sinusH), d, sinusOff)),
    sinusR: sinusR * dil,
    root: { a: mad(A, U_AO, -0.05), b: stj, r1: AO.annR, r2: AO.stjR * dil },
    asc: { a: stj, b: AO_ASC_TOP, r1: AO.stjR * dil, r2: AO.ascR * dil },
    arch: ARCH_SEGS,
    dta: { a: DTA_TOP, b: DTA_BOT, r: DTA_R },
    wall: 0.2,
  };
}

// ---- coronary arteries ------------------------------------------------------
// Live endocardial radius of the (axisymmetric) LV body at height y, and the
// epicardial radius including the wall — used to lay vessels on the surface.
function lvRadiusAt(y, L) {
  const s = (y - L.apexY) / L.sL;
  if (s <= 0) return 0;
  let r;
  if (s <= LVP.sEq) { const q = (s - LVP.sEq) / LVP.sEq; r = LVP.b * Math.sqrt(Math.max(0, 1 - q * q)); }
  else { const u = Math.min(1, (s - LVP.sEq) / (LVP.sTop - LVP.sEq)); r = LVP.b - (LVP.b - LVP.rTop) * Math.pow(u, 2.5); }
  const grow = L.lvScale === 1 ? 1 : 1 + (L.lvScale - 1) * (1 - smoothstep(LVP.sEq - 0.5, LVP.sTop, s));
  return r * L.sS * grow;
}
function lvEpiRadiusAt(y, L) {
  const s = (y - L.apexY) / L.sL;
  return lvRadiusAt(y, L) + L.wall * (LVP.apexWallFrac + (1 - LVP.apexWallFrac) * smoothstep(0, 3.5, s));
}
// The coronary arteries run on the epicardial surface in fat-filled grooves: the
// LEFT MAIN leaves the left-coronary sinus and passes behind the PA trunk to
// divide into the LAD (anterior interventricular groove, to the apex) and the
// CIRCUMFLEX (left AV groove, under the LAA); the RCA leaves the right-coronary
// sinus (the anterior cusp) and runs round the right AV groove to the crux, where
// the PDA descends the inferior interventricular groove. On echo the ostia are
// the landmark of the PSAX-AV view (LM at ~4 o'clock off the left cusp, RCA at
// ~10-11 o'clock off the right cusp). Calibres: LM ~4.5 mm, LAD / RCA ~3.5 mm.
function buildCoronaries(L, ao, T, tvN) {
  const segs = [];
  // Each vessel is a Catmull-Rom curve through its control points, subdivided so
  // it follows the curved grooves smoothly instead of cutting chords across them.
  const SUB = 3;
  const cr = (p0, p1, p2, p3, t) => {
    const t2 = t * t, t3 = t2 * t;
    return [0, 1, 2].map((k) => 0.5 * (2 * p1[k] + (p2[k] - p0[k]) * t +
      (2 * p0[k] - 5 * p1[k] + 4 * p2[k] - p3[k]) * t2 + (3 * p1[k] - p0[k] - 3 * p2[k] + p3[k]) * t3));
  };
  const chain = (ctl, r0, r1) => {
    const pts = [ctl[0]];
    for (let i = 0; i < ctl.length - 1; i++) {
      const p0 = ctl[Math.max(0, i - 1)], p3 = ctl[Math.min(ctl.length - 1, i + 2)];
      for (let j = 1; j <= SUB; j++) pts.push(cr(p0, ctl[i], ctl[i + 1], p3, j / SUB));
    }
    for (let i = 0; i < pts.length - 1; i++) {
      const f0 = i / (pts.length - 1), f1 = (i + 1) / (pts.length - 1);
      const a = pts[i], b = pts[i + 1], ra = r0 + (r1 - r0) * f0, rb = r0 + (r1 - r0) * f1;
      const m = [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2];
      segs.push({ a, b, r1: ra, r2: rb, m, R: Math.hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]) / 2 + Math.max(ra, rb) });
    }
  };
  const onLV = (angDeg, y, lift) => {
    const a = angDeg * Math.PI / 180, r = lvEpiRadiusAt(y, L) + lift;
    return [Math.cos(a) * r, y, Math.sin(a) * r];
  };
  const axis = mad(ao.valve, U_AO, 0.85);
  const dirAt = (deg) => { const a = deg * Math.PI / 180; return add(mul(E_SCR, Math.cos(a)), mul(E_ANT, Math.sin(a))); };
  // left main from the left-coronary sinus (330 deg in the PSAX-AV frame)
  const lmO = mad(axis, dirAt(330), 0.62 + ao.sinusR - 0.12);
  const bif = mad(lmO, unit(add(add(mul(E_SCR, 0.75), mul(E_ANT, 0.45)), mul(U_AO, -0.25))), 1.1);
  chain([lmO, bif], 0.23, 0.22);
  const yb = L.base;
  // LAD: anterior interventricular groove (RV anterior insertion, ~75 deg) to the apex
  const lad = [bif, onLV(72, yb - 0.9, 0.25)];
  for (const f of [0.22, 0.45, 0.68, 0.86]) lad.push(onLV(74 + 8 * f, yb + (L.apexY - yb) * f, 0.24));
  lad.push(onLV(88, L.apexY + 0.5, 0.2));
  chain(lad, 0.2, 0.11);
  // circumflex: left AV groove, lateral then posterolateral
  const cx = [bif];
  for (const ang of [50, 25, 0, -25, -50]) cx.push(onLV(ang, yb - 0.25, 0.28));
  chain(cx, 0.18, 0.12);
  // RCA from the right-coronary sinus (anterior, 90 deg), round the right AV groove
  const rcaO = mad(axis, dirAt(90), 0.62 + ao.sinusR - 0.12);
  const ea = unit(sub([0, 0, 1], mul(tvN, dot([0, 0, 1], tvN))));          // anterior, in the TV plane
  const el = unit(sub([-1, 0, 0], mul(tvN, dot([-1, 0, 0], tvN))));        // lateral (right)
  const rg = TV_R + 0.75;
  const groove = (deg) => { const a = deg * Math.PI / 180; return mad(mad(mad(T, ea, Math.cos(a) * rg), el, Math.sin(a) * rg), tvN, 0.2); };
  const rca = [rcaO, mad(rcaO, unit(add(mul(E_ANT, 0.55), mul(E_SCR, -0.8))), 1.0), groove(15), groove(60), groove(105), groove(150)];
  const crux = onLV(-128, yb - 0.9, 0.3);
  rca.push(crux);
  chain(rca, 0.19, 0.15);
  // PDA down the inferior interventricular groove
  const pda = [crux];
  for (const f of [0.25, 0.5, 0.7]) pda.push(onLV(-126, yb + (L.apexY - yb) * f, 0.24));
  chain(pda, 0.13, 0.09);
  return segs;
}
// Distance to the nearest coronary; segments farther than `reach` are skipped
// (bounding sphere), so the result is exact only below `reach`.
function dCoronary(x, y, z, A, inflate, reach) {
  let d = 1e9;
  const c = A.cor;
  for (let i = 0; i < c.length; i++) {
    const s = c[i];
    const mx = x - s.m[0], my = y - s.m[1], mz = z - s.m[2], rr = s.R + reach;
    if (mx * mx + my * my + mz * mz > rr * rr) continue;
    const e = sdCapsule(x, y, z, s.a[0], s.a[1], s.a[2], s.b[0], s.b[1], s.b[2], s.r1 + inflate, s.r2 + inflate);
    if (e < d) d = e;
  }
  return d;
}

// ---- regional wall-motion abnormality (RWMA) -----------------------------
// A per-segment contractility defect (ischemia / infarct). Each named segment
// maps to an angular sector about the LV long axis (short-axis x–z plane), or to
// the apical cap. `severity` 0..1: 0 = normal, 1 = akinetic (the wall keeps its
// end-diastolic position and does not thicken). path.rwma may be `true`, a region
// name string, or { region, severity }.
//
// Territories follow the ASE 17-segment coronary map. The LAD regions (septal,
// anteroseptal, anterior) are longitudinal: at the base they cover only the
// anteroseptal wall (the basal inferoseptum is RCA/PDA and keeps thickening),
// and they widen toward the apex, where the LAD wraps the tip (apical cap).
// 'posterior' is the circumflex (inferolateral, ~4-5 o'clock) territory.
const RWMA_REGIONS = {
  septal:        { ang: Math.PI,        half: 1.15, lad: { ang: Math.PI * 0.8,  half: 0.75 } },
  anteroseptal:  { ang: Math.PI * 0.75, half: 1.20, lad: { ang: Math.PI * 0.75, half: 0.8 } },
  anterior:      { ang: Math.PI / 2,    half: 1.15, lad: { ang: Math.PI * 0.55, half: 0.8 } },
  lateral:       { ang: 0,              half: 1.15 },
  inferior:      { ang: -Math.PI / 2,   half: 1.15 },
  inferolateral: { ang: -Math.PI / 4,   half: 0.9 },
  apical:        { apical: true },
};
RWMA_REGIONS.posterior = RWMA_REGIONS.inferolateral;
export function buildRwma(path) {
  if (!path.rwma) return null;
  const spec = typeof path.rwma === 'object' ? path.rwma : {};
  const name = spec.region || (typeof path.rwma === 'string' ? path.rwma : 'septal');
  const base = RWMA_REGIONS[name] || RWMA_REGIONS.septal;
  const sev = spec.severity != null ? Math.max(0, Math.min(1, spec.severity)) : 0.85;
  return Object.assign({ sev, region: name }, base);
}
// Affected fraction (0..1, before severity) of the wall at short-axis angle `ang`.
// `f` is the longitudinal position (0 base .. 1 apex) and `cap` the apical-cap
// weight (0..1); both only matter for the LAD regions.
export function rwmaWeight(rw, ang, f, cap) {
  let cAng = rw.ang, half = rw.half;
  if (rw.lad) { cAng = rw.lad.ang + (rw.ang - rw.lad.ang) * f; half = rw.lad.half + (rw.half - rw.lad.half) * f; }
  let d = ang - cAng;
  while (d > Math.PI) d -= 2 * Math.PI;
  while (d < -Math.PI) d += 2 * Math.PI;
  const t = d / half;
  const w = Math.exp(-t * t);
  return rw.lad ? Math.max(w, 0.9 * cap) : w;
}
// Blend weight g in [0,1]: how "diastolic" (non-contracting) this point's wall is.
export function rwmaBlend(x, y, z, A) {
  const rw = A.lv.rwma;
  if (!rw) return 0;
  if (rw.apical) return rw.sev * Math.max(0, Math.min(1, (A.lv.apexY + 4.5 - y) / 3.0));
  const f = 1 - smoothstep(A.lv.apexY + 2.0, A.lv.apexY + 6.0, y);
  const cap = Math.max(0, Math.min(1, (A.lv.apexY + 3.6 - y) / 1.8));   // the LAD-wrapped tip only
  return rw.sev * rwmaWeight(rw, Math.atan2(z, x), f, cap);
}

// ---- LV ------------------------------------------------------------------
// Axisymmetric ED profile in (rho, s): rho = radius from the long axis, s =
// height above the apex. Lower part is a half prolate ellipsoid (semi-axes b,
// sEq); above the equator the radius holds then closes in to rTop with a rounded
// shoulder. Returns an (approximate) signed distance in ED units.
function lvProfile(rho, s) {
  const P = LVP;
  if (s <= P.sEq) return sdEllipsoid(rho, s, 0, 0, P.sEq, 0, P.b, P.sEq, P.b);
  const span = P.sTop - P.sEq;
  const u = Math.min(1, (s - P.sEq) / span);
  const r = P.b - (P.b - P.rTop) * Math.pow(u, 2.5);
  const dr = u < 1 ? (P.b - P.rTop) * 2.5 * Math.pow(u, 1.5) / span : 0;
  const side = (rho - r) / Math.sqrt(1 + dr * dr);
  return smax(side, s - P.sTop, 0.7);
}
// Trabeculation: the apical third of the LV carries fine trabeculae and the RV
// coarse ones; on echo they break up the endocardial border near the apex (and
// are why apical thrombus / non-compaction are judged there). A cheap bounded
// ridge pattern perturbs the cavity surface only close to it.
function trabecular(x, y, z, amp) {
  const a = Math.sin(x * 5.3 + y * 1.7) * Math.sin(z * 4.9 - y * 2.3) + 0.5 * Math.sin(y * 6.1 + x * 2.9 + z * 3.3);
  return amp * (a > 0.35 ? (a - 0.35) : 0);
}
// LV body cavity only (no funnels), world space.
function lvBody(x, y, z, A, g) {
  const L = A.lv;
  let sL = L.sL;
  // inward apical thickening: the endocardial apex rises by apexLift while the
  // base stays put (an akinetic segment does not thicken)
  const dA = L.apexLift * (1 - g);
  const s = (y - L.apexY - dA) / (sL - dA / LVP.len);
  // remodelling scale is full from the equator down and fades to 1 at the base,
  // which the fibrous skeleton holds (the annuli do not balloon with the cavity)
  const grow = L.lvScale === 1 ? 1 : 1 + (L.lvScale - 1) * (1 - smoothstep(LVP.sEq - 0.5, LVP.sTop, s));
  let sR = L.sS * grow;
  if (g > 0) { sR += (grow - sR) * g; sL += (L.sLDia - sL) * g; }
  const rho = Math.sqrt(x * x + z * z) / sR;
  let d = lvProfile(rho, s) * (sR < sL ? sR : sL);
  if (L.septFlat) {                                        // D-shaped LV (RV overload)
    const xp = -LVP.b * sR * (1 - L.septFlat);
    d = Math.max(d, xp - x);
  }
  return d;
}
// Full LV lumen: body + mitral inflow funnel + LV outflow tract.
function dLVlumen(x, y, z, A) {
  const g = A.lv.rwma ? rwmaBlend(x, y, z, A) : 0;
  let d = lvBody(x, y, z, A, g);
  const L = A.lv, M = L.M, Av = L.A, lo = L.lvot;
  // quick reject: the funnels live in the top 3 cm of the LV
  if (y > L.apexY + LVP.len * L.sL - 3.4) {
    const ym = y - mitralLift(x, y, z, A);                  // angle-dependent annular level
    d = smin(d, sdFrustum(x, ym, z, M[0] * L.sS, M[1] - 1.3, (M[2] + 0.2) * L.sS, M[0], M[1] + 0.03, M[2], 1.8 * L.sS, L.mvR), 0.4);
    d = smin(d, sdFrustum(x, y, z, lo[0], lo[1], lo[2], Av[0] + U_AO[0] * 0.03, Av[1] + U_AO[1] * 0.03, Av[2] + U_AO[2] * 0.03, L.lvotR, AO.annR * 0.97), 0.35);
  }
  // apical trabeculation roughens the endocardium in the apical third
  if (A.trab && y < L.apexY + 3.2 && d > -0.3 && d < 0.3) d += trabecular(x, y, z, 0.16);
  return d;
}
function lvWallAt(y, A, g) {
  const L = A.lv;
  const w = L.wall + (L.wallDia - L.wall) * g;
  const dA = L.apexLift * (1 - g);                          // same apical frame as lvBody
  const s = (y - L.apexY - dA) / (L.sL - dA / LVP.len);
  return w * (LVP.apexWallFrac + (1 - LVP.apexWallFrac) * smoothstep(0.0, 3.5, s));
}
function lvEpi(x, y, z, A) {
  const g = A.lv.rwma ? rwmaBlend(x, y, z, A) : 0;
  const L = A.lv, M = L.M, Av = L.A, lo = L.lvot;
  const w = lvWallAt(y, A, g);
  let d = lvBody(x, y, z, A, g) - w;
  if (y > L.apexY + LVP.len * L.sL - 3.6) {
    const wb = L.wall * 0.8;
    const ym = y - mitralLift(x, y, z, A);
    d = smin(d, sdFrustum(x, ym, z, M[0] * L.sS, M[1] - 1.3, (M[2] + 0.2) * L.sS, M[0], M[1] + 0.12, M[2], 1.8 * L.sS + wb, L.mvR + 0.3), 0.5);
    d = smin(d, sdFrustum(x, y, z, lo[0], lo[1], lo[2], Av[0], Av[1], Av[2], L.lvotR + wb, AO.annR + 0.3), 0.5);
  }
  return d;
}

// ---- RV ------------------------------------------------------------------
// A crescent: a body ellipsoid offset toward the right-anterior (10 o'clock in
// PSAX) CARVED by the LV epicardium, so the septum is the LV's own wall and the
// RV wraps it from the anterior to the inferior interventricular groove. An
// inflow lobe sits under the tricuspid annulus; the infundibulum (RVOT) arcs
// over the front of the aortic root to the pulmonary valve. The apex tapers
// and stops short of the LV apex (the LV forms the cardiac apex).
function rvBodyD(x, y, z, A, inflate) {
  const c = A.rv.c, r = A.rv.r;
  const u = (y - c[1]) / r[1];
  const s = u < 0 ? Math.min(1, -u / 0.9) : 0;
  const t = 1 - 0.3 * s * s;
  let d = sdEllipsoid(x, y, z, c[0], c[1], c[2], r[0] * t + inflate, r[1] + inflate, r[2] * t + inflate);
  const f = A.inflow;
  d = smin(d, sdEllipsoid(x, y, z, f.c[0], f.c[1], f.c[2], f.r[0] + inflate, f.r[1] + inflate, f.r[2] + inflate), 0.9);
  const T = A.valves.tricuspid.c, n = A.valves.tricuspid.n;
  d = smin(d, sdFrustum(x, y, z, T[0] - n[0] * 1.1, T[1] - n[1] * 1.1, T[2] - n[2] * 1.1,
    T[0] + n[0] * 0.03, T[1] + n[1] * 0.03, T[2] + n[2] * 0.03, 2.05 + inflate, TV_R + inflate), 0.4);
  // the RV sinus ends at the tricuspid annular plane (the AV junction); only the
  // infundibulum below rises above it, to the pulmonary valve
  d = smax(d, (x - T[0]) * n[0] + (y - T[1]) * n[1] + (z - T[2]) * n[2] - inflate, 0.35);
  const o = A.rvot, oR = A.rvotR;
  d = smin(d, sdCapsule(x, y, z, o[0][0], o[0][1], o[0][2], o[1][0], o[1][1], o[1][2], oR[0] + inflate, oR[1] + inflate), 0.8);
  d = smin(d, sdCapsule(x, y, z, o[1][0], o[1][1], o[1][2], o[2][0], o[2][1], o[2][2], oR[1] + inflate, oR[2] + inflate), 0.5);
  return d;
}
function dRVlumen(x, y, z, A) {
  let d = rvBodyD(x, y, z, A, 0);
  if (d > 1.2) return d;
  // the RV ends at the tricuspid annular plane: above it (within the annulus
  // footprint) is right atrium
  const T = A.valves.tricuspid.c, n = A.valves.tricuspid.n;
  const tx = x - T[0], ty = y - T[1], tz = z - T[2];
  const h = tx * n[0] + ty * n[1] + tz * n[2];
  if (h > -0.3) {
    const rr = Math.sqrt(Math.max(0, tx * tx + ty * ty + tz * tz - h * h));
    if (rr < TV_R + 0.5) d = smax(d, h, 0.1);
  }
  d = ssub(d, lvEpi(x, y, z, A), A.rvCarve);
  // the infundibulum is separated from the aortic root by the aortic wall and
  // the muscular outlet (infundibular) septum, ~0.5 cm in all
  d = ssub(d, dAOroot(x, y, z, A) - A.ao.wall - 0.35, 0.2);
  // coarse RV trabeculation toward the apex
  if (A.trab && d > -0.3 && d < 0.3 && y < A.rv.c[1] - 1.0) d += trabecular(x * 0.7, y * 0.7, z * 0.7, 0.2);
  return d;
}

// ---- atria -----------------------------------------------------------------
function sdAtrium(x, y, z, E, inflate) {
  if (E.g && E.g < 0.999) {
    // roof descent: points above the annulus map back to (height / g) along the
    // long axis, so the annulus stays fixed while the roof comes down
    const a = E.a, h = (x - E.fl[0]) * a[0] + (y - E.fl[1]) * a[1] + (z - E.fl[2]) * a[2];
    if (h > 0) {
      const k = h * (1 / E.g - 1);
      const d = sdAtriumFree(x + a[0] * k, y + a[1] * k, z + a[2] * k, E, 0);
      return inflate > 0 ? d * E.g - inflate : d;
    }
  }
  return sdAtriumFree(x, y, z, E, inflate);
}
function sdAtriumFree(x, y, z, E, inflate) {
  if (E.f && E.f < 0.999 && E.sn) {
    // free-wall emptying: points on the chamber's own side of the septal plane map
    // back to (distance / f) from the plane, so the septal face stays fixed. The
    // compression applies to the lumen shape only; wall thickness (inflate) is
    // added afterwards so the walls keep their true thickness.
    const u = E.sn, sOwn = -((x - IAS_P[0]) * u[0] + (y - IAS_P[1]) * u[1] + (z - IAS_P[2]) * u[2]);
    if (sOwn > 0) {
      const k = sOwn * (1 / E.f - 1);
      const d = sdAtriumRaw(x - u[0] * k, y - u[1] * k, z - u[2] * k, E, 0);
      return inflate > 0 ? d * E.f - inflate : d;
    }
  }
  return sdAtriumRaw(x, y, z, E, inflate);
}
function sdAtriumRaw(x, y, z, E, inflate) {
  let dx = x - E.c[0], dy = y - E.c[1], dz = z - E.c[2];
  if (E.sw > 0) {                          // stretched toward the septum (anchorToSeptum)
    const n = E.sn;
    let t = dx * n[0] + dy * n[1] + dz * n[2];
    t = t < 0 ? 0 : t > E.sw ? E.sw : t;
    dx -= n[0] * t; dy -= n[1] * t; dz -= n[2] * t;
  }
  const u = dx * E.e1[0] + dy * E.e1[1] + dz * E.e1[2];
  const v = dx * E.a[0] + dy * E.a[1] + dz * E.a[2];
  const w = dx * E.e2[0] + dy * E.e2[1] + dz * E.e2[2];
  return sdEllipsoid(u, v, w, 0, 0, 0, E.r1 + inflate, E.rl + inflate, E.r2 + inflate);
}
// signed distance to the interatrial septal plane (+ toward the RA) and the
// local half-thickness (thin fossa ovalis, thicker limbus; negative = ASD hole).
const _ias = { s: 0, t: 0 };
function iasAt(x, y, z, A) {
  const s = A.ias;
  const dx = x - s.p[0], dy = y - s.p[1], dz = z - s.p[2];
  const sd = dx * s.n[0] + dy * s.n[1] + dz * s.n[2];
  const fx = x - s.fc[0], fy = y - s.fc[1], fz = z - s.fc[2];
  const fn = fx * s.n[0] + fy * s.n[1] + fz * s.n[2];
  const rad = Math.sqrt(Math.max(0, fx * fx + fy * fy + fz * fz - fn * fn));
  const f = rad <= s.rFossa ? 0 : Math.min(1, (rad - s.rFossa) / 0.8);
  _ias.s = sd;
  _ias.t = s.tFossa + (s.tLimbus - s.tFossa) * f * f * (3 - 2 * f);
  return _ias;
}
function dAppendage(x, y, z, segs, k) {
  let d = 1e9;
  for (let i = 0; i < segs.length; i++) {
    const s = segs[i];
    d = smin(d, sdRoundCone(x, y, z, s.a[0], s.a[1], s.a[2], s.b[0], s.b[1], s.b[2], s.r1, s.r2), k);
  }
  return d;
}
function dPV(x, y, z, A) {
  let d = 1e9;
  const v = A.pv.veins;
  for (let i = 0; i < v.length; i++) {
    const s = v[i];
    d = smin(d, sdCapsule(x, y, z, s.a[0], s.a[1], s.a[2], s.b[0], s.b[1], s.b[2], s.r), 0.3);
  }
  return d;
}
function dLAlumen(x, y, z, A) {
  const e = sdAtrium(x, y, z, A.la, 0);
  if (e > 3.2) {                                                     // far away: only the vein tubes matter
    const v = Math.min(e, dPV(x, y, z, A));
    return v > 1.5 ? v : ssub(v, dAOarchDTA(x, y, z, A) - A.ao.wall - 0.15, 0.2);
  }
  let d = smin(e, dAppendage(x, y, z, A.la.aa, 0.22), 0.14);       // NARROW neck
  d = smin(d, dPV(x, y, z, A), 0.4);                                 // veins flare in
  const M = A.lv.M;                                                  // mitral orifice
  // the LA opens widely onto the mitral orifice: the inflow runs from the annulus
  // all the way into the LA body, so no partition separates them at any phase
  const inl = Math.max(0.8, Math.hypot(A.la.c[0] - M[0], A.la.c[1] - M[1], A.la.c[2] - M[2]));
  const ym = y - mitralLift(x, y, z, A);                             // (saddle deepens in systole)
  d = smin(d, sdCylinder(x, ym, z, M[0], M[1] - 0.03, M[2], M[0] + A.la.a[0] * inl, M[1] + A.la.a[1] * inl, M[2] + A.la.a[2] * inl, A.lv.mvR * 0.98), 0.6);
  d = clipAtAnnulus(d, x, ym, z, M, MV_N, MV_R + 0.6);               // LA ends at the mitral annulus
  const s = iasAt(x, y, z, A);
  d = smax(d, s.s + s.t, 0.12);                                      // cut by the septum
  // the aortic root sits in front of the LA across the transverse sinus: the
  // aortic wall and the LA wall lie between the two blood pools
  d = ssub(d, dAOroot(x, y, z, A) - A.ao.wall - 0.22, 0.15);
  // an enlarged LA abuts but never overruns the descending aorta behind it
  d = ssub(d, dAOarchDTA(x, y, z, A) - A.ao.wall - 0.15, 0.2);
  return d;
}
function dRAlumen(x, y, z, A) {
  const e = sdAtrium(x, y, z, A.ra, 0);
  if (e > 4.0) {                                                     // far away: only the caval tubes matter
    const v = A.svc, w = A.ivc;
    const t = Math.min(e, sdCapsule(x, y, z, v.a[0], v.a[1], v.a[2], v.b[0], v.b[1], v.b[2], v.r),
      sdCapsule(x, y, z, w.a[0], w.a[1], w.a[2], w.b[0], w.b[1], w.b[2], w.r));
    return t > 1.5 ? t : ssub(t, dAOroot(x, y, z, A) - A.ao.wall - 0.22, 0.15);   // SVC abuts, never enters, the aorta
  }
  let d = smin(e, dAppendage(x, y, z, A.ra.aa, 0.3), 0.7);          // BROAD mouth
  const v = A.svc, w = A.ivc;
  d = smin(d, sdCapsule(x, y, z, v.a[0], v.a[1], v.a[2], v.b[0], v.b[1], v.b[2], v.r), 0.9);
  d = smin(d, sdCapsule(x, y, z, w.a[0], w.a[1], w.a[2], w.b[0], w.b[1], w.b[2], w.r), 0.6);
  // sinus venarum: the smooth-walled venous back of the RA between the two caval
  // orifices, so SVC -> RA -> IVC is one continuous channel (the ME bicaval view)
  // It bulges away from the septum between the orifices, so the RA reads as a
  // chamber there and not as a tube of caval calibre (the ME bicaval view).
  d = smin(d, sdRoundCone(x, y, z, w.a[0], w.a[1], w.a[2], v.a[0], v.a[1], v.a[2], w.r * 1.1, v.r * 0.9), 0.6);
  const sv = A.ra.sv;
  d = smin(d, sdCapsule(x, y, z, sv.a[0], sv.a[1], sv.a[2], sv.b[0], sv.b[1], sv.b[2], sv.r), 0.6);
  // the SVC widens over its last cm into a flared orifice, its stream directed
  // at the tricuspid valve, and opens straight into its body (a funnel toward the RA
  // centre), whatever the phase
  const c = A.ra.c, fl = A.ra.svcFlare;
  d = smin(d, sdRoundCone(x, y, z, v.a[0], v.a[1], v.a[2], c[0], c[1], c[2], v.r * 0.85, 1.0), 0.5);
  d = smin(d, sdRoundCone(x, y, z, fl.a[0], fl.a[1], fl.a[2], fl.b[0], fl.b[1], fl.b[2], v.r, CAVA.flareR), 0.4);
  const T = A.valves.tricuspid.c, n = A.valves.tricuspid.n;          // tricuspid orifice
  d = smin(d, sdCylinder(x, y, z, T[0] - n[0] * 0.03, T[1] - n[1] * 0.03, T[2] - n[2] * 0.03, T[0] + n[0] * 0.8, T[1] + n[1] * 0.8, T[2] + n[2] * 0.8, TV_R * 0.98), 0.3);
  d = clipAtAnnulus(d, x, y, z, T, n, TV_R + 0.6);                   // RA ends at the tricuspid plane
  const s = iasAt(x, y, z, A);
  d = smax(d, -s.s + s.t, 0.12);
  // the aortic root (non-coronary sinus) bulges against the RA: its wall and the
  // RA wall lie between the two blood pools, never a shared lumen
  return ssub(d, dAOroot(x, y, z, A) - A.ao.wall - 0.22, 0.15);
}
// Vertical offset of the live mitral annulus (and everything hinged on it) from
// the level of its centre, at the point's angle round the annulus: each sector
// descends by its own share of the excursion (mvWeight), so in systole the
// aorto-mitral curtain rides high with the root and the posterolateral annulus
// low. Zero at the centre (growing linearly to the rim) and fading within ~2 cm
// above and below the annulus, so only the hinge-level primitives are bent.
export function mitralLift(x, y, z, A) {
  const L = A.lv, k = L.mvK;
  if (!k) return 0;
  const c = L.M, dx = x - c[0], dz = z - c[2];
  const r = Math.sqrt(dx * dx + dz * dz);
  if (r < 1e-6) return 0;
  const h = Math.abs(y - c[1]);
  if (h > 2.2) return 0;
  const w = mvWeight((dx * MV_AP[0] + dz * MV_AP[2]) / r);
  return k * (w - MV_W_MEAN) * Math.min(1, r / MV_R) * (1 - smoothstep(0.8, 2.2, h));
}
// An atrium ends at its AV-valve plane: below the plane (toward the ventricle),
// within the annulus footprint, is ventricle. n points from ventricle to atrium.
const MV_N = unit([0.02, 1, -0.08]);
function clipAtAnnulus(d, x, y, z, C, n, R) {
  const tx = x - C[0], ty = y - C[1], tz = z - C[2];
  const h = tx * n[0] + ty * n[1] + tz * n[2];
  if (h > 0.3) return d;
  const rr = Math.sqrt(Math.max(0, tx * tx + ty * ty + tz * tz - h * h));
  if (rr > R + 0.8) return d;
  // full cut inside the annulus footprint, fading out beyond it so the chamber
  // outline stays rounded rather than squared off at the plane
  const w = 1 - smoothstep(R - 0.4, R + 0.8, rr);
  return smax(d, -h * w - (1 - w) * 0.6, 0.15);
}
function dCrista(x, y, z, A) {
  let d = 1e9;
  const segs = A.crista;
  for (let i = 0; i < segs.length; i++) {
    const s = segs[i];
    d = smin(d, sdCapsule(x, y, z, s.a[0], s.a[1], s.a[2], s.b[0], s.b[1], s.b[2], s.r), 0.25);
  }
  return d;
}
function dEustachian(x, y, z, A) {
  const e = A.eustachian;
  return sdCapsule(x, y, z, e.a[0], e.a[1], e.a[2], e.b[0], e.b[1], e.b[2], e.r);
}
function dCS(x, y, z, A, inflate) {
  let d = 1e9;
  const c = A.cs;
  for (let i = 0; i < c.length; i++) {
    const s = c[i];
    d = smin(d, sdCapsule(x, y, z, s.a[0], s.a[1], s.a[2], s.b[0], s.b[1], s.b[2], s.r + inflate), 0.2);
  }
  return d;
}

// ---- great vessels ---------------------------------------------------------
function dAOroot(x, y, z, A) {
  const ao = A.ao;
  let d = sdRoundCone(x, y, z, ao.root.a[0], ao.root.a[1], ao.root.a[2], ao.root.b[0], ao.root.b[1], ao.root.b[2], ao.root.r1, ao.root.r2);
  for (let i = 0; i < 3; i++) { const s = ao.sinus[i]; d = smin(d, sdSphere(x, y, z, s[0], s[1], s[2], ao.sinusR), 0.3); }
  d = smin(d, sdRoundCone(x, y, z, ao.asc.a[0], ao.asc.a[1], ao.asc.a[2], ao.asc.b[0], ao.asc.b[1], ao.asc.b[2], ao.asc.r1, ao.asc.r2), 0.4);
  return d;
}
function dAOarchDTA(x, y, z, A) {
  const ao = A.ao;
  let d = 1e9;
  for (let i = 0; i < ao.arch.length; i++) {
    const s = ao.arch[i];
    d = smin(d, sdCapsule(x, y, z, s[0][0], s[0][1], s[0][2], s[1][0], s[1][1], s[1][2], s[2]), 0.6);
  }
  const t = ao.dta;
  return smin(d, sdCapsule(x, y, z, t.a[0], t.a[1], t.a[2], t.b[0], t.b[1], t.b[2], t.r), 0.6);
}
function dAOlumen(x, y, z, A) { return Math.min(dAOroot(x, y, z, A), dAOarchDTA(x, y, z, A)); }
function dPAlumen(x, y, z, A) {
  const pa = A.pa;
  let d = sdRoundCone(x, y, z, pa.main.a[0], pa.main.a[1], pa.main.a[2], pa.main.b[0], pa.main.b[1], pa.main.b[2], pa.main.r1, pa.main.r2);
  for (const seg of pa.branch) d = smin(d, sdCapsule(x, y, z, seg[0][0], seg[0][1], seg[0][2], seg[1][0], seg[1][1], seg[1][2], seg[2]), 0.4);
  // the two great arteries touch but never share a lumen (walls between them)
  if (d < 1.0) d = ssub(d, dAOlumen(x, y, z, A) - A.ao.wall - pa.wall, 0.1);
  return d;
}
function dPap(x, y, z, A) {
  let d = 1e9;
  for (const p of A.pap) d = smin(d, sdRoundCone(x, y, z, p.a[0], p.a[1], p.a[2], p.b[0], p.b[1], p.b[2], p.rb, p.rt), 0.35);
  return d;
}
function dMod(x, y, z, A) {
  const m = A.mod;
  let d = sdCapsule(x, y, z, m.a[0], m.a[1], m.a[2], m.b[0], m.b[1], m.b[2], m.r);
  const p = A.rvPap;
  if (p) d = smin(d, sdCapsule(x, y, z, p.a[0], p.a[1], p.a[2], p.b[0], p.b[1], p.b[2], p.r, p.r * 0.6), 0.35);
  return d;
}

// ---- walls -------------------------------------------------------------------
// Myocardial + atrial shell only (no great vessels). Used as the 3D muscle
// surface so the aorta/PA render as distinct vessel meshes. Computed in two
// parts so the classifier can also use the CORE (chambers and appendages,
// without the venous tubes) as the surface the pericardial sac wraps.
const _myo = { core: 0, all: 0 };
function myoParts(x, y, z, A) {
  let d = lvEpi(x, y, z, A);
  // RV free wall = RV lumen inflated (the septum is the LV's own wall)
  d = smin(d, rvBodyD(x, y, z, A, A.rvWall), 0.3);
  // atrial walls (thin), with their appendages, veins and the septum — each behind
  // a cheap bound on its ellipsoid since this runs for every sample.
  // Each atrial wall stops at the interatrial septal plane (the other atrium's
  // wall — or the septum itself — continues across it), so the part of an
  // atrial ellipsoid that the septum cuts away leaves no phantom muscle behind.
  let tubes = 1e9;
  const la = sdAtrium(x, y, z, A.la, 0.22), ra = sdAtrium(x, y, z, A.ra, 0.22);
  const sep = (la < 3.2 || ra < 4.0) ? iasAt(x, y, z, A).s : 0;
  if (la < 3.2) {
    d = smin(d, smin(smax(la, sep - 0.02, 0.1), dAppendage(x, y, z, A.la.aa, 0.22) - 0.15, 0.35), 0.4);
    tubes = dPV(x, y, z, A) - 0.1;
  } else d = Math.min(d, la);
  if (ra < 4.0) {
    d = smin(d, smin(smax(ra, -sep - 0.02, 0.1), dAppendage(x, y, z, A.ra.aa, 0.3) - 0.15, 0.4), 0.4);
    const v = A.svc, iv = A.ivc;
    tubes = Math.min(tubes, sdCapsule(x, y, z, v.a[0], v.a[1], v.a[2], v.b[0], v.b[1], v.b[2], v.r + 0.1));
    tubes = Math.min(tubes, sdCapsule(x, y, z, iv.a[0], iv.a[1], iv.a[2], iv.b[0], iv.b[1], iv.b[2], iv.r + 0.1));
    // (the sinus venarum and the caval funnels carry the RA wall with them)
    const sv = A.ra.sv, c = A.ra.c;
    let sw = sdRoundCone(x, y, z, iv.a[0], iv.a[1], iv.a[2], v.a[0], v.a[1], v.a[2], iv.r * 1.1 + 0.18, v.r * 0.9 + 0.18);
    sw = Math.min(sw, sdCapsule(x, y, z, sv.a[0], sv.a[1], sv.a[2], sv.b[0], sv.b[1], sv.b[2], sv.r + 0.18));
    sw = Math.min(sw, sdRoundCone(x, y, z, v.a[0], v.a[1], v.a[2], c[0], c[1], c[2], v.r * 0.85 + 0.18, 1.18));
    const fl = A.ra.svcFlare;
    sw = Math.min(sw, sdRoundCone(x, y, z, fl.a[0], fl.a[1], fl.a[2], fl.b[0], fl.b[1], fl.b[2], v.r + 0.18, CAVA.flareR + 0.18));
    tubes = Math.min(tubes, smax(sw, -sep - 0.02, 0.1));
  } else d = Math.min(d, ra);
  // coronary sinus wall in the posterior AV groove
  tubes = Math.min(tubes, dCS(x, y, z, A, 0.1));
  _myo.core = d;
  _myo.all = smin(d, tubes, 0.4);
  return _myo;
}
export function myoDist(x, y, z, A) { return myoParts(x, y, z, A).all; }

// Whole-heart epicardial surface including the great-vessel walls up to the
// pericardial reflection (used for the pericardium and effusion).
export function epiDist(x, y, z, A) {
  let d = myoParts(x, y, z, A).core;
  d = smin(d, dAOroot(x, y, z, A) - A.ao.wall, 0.4);
  d = smin(d, dPAlumen(x, y, z, A) - A.pa.wall, 0.4);
  return d;
}

// The visceral envelope an effusion is measured from: the epicardium plus the
// coronary sinus with its wall and a sheath of AV-groove fat. The CS lies deep to
// the visceral pericardium against the LV and LA, so fluid starts outside it and
// never opens up between the CS and the heart.
export function visceralDist(x, y, z, A) {
  return smin(epiDist(x, y, z, A), dCS(x, y, z, A, 0.35), 0.3);
}

// Distance to a named blood pool (for 3D chamber surfaces).
export function lumenDist(x, y, z, A, which) {
  switch (which) {
    case 'LV': return dLVlumen(x, y, z, A);
    case 'RV': return dRVlumen(x, y, z, A);
    case 'LA': return dLAlumen(x, y, z, A);
    case 'RA': return dRAlumen(x, y, z, A);
    case 'AO': return dAOlumen(x, y, z, A);
    case 'AOROOT': return dAOroot(x, y, z, A);
    case 'PA': return dPAlumen(x, y, z, A);
    case 'CS': return dCS(x, y, z, A, 0);
  }
  return 1e9;
}

// ---- extracardiac: diaphragm, liver, hepatic veins -------------------------
// Height along the body-inferior axis of the diaphragm below the heart, with
// its two cupolae. Positive = below the diaphragm. The central tendon is pressed
// down by the heart; the domes rise on either side of it, the right (over the
// liver) higher than the left.
function domeRise(l) { return (l < 0 ? 0.02 : 0.012) * l * l; }
function belowDiaphragm(x, y, z) {
  const h = x * DIAPH_N[0] + y * DIAPH_N[1] + z * DIAPH_N[2] - diaphH();
  const l = x * BL[0] + y * BL[1] + z * BL[2] - MIDLINE;
  return h + domeRise(l);                                       // a lateral point sits deeper under the dome
}
// Diaphragm lift. The diaphragm is one continuous sheet, so its local rise
// under the heart is a HEIGHT FIELD over the diaphragm plane: for each column,
// how far above the dome the (end-diastolic) pericardial sac bottoms out, when
// that is within DIAPH_WRAP. The field is smoothed so the sheet curves up to
// the heart and back down to its dome. Built once, on first use.
let _lift = null;
const LIFT = { n: 41, span: 12.0 };
function buildLift() {
  const A = anatomyParams(0, 0, {}, null, 0);
  const N = DIAPH_N, u1 = unit(cross(N, Math.abs(N[1]) < 0.9 ? [0, 1, 0] : [1, 0, 0])), u2 = cross(N, u1);
  const n = LIFT.n, step = 2 * LIFT.span / (n - 1);
  const c0 = mad([-1.0, -3.0, 0], N, diaphH() - dot([-1.0, -3.0, 0], N)); // heart centre on the plane
  const raw = new Float32Array(n * n);
  const sacAt = (p) => epiDist(p[0], p[1], p[2], A) - (PERI.gap + PERI.thick + 0.02);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const base = mad(mad(c0, u1, -LIFT.span + i * step), u2, -LIFT.span + j * step);
      let lift = 0;
      for (let t = 0; t <= DIAPH_WRAP; t += 0.1) {     // up from the dome plane (-N is up)
        const p = mad(base, N, -t);
        const l = dot(p, BL) - MIDLINE;
        const dome = domeRise(l);                      // the dome's own rise here
        if (t < dome) continue;
        if (sacAt(p) < 0) { lift = t - dome + 0.1; break; }    // just meets the sac
      }
      raw[i * n + j] = lift;
    }
  }
  // bridge narrow clefts (the interatrial groove, the AV grooves): a column may not
  // rise far above its neighbours, so the sheet meets the heart's broad inferior
  // surface without climbing up into the grooves between the chambers
  {
    const rr = 2, cp = raw.slice();
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
      if (!(cp[i * n + j] > 0)) continue;
      let mn = cp[i * n + j];
      for (let a = -rr; a <= rr; a++) for (let b = -rr; b <= rr; b++) {
        const ii = i + a, jj = j + b;
        if (ii < 0 || jj < 0 || ii >= n || jj >= n) continue;
        const v = cp[ii * n + jj];
        if (v > 0 && v < mn) mn = v;
      }
      raw[i * n + j] = Math.min(cp[i * n + j], mn + 3.0);
    }
  }
  // smooth: separable Gaussian (sigma ~1.2 cm), keeping full contact under the heart
  const sig = 1.2 / step, R = Math.ceil(2.5 * sig), k = [];
  for (let t = -R; t <= R; t++) k.push(Math.exp(-t * t / (2 * sig * sig)));
  const ks = k.reduce((a, b) => a + b, 0);
  const tmp = new Float32Array(n * n), sm = new Float32Array(n * n);
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    let a = 0; for (let t = -R; t <= R; t++) { const jj = Math.min(n - 1, Math.max(0, j + t)); a += raw[i * n + jj] * k[t + R]; }
    tmp[i * n + j] = a / ks;
  }
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
    let a = 0; for (let t = -R; t <= R; t++) { const ii = Math.min(n - 1, Math.max(0, i + t)); a += tmp[ii * n + j] * k[t + R]; }
    sm[i * n + j] = Math.max(raw[i * n + j], a / ks * 1.6 > raw[i * n + j] ? Math.min(a / ks * 1.6, DIAPH_WRAP) : raw[i * n + j]);
  }
  _lift = { u1, u2, c0, step, map: sm };
}
function diaphragmLift(x, y, z) {
  if (!_lift) buildLift();
  const L = _lift, n = LIFT.n;
  const dx = x - L.c0[0], dy = y - L.c0[1], dz = z - L.c0[2];
  const fi = (dx * L.u1[0] + dy * L.u1[1] + dz * L.u1[2] + LIFT.span) / L.step;
  const fj = (dx * L.u2[0] + dy * L.u2[1] + dz * L.u2[2] + LIFT.span) / L.step;
  if (fi < 0 || fj < 0 || fi >= n - 1 || fj >= n - 1) return 0;
  const i = fi | 0, j = fj | 0, a = fi - i, b = fj - j, m = L.map;
  return (m[i * n + j] * (1 - a) + m[(i + 1) * n + j] * a) * (1 - b) +
         (m[i * n + j + 1] * (1 - a) + m[(i + 1) * n + j + 1] * a) * b;
}
// Height below the (lifted) diaphragm sheet, + inferior: for audits and tools.
export function diaphragmBelow(x, y, z) { return belowDiaphragm(x, y, z) + diaphragmLift(x, y, z); }
// Stomach (fundus and upper body): signed distance from its lumen surface, in body
// left / posterior / below-diaphragm coordinates (the gastric wall is the shell
// 0 < d < GASTRIC.t). The same ellipsoid hollows the liver's gastric impression.
// h is measured below the dome sheet WITHOUT its local lift round the sac (that
// drape falls off steeply beside the heart and would bend the organ into a
// straight-sided wedge); the stomach is classified only below the lifted sheet.
function dStomach(l, p, h) {
  const S = LIV_LOBE.stomach;
  return sdEllipsoid(l, p, h, S.c[0], S.c[1], S.c[2], S.r[0], S.r[1], S.r[2]);
}
export function stomachDist(x, y, z) {
  return dStomach(x * BL[0] + y * BL[1] + z * BL[2] - MIDLINE, x * BP[0] + y * BP[1] + z * BP[2], belowDiaphragm(x, y, z));
}

// Liver: right-upper-quadrant solid organ beneath the diaphragm whose left lobe
// crosses the midline under the heart — the subcostal acoustic window.
function dLiver(x, y, z, h, hRaw) {
  if (h < -1) return -h;
  const l = x * BL[0] + y * BL[1] + z * BL[2] - MIDLINE;        // + left
  const p = x * BP[0] + y * BP[1] + z * BP[2];                   // + posterior
  // A smooth organ: its dome is the diaphragm; the rest is bounded by a large
  // ellipsoid centred right of the midline — ~20 cm across (the left lobe
  // tapering as it crosses under the heart), ~16 cm AP, ~14 cm tall.
  const ell = sdEllipsoid(l, p, h, -1.0, -1.0, 6.0, 11.5, 8.0, 7.0);
  // The left lobe is a thin anterior wedge. Left of the midline its visceral
  // face is hollowed by the stomach (the gastric impression) and its lower edge
  // rises toward the dome, so it tapers to a rounded tip under the apex instead
  // of filling the space behind the heart (seen edge-on from the parasternal
  // windows as a straight wall of liver). All faces are curved, so no imaging
  // plane cuts them in a straight line.
  const Q = LIV_LOBE, lp = l > 0 ? l : 0;
  const stom = dStomach(l, p, hRaw) - GASTRIC.t - 0.25;            // outside the wall and a gap
  const back = (p + Q.pk * l + Q.pc * lp * lp - Q.p0) / Math.hypot(1, Q.pk + 2 * Q.pc * lp);
  const low = (h + Q.hk * l + Q.hc * lp * lp - Q.h0) / Math.hypot(1, Q.hk + 2 * Q.hc * lp);
  return smax(smax(smax(smax(0.1 - h, ell, 1.5), back, 1.2), low, 1.2), -stom, 0.8);
}
// left-lobe faces (body l / p / below, cm): posterior p < p0 - pk*l - pc*l^2 and
// inferior h < h0 - hk*l - hc*l^2 (l > 0 for the quadratic terms); the stomach
// (fundus and body under the left dome) as an ellipsoid the lobe wraps around
const LIV_LOBE = {
  p0: 5.0, pk: 0.6, pc: 0.06, h0: 7.0, hk: 0.5, hc: 0.05,
  // lumen: fundus/upper body ~5 x 6 x 6.5 cm, entirely left of the midline (wall
  // included); its upper-right surface carries the TG contact (GASTRIC.p, with
  // the lumen normal there along GASTRIC.n). h: below the unlifted dome sheet
  stomach: { c: [2.93, 1.56, 2.99], r: [2.6, 3.0, 3.2] },
};
function dHepaticVeins(x, y, z, A) {
  const ivc = A.ivc;
  // three hepatic veins fanning into the IVC just below the diaphragm, their
  // confluence on the IVC axis 1 cm below the RA junction (so they open through
  // its wall ~1-2.5 cm below it). Each runs out through the liver obliquely:
  // the right vein laterally, the middle and left veins forward into the
  // anterior liver, so the middle vein lies in the subcostal IVC plane and is
  // seen joining the IVC's anterior wall there.
  const j = mad(mad(ivc.a, BODY_AX.I, HV_CONF), BL, HV_LEFT);   // (on the IVC's right-anterior wall)
  let d = 1e9;
  // the middle vein, the subcostal IVC landmark, is the widest (~1 cm at the
  // confluence) and runs forward, down and a little right (A/I ~ 0.6/0.8), in
  // the subcostal IVC plane (rolled toward the patient's right), ~4 cm in plane
  for (const [l, a, i, r1, r2, len] of HV_DEF) {
    const dir = unit(add(add(mul(BL, l), mul(BODY_AX.A, a)), mul(BODY_AX.I, i)));
    const e = mad(j, dir, len);
    d = Math.min(d, sdCapsule(x, y, z, j[0], j[1], j[2], e[0], e[1], e[2], r1, r2));
  }
  return d;
}
const HV_CONF = 1.0, HV_LEFT = -0.25;
// right / middle / left hepatic veins: direction (L, A, I), radius at the
// confluence and at the far end, length (cm)
const HV_DEF = [[-0.7, 0.1, 0.5, 0.45, 0.22, 5.5], [-0.55, 0.6, 0.8, 0.5, 0.3, 5.5], [0.55, 0.4, 0.5, 0.45, 0.22, 5.5]];

// Tissue codes must match cardiac-model TISSUE (imported there). We return a
// small code + echogenicity; cardiac-model maps codes to its TISSUE enum.
export const BODY = {
  OUT: 0, MYO: 1, LV: 2, RV: 3, LA: 4, RA: 5, AO: 6, PAP: 10,
  PERI: 11, LIVER: 12, VESSELWALL: 13, FAT: 14, VEIN: 15, LUNG: 16,
};

// The pericardial sac: the parietal layer is a bright fibrous line hugging the
// heart and the proximal great vessels, ~2 mm thick, separated from the
// epicardium by a virtual space (which an effusion opens up). Epicardial fat sits
// in the AV and interventricular grooves.
export const PERI = { gap: 0.1, thick: 0.15 };

// Classify a point as blood / myocardium / vessel / surrounding tissue (valves &
// effusion handled by the caller). Returns { code, echo }. `periOff` (>=0, cm)
// pushes the parietal pericardium outward — the effusion.
const _bc = { code: 0, echo: 0 };
function bc(code, echo) { _bc.code = code; _bc.echo = echo; return _bc; }
export function bodyClassify(x, y, z, A, periOff = 0) {
  const dlv = dLVlumen(x, y, z, A);
  const drv = dRVlumen(x, y, z, A);

  // papillary muscles & moderator band read as muscle inside the blood pool
  if (dlv < 0.2 && dPap(x, y, z, A) < 0) return bc(BODY.PAP, 0.7);
  if (drv < 0.2 && dMod(x, y, z, A) < 0) return bc(BODY.PAP, 0.66);

  const dla = dLAlumen(x, y, z, A);
  const dra = dRAlumen(x, y, z, A);
  // right-atrial ridges standing proud of the wall: the crista terminalis reads
  // as muscle (a classic pseudo-mass), the Eustachian valve as brighter fibrous
  // tissue. Gated on being in/near the RA so the tests cost nothing elsewhere.
  if (dra < 0.25) {
    // (only within the wall shell: the ridge can never read as a free bar)
    if (dra > -CRISTA.shell && dCrista(x, y, z, A) < 0) return bc(BODY.PAP, 0.68);
    if (dEustachian(x, y, z, A) < 0) return bc(BODY.PAP, 0.82);
  }
  // the LSPV/LAA ("warfarin") ridge stands into the LA — reads as bright tissue
  if (dla < 0.25) {
    const g = A.pv.ridge;
    if (sdCapsule(x, y, z, g.a[0], g.a[1], g.a[2], g.b[0], g.b[1], g.b[2], g.r) < 0) return bc(BODY.PAP, 0.78);
  }
  const dao = dAOlumen(x, y, z, A);
  const dpa = dPAlumen(x, y, z, A);
  const dcs = dCS(x, y, z, A, 0);

  // nearest blood pool that we are inside
  let best = 0.0, code = BODY.OUT;
  if (dlv < best) { best = dlv; code = BODY.LV; }
  if (drv < best) { best = drv; code = BODY.RV; }
  if (dla < best) { best = dla; code = BODY.LA; }
  if (dra < best) { best = dra; code = BODY.RA; }
  if (dao < best) { best = dao; code = BODY.AO; }
  if (dpa < best) { best = dpa; code = BODY.AO; }
  if (dcs < best) { best = dcs; code = BODY.RA; }
  if (code !== BODY.OUT) return bc(code, 0.03);

  // myocardium / vessel walls: inside the epicardial body but outside every lumen
  const m = myoParts(x, y, z, A);
  // epicardial coronary arteries: arterial lumen, bright wall, a sheath of fat
  if (m.all < 1.0 && m.all > -0.35) {
    const dc = dCoronary(x, y, z, A, 0, 0.45);
    if (dc < 0.45) {
      if (dc < 0) return bc(BODY.AO, 0.03);
      if (dc < 0.06) return bc(BODY.VESSELWALL, 0.5);
      if (m.all > -0.05) return bc(BODY.FAT, 0.5);
    }
  }
  // venous tubes (SVC, IVC, pulmonary veins, CS) have thin fibrous walls, not
  // muscle: a grazing cut through one must read as a faint line, not a solid mass
  if (m.all < 0) {
    if (m.core < 0.02) return bc(BODY.MYO, 0.55);
    // the hepatic veins open through the IVC wall (their confluence), not behind it
    if (dHepaticVeins(x, y, z, A) < 0) return bc(BODY.VEIN, 0.03);
    return bc(BODY.VESSELWALL, 0.35);
  }
  if (dao < A.ao.wall || dpa < A.pa.wall) return bc(BODY.VESSELWALL, 0.3);

  // pericardium: the parietal layer wraps the chambers, the aortic root and the
  // PA trunk (not the venous tubes, which pierce it); epicardial fat fills the
  // AV groove around the coronary sinus
  // The heart rests on the diaphragm: the fibrous pericardium is fused to its
  // central tendon, so the diaphragm (and the liver beneath it) rises to meet
  // the pericardial sac instead of leaving an echo-free wedge that would mimic
  // an effusion on subcostal views (see diaphragmLift).
  // the interatrial (Waterston's) groove between the two atria holds fat, never
  // diaphragm or liver, however deep it is at end-diastole
  if (m.core > 0 && dla + dra < 0.9) return bc(BODY.FAT, 0.45);
  const bRaw = belowDiaphragm(x, y, z), below = bRaw + diaphragmLift(x, y, z);
  // Anterior epicardial fat pad: a few mm of fat between the RV free wall and
  // the parietal pericardium, thickest anteriorly — the classic mimic of an
  // anterior effusion (it is granular, not echo-free, and does not track
  // posteriorly), so it lifts the pericardium off the heart only in front.
  let fp = 0;
  if (periOff === 0) {
    const ab = -((x - _HC[0]) * BP[0] + (y - _HC[1]) * BP[1] + (z - _HC[2]) * BP[2]);
    fp = FAT_PAD * smoothstep(1.2, 3.0, ab);
  }
  const off = PERI.gap + periOff + fp;
  if (fp > 0.02 && m.core > 0 && m.core < off - 0.02 && below < 0.1) return bc(BODY.FAT, 0.35);
  // With an effusion the coronary sinus stays tethered in its AV-groove fat: the
  // visceral envelope includes the CS + fat sheath, so the parietal layer (and the
  // fluid inside it) starts outside it and never runs between the CS and the heart
  if (periOff > 0 && dcs < 0.35) return bc(BODY.FAT, 0.62);
  const mCore = periOff > 0 ? Math.min(m.core, dcs - 0.35) : m.core;   // (m is shared scratch: keep it)
  if (below < 0.1 + periOff && mCore < off + PERI.thick + 0.3) {   // (an effusion sac bulges the tendon down)
    const dr = dAOroot(x, y, z, A) - A.ao.wall;
    const epi = Math.min(mCore, dr, dpa - A.pa.wall);
    // (not in the transverse sinus — the potential space squeezed between the
    // great arteries and the atria, where no parietal layer is seen)
    // nor in the interatrial groove between the two atria (only fat enters it)
    const squeezed = periOff === 0 && (Math.max(mCore, Math.min(dr, dpa - A.pa.wall)) < 0.7 ||
      (dla < 2.0 && dra < 2.0));
    // The layer is a fixed ~1.5 mm in SPACE. The blended epicardial field is not
    // a true distance where two parts of the heart meet (its gradient drops on
    // the ridge between them), and a band of field values there is a thick slab
    // that swells and shrinks as the chambers beat — so the offset is measured
    // along the gradient (value / |grad|) before it is tested. (An effusion's
    // sac, pushed well off the heart, is left as it is: there the field is
    // smooth and the ridges would only break the parietal line.)
    if (epi > 0 && epi < off + PERI.thick && !squeezed) {
      const en = periOff > 0 ? epi : epi / epiGrad(x, y, z, A, epi);
      if (en < off + PERI.thick && en > off - 0.02) return bc(BODY.PERI, 0.3);
    }
  }
  if (periOff === 0 && mCore < off && dcs < 0.55) return bc(BODY.FAT, 0.62);

  // below the diaphragm: the gastric fundus wall against the transgastric probe
  // (mucosa, dark muscularis, bright serosa), then liver parenchyma with its
  // hepatic veins and the IVC
  if (below > -0.4) {
    const liv = dLiver(x, y, z, below, bRaw);
    if (liv >= 0) {
      // the stomach, outside the liver and left of the midline only
      const lb = x * BL[0] + y * BL[1] + z * BL[2] - MIDLINE;
      if (lb > 0) {
        const gi = dStomach(lb, x * BP[0] + y * BP[1] + z * BP[2], bRaw);
        if (gi < 0) return bc(BODY.LUNG, 0.9);                     // swallowed gas in the lumen
        // the gut signature: echogenic mucosa, hypoechoic muscularis, echogenic serosa
        if (gi < GASTRIC.t) return gi < 0.12 || gi > 0.33 ? bc(BODY.VESSELWALL, gi < 0.12 ? 0.55 : 0.65) : bc(BODY.FAT, 0.3);
      }
    } else {
      const iv = A.ivc;
      if (sdCapsule(x, y, z, iv.a[0], iv.a[1], iv.a[2], iv.b[0], iv.b[1], iv.b[2], iv.r) < 0) return bc(BODY.VEIN, 0.03);
      if (dHepaticVeins(x, y, z, A) < 0) return bc(BODY.VEIN, 0.03);
      // Glisson capsule + diaphragm: one bright line at the liver surface,
      // fused with the parietal pericardium where the liver meets the heart.
      // Under the dome the line is the diaphragm itself (~2.5 mm); on the
      // visceral faces (gastric impression, lower edge) only the thin capsule.
      const cap = below < 0.4 ? 0.25 : below > 1.2 ? 0.1 : 0.25 - 0.15 * (below - 0.4) / 0.8;
      if (liv > -cap) return bc(BODY.PERI, 0.8);
      const fuse = PERI.gap + PERI.thick + 0.05;
      if (mCore < fuse && mCore / epiGrad(x, y, z, A, Math.min(mCore, dAOroot(x, y, z, A) - A.ao.wall, dpa - A.pa.wall)) < fuse) return bc(BODY.PERI, 0.8);
      return bc(BODY.LIVER, 0.42);
    }
  }
  // aerated lung: everywhere above the diaphragm beyond the mediastinal soft
  // tissue that invests the heart, except the posterior mediastinum (aorta,
  // oesophagus) and the two acoustic windows (the cardiac notch under the
  // parasternal window and the intercostal window over the apex)
  if (below < -0.5 && periOff === 0) {
    const sac = Math.min(mCore, dAOroot(x, y, z, A) - A.ao.wall, dpa - A.pa.wall);
    // (the lingula and left lung lie closer against the sac laterally, behind the
    // apex and the LV free wall, than in front of the heart)
    const gl = (x - _HC[0]) * BL[0] + (y - _HC[1]) * BL[1] + (z - _HC[2]) * BL[2];
    const gap = gl < 0 ? LUNG_GAP : gl > 2.5 ? LUNG_GAP_LAT : LUNG_GAP - (LUNG_GAP - LUNG_GAP_LAT) * gl / 2.5;
    if (sac > gap && isLung(x, y, z)) return bc(BODY.LUNG, 0.9);
  }
  return bc(BODY.OUT, 0);
}
// |grad| of the whole-heart epicardial field (chambers, aortic root, PA trunk)
// by forward differences, given its value e at the point; floored so a flat
// ridge reads as far from the surface rather than dividing by ~0.
function epiGrad(x, y, z, A, e) {
  const h = 0.03;
  const f = (px, py, pz) => Math.min(myoParts(px, py, pz, A).core,
    dAOroot(px, py, pz, A) - A.ao.wall, dPAlumen(px, py, pz, A) - A.pa.wall);
  const gx = f(x + h, y, z) - e, gy = f(x, y + h, z) - e, gz = f(x, y, z + h) - e;
  return Math.max(0.2, Math.sqrt(gx * gx + gy * gy + gz * gz) / h);
}
const FAT_PAD = 0.35;                                    // anterior epicardial fat pad (cm)
const AORTA_PLEURA = 0.5;                                // pleura this far behind the aortic wall (cm)
const LUNG_GAP = 1.0;                                    // mediastinal tissue around the sac (cm)
const RIGHT_LUNG_L = -3.5, RIGHT_LUNG_A = 0.0;           // right pleura: this far right of the midline, in front of this (cm)
const LUNG_GAP_LAT = 0.3;                                // ... lateral to the LV
const _HC = [-1.0, -1.5, 0];                             // mid-heart
const LUNG_WINDOWS = (() => {
  const apex = [0, LVP.apexY - LVP.wall * LVP.apexWallFrac, 0];
  const crux = lerp3(M0, T0, 0.5);
  const wa = mad(apex, unit(sub(apex, crux)), 2.5);                   // apical window (skin)
  const wp = mad([-0.35, 0.3, 0.35], BODY_AX.A, 7.0);                 // left parasternal window
  return [wa, wp].map((w) => ({ w, u: unit(sub(_HC, w)), tan: Math.tan(32 * Math.PI / 180) }));
})();
function isLung(x, y, z) {
  const lb = x * BL[0] + y * BL[1] + z * BL[2] - MIDLINE;
  const ab = -(x * BP[0] + y * BP[1] + z * BP[2]);                     // + anterior
  // the right lung covers the anterolateral wall of the RA (there is no acoustic
  // window there: the cardiac notch is on the left)
  if (lb < RIGHT_LUNG_L && ab > RIGHT_LUNG_A) return true;
  // general lung boundary (negative = aerated): everywhere outside the two
  // acoustic-window cones and, posteriorly, outside the posterior mediastinum
  let g = ab < 0 ? 2.8 - Math.hypot(lb, Math.min(0, ab + 1.5)) : -1e9;   // posterior mediastinum (rounded)
  for (const W of LUNG_WINDOWS) {
    const vx = x - W.w[0], vy = y - W.w[1], vz = z - W.w[2];
    const al = vx * W.u[0] + vy * W.u[1] + vz * W.u[2];
    if (al <= 0) continue;
    const px = vx - W.u[0] * al, py = vy - W.u[1] * al, pz = vz - W.u[2] * al;
    const c = al * W.tan + 1.2 - Math.sqrt(px * px + py * py + pz * pz);
    if (c > g) g = c;
  }
  if (ab < 0) {
    // the left pleura/lung wraps the descending aorta posterolaterally: aerated lung
    // from AORTA_PLEURA behind the aortic wall (the bright pleural line behind it
    // on the descending-aorta short axis), whatever the acoustic windows say.
    // The shell and the general boundary are one pleural surface: a smooth union,
    // so it curves round the aorta without cusps where the two meet.
    const dx = x - DTA_P[0], dy = y - DTA_P[1], dz = z - DTA_P[2], s = dx * BS[0] + dy * BS[1] + dz * BS[2];
    const ox = dx - BS[0] * s, oy = dy - BS[1] * s, oz = dz - BS[2] * s;
    const od = Math.hypot(ox, oy, oz);
    const pl = ox * (BL[0] * 0.5 + BP[0] * 0.85) + oy * (BL[1] * 0.5 + BP[1] * 0.85) + oz * (BL[2] * 0.5 + BP[2] * 0.85);
    const s1 = Math.max(DTA_R + AORTA_PLEURA - od, 0.45 * od - pl, od - 6);
    return smin(s1, g, 2.5) < 0;
  }
  return g < 0;
}
