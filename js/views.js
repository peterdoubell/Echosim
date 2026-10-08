// views.js — the standard echocardiographic windows, derived from the anatomy.
//
// Every view is built from anatomical landmarks (anatomy.js LM) rather than
// tuned by hand, and each one fixes the SCREEN-RIGHT direction explicitly so the
// image is displayed in the ASE/EACVI orientation a sonographer expects:
//   PLAX    aorta / LA on the right, LV apex to the left, RV(OT) in the near field
//   PSAX    viewed from the apex: RV on the left and anterior, LV lateral wall right
//   PSAX-AV aortic valve central, RVOT arching over it to the pulmonary valve (right)
//   A4C     apex at the top, LV + LA on the right, RV + RA on the left
//   A2C     anterior wall on the right, inferior wall on the left
//   SUBCOSTAL liver in the near field, RV nearest, apex pointing to the right
// TEE views follow the TEE convention (probe at the top of the sector, patient's
// left on the right at 0 deg).
//
// A probe is { pos, dir, lat, normal, target }: dir = beam centre, lat = screen
// right, normal = dir x lat (right-handed). Depths are real clinical settings.

import { LM, BODY_AX, anatomyParams, epiDist, lumenDist, PERI } from './anatomy.js';
import { norm, vadd, vsub, vscale, vdot, vcross, vrot } from './mathutils.js';

// Build a probe that looks along `dir` at `target`, with screen-right as close
// to `right` as the imaging plane allows, `standoff` cm back from the target.
export function probeFrom(target, dir, right, standoff) {
  const d = norm(dir);
  const lat = norm(vsub(right, vscale(d, vdot(right, d))));
  const normal = vcross(d, lat);
  return { pos: vsub(target, vscale(d, standoff)), dir: d, lat, normal, target: target.slice() };
}
// Parasternal windows look through ~2 cm of chest wall (skin, fat, pectoralis,
// intercostal muscle) before the first cardiac structure. Rather than a fixed
// standoff, back the probe off along the beam until the shallowest point of the
// heart or pericardium ANYWHERE in the sector lies at least `wall` cm deep.
const _A0 = anatomyParams(0, 0, {}, null, 0);
function firstHeartDepth(pos, dir, lat, halfAngle, maxDepth) {
  let best = 1e9;
  for (let a = -halfAngle; a <= halfAngle + 1e-9; a += halfAngle / 8) {
    const c = Math.cos(a), s = Math.sin(a);
    const d = [dir[0] * c + lat[0] * s, dir[1] * c + lat[1] * s, dir[2] * c + lat[2] * s];
    for (let t = 0; t < Math.min(best, maxDepth); t += 0.05) {
      const p = [pos[0] + d[0] * t, pos[1] + d[1] * t, pos[2] + d[2] * t];
      if (epiDist(p[0], p[1], p[2], _A0) < PERI.gap + PERI.thick) { best = t; break; }
    }
  }
  return best;
}
export function chestWallProbe(target, dir, right, standoff, wall) {
  let pr = probeFrom(target, dir, right, standoff);
  for (let i = 0; i < 3; i++) {
    const d0 = firstHeartDepth(pr.pos, pr.dir, pr.lat, 0.66, standoff + 4);
    if (Math.abs(d0 - wall) < 0.05) break;
    standoff += wall - d0;
    pr = probeFrom(target, dir, right, standoff);
  }
  return pr;
}
// Legacy constructor kept for the debugging hook: plane normal given explicitly
// (screen-right = normal x dir).
export function makeProbe(target, dirGuess, normal, standoff) {
  const n = norm(normal);
  const dir = norm(vsub(dirGuess, vscale(n, vdot(dirGuess, n))));
  const lat = norm(vcross(n, dir));
  return { pos: vsub(target, vscale(dir, standoff)), dir, lat, normal: n, target: target.slice() };
}

const M = LM.M, T = LM.T, A = LM.A;
const Y = [0, 1, 0];
// the true (epicardial) apex, where the apical probe sits
const APEX = [0, LM.LV.apexY - LM.LV.wall * LM.LV.apexWallFrac, 0];
const CRUX = vscale(vadd(M, T), 0.5);
const CRUX0 = CRUX;
// PLAX / APLAX plane: contains the LV long axis and the aortic root.
const D_AO = norm([A[0], 0, A[2]]);                 // horizontal direction LV axis -> aortic valve
// chest-wall thickness between the skin and the first cardiac structure
const WALL = 2.0;

// The A4C is aimed just posterior (inferior, heart -z) of the crux: tilted off the
// LVOT and the aortic root, so the upper atrial septum is not cut tangentially
// against the non-coronary sinus (which would read as a solid crux 'mass').
const A4C_TILT = LM.A4C_TILT;
function a4cProbe() {
  const CRUX = vadd(CRUX0, [0, 0, -A4C_TILT]);
  const dir = norm(vsub(CRUX, APEX));
  // screen right = LV side (from the tricuspid toward the mitral)
  return probeFrom(CRUX, dir, vsub(M, T), vdot(vsub(CRUX, APEX), dir) + 1.4);
}
function a2cProbe() {
  // plane through the long axis and the anterior / inferior walls, rotated ~15
  // deg toward lateral from 12 o'clock — as a sonographer rotates from the A4C
  // until the RV and RVOT drop out of the image (anterior wall on the right)
  const tgt = [0, M[1] - 0.2, 0];
  const dir = norm(vsub(tgt, APEX));
  return probeFrom(tgt, dir, [0.26, 0, 0.97], vdot(vsub(tgt, APEX), dir) + 1.4);
}
function a5cProbe() {
  const tgt = vscale(vadd(A, CRUX), 0.5);
  const dir = norm(vsub(tgt, APEX));
  return probeFrom(tgt, dir, vsub(M, T), vdot(vsub(tgt, APEX), dir) + 1.4);
}
function a3cProbe() {
  // APLAX: the same plane as PLAX (apex + mitral centre + aortic valve centre),
  // imaged from the apex — so the RV drops out; aortic outflow on the right
  const tgt = vscale(vadd(A, M), 0.5);
  const dir = inPlane(vsub(tgt, APEX), PLAX_N);
  return probeFrom(tgt, dir, inPlane(D_AO, PLAX_N), vdot(vsub(tgt, APEX), dir) + 1.4);
}
// The parasternal window is on the anterior chest wall, so the parasternal
// beams run POSTERIORLY: body-anterior projected into each imaging plane.
// The textbook PLAX plane contains the LV apex and the centres of BOTH the mitral
// (through A2/P2) and the aortic valve (through the RCC / NCC).
const PLAX_N = (() => {
  const n = norm(vcross(vsub(M, APEX), vsub(A, APEX)));
  return vdot(n, vcross(Y, D_AO)) < 0 ? vscale(n, -1) : n;
})();
const inPlane = (v, n) => norm(vsub(v, vscale(n, vdot(v, n))));
function plaxProbe() {
  // aorto-mitral curtain, image centre (projected into the plane)
  const t0 = [-0.35, 0.3, 0.35];
  const tgt = vsub(t0, vscale(PLAX_N, vdot(vsub(t0, M), PLAX_N)));
  // Beam runs posteriorly from the chest wall, tipped ~20 deg toward the base —
  // the sonographer's adjustment (one interspace up / probe angulation) that lays
  // the LV long axis near-horizontal on screen with the Ao/LA stacked at right.
  const post = vscale(inPlane(BODY_AX.A, PLAX_N), -1);
  const yIn = inPlane(Y, PLAX_N);                   // long axis within the plane
  const perp = inPlane(vsub(post, vscale(yIn, vdot(post, yIn))), PLAX_N); // posterior, perpendicular to it
  const dir = norm(vadd(vscale(perp, Math.cos(0.35)), vscale(yIn, Math.sin(0.35))));
  // base / aorta to the right of the screen (ASE)
  return chestWallProbe(tgt, dir, yIn, 5.0 + WALL, WALL);
}
function psaxProbe(y) {
  const tgt = [0, y, 0];
  // beam posterior in the short-axis plane; the probe is rotated (index toward the
  // left shoulder) until the LV reads in its standard orientation — anterior wall at
  // 12, papillary muscles at ~4 and ~8 o'clock — i.e. the heart's own anterior axis
  const dir = vscale(inPlane(LM.ANT_HEART, Y), -1);
  // viewed from the apex: patient-left (LV lateral wall, +x) on the right
  const right = vcross(Y, dir);
  return chestWallProbe(tgt, dir, vdot(right, [1, 0, 0]) > 0 ? right : vscale(right, -1), 4.8 + WALL, WALL);
}
function psaxAvProbe() {
  const tgt = vadd(A, vscale(LM.U_AO, 0.75));        // cusp coaptation level
  return chestWallProbe(tgt, vscale(LM.E_ANT, -1), LM.E_SCR, 4.3 + WALL, WALL);
}
const SC_W = [0.7, 0.5, 0.4];
const SC_WALL = 6.5;                                 // subcostal skin -> first cardiac structure (cm)
function subcostalProbe() {
  // Subcostal four-chamber: the four-chamber plane that cuts the interatrial
  // septum at right angles through the fossa and runs to the apex (so the IAS lies
  // perpendicular to the beam — the ASD view), imaged from BELOW through the liver:
  // the sub-xiphoid beam aimed up toward the left shoulder, taken within that
  // plane. Index marker to the patient's left, so the apex points right.
  const S = BODY_AX.S, P = BODY_AX.P, L = BODY_AX.L;
  const fossa = LM.IAS.fossaC;
  // the beam runs from under the xiphoid up, back and to the left — close to the
  // septal normal (RA -> LA), so the septum crosses it near-perpendicular — and
  // the plane is the one through that beam and the apex (all four chambers)
  const raw = vadd(vadd(vscale(S, SC_W[0]), vscale(P, SC_W[1])), vscale(L, SC_W[2]));
  const tgt0 = vadd(fossa, vscale(vsub(APEX, fossa), 0.35));
  const n = norm(vcross(raw, vsub(APEX, tgt0)));
  const dir = norm(raw);
  let right = vcross(n, dir);
  if (vdot(right, L) < 0) right = vscale(right, -1);
  const tgt = vadd(fossa, vscale(vsub(APEX, fossa), 0.35));
  // skin-to-heart depth: ~2 cm of abdominal wall and 3-5 cm of left lobe put
  // the first cardiac structure ~6 cm deep, the whole heart in the mid field
  return chestWallProbe(tgt, dir, right, 14.0, SC_WALL);
}

// apical windows sit over the apex: a remodelled (longer) LV moves its apex, and
// the apical window with it, out toward the left axilla / caudally
const APICAL_TRACK = (A) => [0, A.lv.apexShift || 0, 0];
// A4C: as the root descends in systole the sonographer's small posterior tilt
// keeps the plane off it (as at end-diastole), so the root-clearance shell of the
// atria is never brought into the upper septum as a false crux 'mass'
const A4C_ROOT_TRACK = 0.25;                         // cm of posterior tilt per cm of root descent
let _a4cN = null;
function a4cTrack(A) {
  if (!_a4cN) _a4cN = a4cProbe().normal;
  const back = vscale([0, 0, 1], A4C_ROOT_TRACK * (A.lv.A[1] - A0.lv.A[1]));
  return vadd(APICAL_TRACK(A), vscale(_a4cN, vdot(back, _a4cN)));
}
export const TTE_VIEWS = {
  PLAX: { probe: plaxProbe, depth: 16 },
  PSAX: { probe: () => psaxProbe(-3.6), depth: 15 },   // mid-papillary
  PSAX_AV: { probe: psaxAvProbe, depth: 15, track: (A) => [0, 0.8 * (A.lv.A[1] - A0.lv.A[1]), 0] },   // follows the root's descent
  A4C: { probe: a4cProbe, depth: 17, track: a4cTrack },
  A2C: { probe: a2cProbe, depth: 17, track: APICAL_TRACK },
  SUBCOSTAL: { probe: subcostalProbe, depth: 18 },
};
// RV inflow (parasternal, tilted from the PLAX toward the right hip): RA, the
// tricuspid (anterior and posterior leaflets) and the RV inflow, with the
// coronary sinus / IVC entering the RA. RV in the near field, base on the right.
// From the PLAX window aimed at the tricuspid, the plane is rolled RVIT_ROLL about
// the beam and tilted RVIT_TILT, then slid until the posterior RA floor is in
// plane: the coronary-sinus ostium and the Eustachian valve at the IVC orifice,
// with the tricuspid and the RV inflow, and not the SVC. The CS ostium and the
// IVC lie ~1.3 cm apart (about the sum of their radii) across any plane from
// this window that keeps the RV clear of the septum, so the plane is set
// RVIT_CS cm from the CS ostium — the CS has priority, as it is what the view is
// for — and the IVC tube itself lies just off the plane.
const RVIT_ROLL = -12, RVIT_TILT = -8, RVIT_CS = 0.3;     // deg, deg, cm
const csOstium = (A) => A.cs[A.cs.length - 1].b;
let _rvit = null;
function rvInflowProbe() {
  if (!_rvit) {
    const pl = plaxProbe();
    const tgt = vadd(T, [0, 0.3, 0]);
    const d0 = norm(vsub(tgt, pl.pos));
    const raToApex = norm(vsub(A0.ra.c, [LM.T[0] - 0.5, LM.LV.apexY + 2.5, LM.T[2] + 1.0]));
    const base = probeFrom(tgt, d0, raToApex, vdot(vsub(tgt, pl.pos), d0));
    const lat = vrot(base.lat, base.dir, RVIT_ROLL * Math.PI / 180);
    const dir = vrot(base.dir, lat, RVIT_TILT * Math.PI / 180);
    const n = vcross(dir, lat);
    const shift = vscale(n, vdot(vsub(csOstium(A0), tgt), n) - RVIT_CS);
    const t2 = vadd(tgt, shift);
    _rvit = chestWallProbe(t2, dir, lat, vdot(vsub(tgt, pl.pos), d0), WALL);
  }
  return { ..._rvit, pos: _rvit.pos.slice(), target: _rvit.target.slice() };
}
// the RA floor rises and falls with the tricuspid annulus: the plane follows the
// CS ostium across the cycle (the small tilt a sonographer makes to hold it)
function rvInflowTrack(A) {
  const n = rvInflowProbe().normal;
  return vscale(n, vdot(vsub(csOstium(A), csOstium(A0)), n));
}
// Subcostal IVC long axis: the IVC running through the liver into the RA, the
// hepatic veins joining it; cranial (RA) to the right of the screen.
const SC_IVC_CRANIAL = 15;                        // deg
const SC_IVC_ROLL = -25;                          // deg
function subcostalIvcProbe() {
  // transducer at the subxiphoid window (the subcostal 4C skin point), aimed at
  // the IVC ~1.5 cm below its RA junction; the plane contains the IVC axis so the
  // vessel runs across the image through the liver into the RA (cranial right)
  const iv = A0.ivc;
  const ax = norm(vsub(iv.a, iv.b));                      // IVC axis, toward the RA
  // sagittal sub-xiphoid probe held back and ~15 deg cranial, so the IVC runs
  // ACROSS the image — caudal (liver) screen-left, entering the RA screen-right —
  // rather than diving into the far field. The cranial tilt puts the transducer
  // below the xiphoid over the liver (abdominal wall, then liver, in the near
  // field) instead of in front of the lung bases.
  const tgt = vadd(iv.a, vscale(ax, -0.8));
  const back = inPlane(vadd(BODY_AX.P, vscale(BODY_AX.S, 0.3)), vcross(ax, BODY_AX.P));
  let perp = norm(vsub(back, vscale(ax, vdot(back, ax))));
  // rolled about the IVC axis toward the patient's right, so above the caval
  // orifice the plane runs up through the RA body, clear of the RV in front of it
  const roll = SC_IVC_ROLL * Math.PI / 180;
  perp = norm(vadd(vscale(perp, Math.cos(roll)), vscale(vcross(ax, perp), Math.sin(roll))));
  const th = SC_IVC_CRANIAL * Math.PI / 180;
  const dir = norm(vadd(vscale(perp, Math.cos(th)), vscale(ax, Math.sin(th))));
  return chestWallProbe(tgt, dir, ax, 9.0, 6.0);
}
// Suprasternal notch: the aortic arch in long axis, the ascending aorta on the
// left of the screen and the descending aorta on the right, with the right
// pulmonary artery in cross-section beneath the arch.
function suprasternalProbe() {
  const asc = LM.AO_ASC_TOP, arch = LM.AO_ARCH_MID, desc = LM.DTA_TOP;
  const n = norm(vcross(vsub(asc, arch), vsub(desc, arch)));
  const tgt = vscale(vadd(vadd(asc, desc), vscale(arch, 2)), 0.25);
  const dir = inPlane(vscale(BODY_AX.S, -1), n);
  return probeFrom(tgt, dir, vsub(desc, asc), 7.0);
}
// additional apical / parasternal / subcostal / suprasternal views
export const EXTRA_VIEWS = {
  A5C: { probe: a5cProbe, depth: 17, track: APICAL_TRACK },
  A3C: { probe: a3cProbe, depth: 17, track: APICAL_TRACK },
  // mitral leaflet level, 0.6 cm below the annulus (closed coaptation line at ED,
  // fish-mouth in diastole): the probe follows the posterolateral annular
  // excursion (the small tilt a sonographer makes), so the plane stays at the
  // leaflets and never climbs into the posterior AV groove (CS / RA) in systole
  PSAX_MV: { probe: () => psaxProbe(-0.6), depth: 15, track: (A) => [0, 0.95 * (A.lv.base - A0.lv.base), 0] },
  RVIT: { probe: rvInflowProbe, depth: 13, track: rvInflowTrack },
  SC_IVC: { probe: subcostalIvcProbe, depth: 18 },
  SSN: { probe: suprasternalProbe, depth: 16 },
};

// ---- TEE -------------------------------------------------------------------
// The probe sits in the oesophagus directly behind the left atrium (just in
// front of the descending aorta), so the atria are in the near field.
const A0 = anatomyParams(0, 0, {}, null, 0);
// oesophagus directly against the LA's posterior wall: step back from the LA centre
// along body-posterior until leaving the drawn LA lumen, then clear its wall
export const ESO = (() => {
  let t = 0;
  while (t < 6 && lumenDist(...vadd(A0.la.c, vscale(BODY_AX.P, t)), A0, 'LA') < 0) t += 0.05;
  // LA posterior wall (~0.2) + oesophageal wall and periaortic fat (~0.6): the LA
  // blood pool starts ~0.8-1 cm below the transducer face
  return vadd(A0.la.c, vscale(BODY_AX.P, t + 0.8));
})();
// every mid-oesophageal probe keeps the oesophageal wall in front of the LA across
// the whole sector: back the transducer off along the beam until no beam reaches
// LA blood within the first 0.7 cm
function teeClear(pr) {
  for (let it = 0; it < 30; it++) {
    let near = 1e9;
    for (let a = -0.5; a <= 0.501; a += 0.125) {
      const d = vadd(vscale(pr.dir, Math.cos(a)), vscale(pr.lat, Math.sin(a)));
      for (let r = 0; r <= 0.7; r += 0.1) {
        if (lumenDist(...vadd(pr.pos, vscale(d, r)), A0, 'LA') < 0) { near = Math.min(near, r); break; }
      }
    }
    if (near > 0.7) break;
    pr.pos = vsub(pr.pos, vscale(pr.dir, 0.1));
  }
  return pr;
}
const ME4C_RETRO = 0.4, ME4C_SWING = 1, ME4C_UP = 2;    // cm, deg, cm
function me4cProbe() {
  // The plane that cuts the interatrial septum at right angles through the fossa
  // ovalis and runs down to the apex (both atria, the septum and any ASD / PFO in
  // profile, both ventricles); the probe sits where the oesophagus meets it.
  // (0 deg): the plane through both AV-valve centres and the true apex, so the
  // crux, the atrial septum and both ventricles lie in plane with their true
  // proportions. Retroflexed slightly: the base end of the plane is tipped
  // ME4C_RETRO cm posterior, off the LVOT, so the non-coronary sinus is not cut
  // tangentially against the upper septum (a crux 'mass') and the plane takes
  // in the posterior (larger) part of the RA and the fossa ovalis.
  const Mp = vadd(M, [0, 0, -ME4C_RETRO]), Tp = vadd(T, [0, 0, -ME4C_RETRO]);
  const n = norm(vcross(vsub(Tp, Mp), vsub(APEX, Mp)));
  // (the transducer withdrawn ME4C_UP cm up the oesophagus from the LA-level
  // reference: the beam then runs further down the LV long axis, so the apex is
  // not foreshortened against the right edge of the sector)
  const eso = vadd(ESO, vscale(BODY_AX.S, ME4C_UP));
  const pos = vsub(eso, vscale(n, vdot(vsub(eso, Mp), n)));
  const crux = vscale(vadd(Mp, Tp), 0.5);
  const tgt = vadd(crux, vscale(vsub(APEX, crux), 0.2));
  let dir = norm(vsub(tgt, pos));
  let right = vcross(n, dir);
  if (vdot(right, vsub(M, T)) < 0) right = vscale(right, -1);   // LV on the right
  // the beam swung ME4C_SWING toward the RA about the transducer (in plane): just
  // enough to keep the RA and the fossa inside the left edge while the LV long axis
  // stays within ~30 deg of the beam and the whole apex well inside the right edge
  const sw = ME4C_SWING * Math.PI / 180;
  dir = norm(vsub(vscale(dir, Math.cos(sw)), vscale(right, Math.sin(sw))));
  const sd = vdot(vsub(tgt, pos), dir);
  return teeClear(probeFrom(vadd(pos, vscale(dir, sd)), dir, right, sd));
}
function melaxProbe() {
  const tgt = vscale(vadd(A, M), 0.5);
  const dir = norm(vsub(tgt, ESO));
  // ME LAX (~130 deg): LVOT / aorta to the right of the screen
  return teeClear(probeFrom(tgt, dir, vadd(Y, vscale(D_AO, 0.6)), vdot(vsub(tgt, ESO), dir)));
}
function tgsaxProbe() {
  // transgastric: probe in the fundus below the heart, looking up through the
  // inferior wall; mid-papillary short axis, lateral wall on the right
  // (the transducer rests on the gastric wall; anatomy.js draws that wall in the
  // first ~0.45 cm, then the diaphragm and the inferior wall)
  const G = LM.GASTRIC;
  return probeFrom(vadd(G.p, vscale(G.n, 4.6)), G.n, [1, 0, 0], 4.6);
}
// ME LAA (~60-90 deg, turned left): the appendage along its length as a finger
// on the screen-right of the LA, the LSPV beside it across the warfarin ridge.
// The plane holds the oesophagus and the LAA ostium; its roll about that line is
// the one that lays the most appendage lumen in plane while keeping the LSPV
// ostium within 0.4 cm of it and the least ventricular myocardium in the sector
// at end-diastole and end-systole (the LV wall swings into a plane grazing the
// AV groove as the base descends).
const A_ES = anatomyParams(1, 0, {}, null, 0.42);
function laaArea(aa, pos, dir, lat, step = 0.15) {
  // in-plane lumen of the appendage segments (round cones), sampled on the plane
  let n = 0;
  const c = aa[1].a;
  const s0 = vdot(vsub(c, pos), dir), l0 = vdot(vsub(c, pos), lat);
  for (let a = -3; a <= 3; a += step) for (let b = -3; b <= 3; b += step) {
    if (s0 + a < 0.3 || Math.abs(l0 + b) > (s0 + a) * 0.75) continue;   // inside the sector
    const q = vadd(pos, vadd(vscale(dir, s0 + a), vscale(lat, l0 + b)));
    for (const g of aa) {
      const ab = vsub(g.b, g.a), t = Math.max(0, Math.min(1, vdot(vsub(q, g.a), ab) / vdot(ab, ab)));
      const r = g.r1 + (g.r2 - g.r1) * t;
      if (Math.hypot(...vsub(q, vadd(g.a, vscale(ab, t)))) < r) { n++; break; }
    }
  }
  return n * step * step;
}
function myoArea(A, pos, dir, lat, depth, step = 0.25) {
  let n = 0;
  for (let d = 0.5; d <= depth; d += step) for (let x = -d * 0.6; x <= d * 0.6; x += step) {
    const q = vadd(pos, vadd(vscale(dir, d), vscale(lat, x)));
    if (epiDist(...q, A) < 0 && lumenDist(...q, A, 'LV') > 0 && lumenDist(...q, A, 'LA') > 0) n++;
  }
  return n * step * step;
}
let _melaa = null;
// the appendage rides down with the AV plane in systole: the operator follows it
// (a slight probe flex), so the plane shifts along its normal with the LAA
function melaaTrack(A) {
  const p = melaaProbe(), n = p.normal;
  const c = vscale(vadd(A.la.aa[0].b, A.pv.veins[0].w), 0.5), c0 = vscale(vadd(A0.la.aa[0].b, A0.pv.veins[0].w), 0.5);
  return vscale(n, vdot(vsub(c, c0), n));
}
function melaaProbe() {
  if (_melaa) return { ..._melaa, pos: _melaa.pos.slice(), target: _melaa.target.slice() };
  // plane through the oesophagus that best holds the appendage's axis and the
  // LSPV's (both leave the lateral wall ~1-1.5 cm apart, the ridge between)
  const aa = A0.la.aa, v = A0.pv.veins[0], o = aa[0].a, tip = aa[2].b;
  const vAx = norm(vsub(v.b, v.w));
  const pts = [[o, 2], [aa[0].b, 2], [aa[1].b, 1], [v.w, 2], [vadd(v.w, vscale(vAx, 1.2)), 1]];
  let best = null;
  for (let th = 0; th < Math.PI; th += Math.PI / 90) for (let ph = 0; ph < Math.PI; ph += Math.PI / 90) {
    const n = [Math.sin(th) * Math.cos(ph), Math.sin(th) * Math.sin(ph), Math.cos(th)];
    let cost = 0;
    for (const [q, w] of pts) cost += w * vdot(vsub(q, ESO), n) ** 2;
    if (!best || cost < best.cost) best = { cost, n };
  }
  // fine roll about the oesophagus -> appendage line: least ventricular
  // myocardium in the sector at end-diastole and end-systole, the LSPV kept in plane
  const ax = norm(vsub(aa[0].b, ESO));
  let pick = null;
  for (let deg = -15; deg <= 15; deg += 3) {
    const n = vrot(best.n, ax, deg * Math.PI / 180);
    if (Math.abs(vdot(vsub(v.w, ESO), n)) > 0.4) continue;
    const tgt = vscale(vadd(aa[1].a, vadd(v.w, vscale(vAx, 1.0))), 0.5);
    const dir = norm(vsub(vsub(tgt, ESO), vscale(n, vdot(vsub(tgt, ESO), n))));
    const lat = vcross(n, dir);
    const myo = Math.max(myoArea(A0, ESO, dir, lat, 10), myoArea(A_ES, ESO, dir, lat, 10));
    const score = laaArea(aa, ESO, dir, lat) - 0.5 * myo;
    if (!pick || score > pick.score) pick = { score, n, dir, tgt };
  }
  const { n, dir, tgt } = pick || { n: best.n, tgt: vscale(vadd(aa[1].a, vadd(v.w, vscale(vAx, 1.0))), 0.5) };
  const d = dir || norm(vsub(vsub(tgt, ESO), vscale(n, vdot(vsub(tgt, ESO), n))));
  let right = vcross(n, d);
  if (vdot(right, vsub(tip, o)) < 0) right = vscale(right, -1);    // appendage tip on screen-right
  _melaa = teeClear(probeFrom(tgt, d, right, vdot(vsub(tgt, ESO), d)));
  return melaaProbe();
}
let _me2c = null;
// ME two-chamber (~60-90 deg): LA near field, LV below, the LAA and anterior
// wall on the right, the inferior wall and coronary sinus on the left.
function me2cProbe() {
  const tgt = vscale(vadd(M, APEX), 0.5);
  const dir0 = norm(vsub(tgt, ESO));
  const l0 = norm(vsub([0.26, 0, 0.97], vscale(dir0, vdot([0.26, 0, 0.97], dir0))));
  const w = vcross(dir0, l0);
  // roll the plane about the beam (from the [0.26, 0, 0.97] reference) to the
  // angle that lays the most appendage lumen in plane at the worse of
  // end-diastole and end-systole: the LAA is then a finger on the right of the
  // screen, while the CS stays on the left; and swing the beam up to 12 deg
  // toward the appendage within the plane (the apex then sits a little left of
  // centre), so the LAA lies inside the sector
  if (_me2c) return { ..._me2c, pos: _me2c.pos.slice(), target: _me2c.target.slice() };
  let right = l0, dir = dir0, best = -1;
  for (let deg = -10; deg <= 50; deg += 2) {
    const th = deg * Math.PI / 180;
    const r = vadd(vscale(l0, Math.cos(th)), vscale(w, Math.sin(th)));
    for (let sw = 0; sw <= 12; sw += 2) {
      const d = norm(vadd(vscale(dir0, Math.cos(sw * Math.PI / 180)), vscale(r, Math.sin(sw * Math.PI / 180))));
      const lat = norm(vsub(r, vscale(d, vdot(r, d))));
      const a = Math.min(laaArea(A0.la.aa, ESO, d, lat), laaArea(A_ES.la.aa, ESO, d, lat)) - 0.02 * sw;
      if (a > best + 1e-6) { best = a; right = lat; dir = d; }
    }
  }
  _me2c = teeClear(probeFrom(vadd(ESO, vscale(dir, vdot(vsub(tgt, ESO), dir))), dir, right, vdot(vsub(tgt, ESO), dir)));
  return me2cProbe();
}
// ME aortic-valve short axis (~30-45 deg): the three cusps en face, the LA in
// the near field, the RVOT in the far field; NCC beside the interatrial septum.
function meAvSaxProbe() {
  const tgt = vadd(A, vscale(LM.U_AO, 0.75));
  const raw = vsub(tgt, ESO);
  const dir = norm(vsub(raw, vscale(LM.U_AO, vdot(raw, LM.U_AO))));
  return teeClear(probeFrom(tgt, dir, LM.E_SCR, vdot(raw, dir)));
}
// ME bicaval (~90-110 deg): the RA with the SVC entering on the right of the
// screen and the IVC on the left, the interatrial septum and fossa in profile
// between the LA (near field) and the RA — the view for PFO / sinus venosus ASD.
// The transducer stays in the oesophagus: it slides along it (BODY_AX.S) or is
// backed off behind it, and the omniplane plane is turned and rolled about the
// beam. Of the planes within BICAVAL_FOSSA cm of the fossa ovalis centre (so the
// fossa membrane is cut and the septum lies between the LA and the RA), the one
// holding the most of the first 3 cm of BOTH cavae is kept; the beam is then
// aimed through the LA at the septum just above the fossa.
const BICAVAL_SLIDE = [-3, 2.5];                 // search along the oesophagus (cm, + = withdrawn)
const BICAVAL_BACK = 1.5;                        // ... and behind it (cm)
const BICAVAL_LAT = 1.0;                         // ... and across it, left/right (cm)
const BICAVAL_FOSSA = 0.45;                      // fossa centre within this of the plane (cm)
const BICAVAL_AIM = 1.2;                         // beam aimed this far above the fossa (cm)
const BICAVAL_HALF = 0.45;                       // (its axis within this share of its radius of the plane)
const BICAVAL_CAVA = 2.2;                        // each cava's first 3 cm: at least this much in plane (cm)
const BICAVAL_SEPT_W = 1.0;                      // septal length traded 1:1 for caval length
let _bicaval = null;
function meBicavalProbe() {
  if (!_bicaval) {
    const along = (V, t) => vadd(V.a, vscale(norm(vsub(V.b, V.a)), t));
    const fos = LM.IAS.fossaC;
    // in-plane share of a cava's first 3 cm (its axis within ~half a radius of
    // the plane); 0 when the orifice itself is out of plane
    const inPl = (V, pos, n) => {
      let l = 0;
      for (let t = 0; t <= 3.001; t += 0.1) {
        if (Math.abs(vdot(vsub(along(V, t), pos), n)) < V.r * BICAVAL_HALF) l += 0.1;
        else if (t < 0.5) return 0;
      }
      return l;
    };
    // septum in plane: length of the line where the plane crosses the septal
    // plane that has LA blood on one side and RA blood on the other
    const N = LM.IAS_N;
    const septLen = (pos, n) => {
      const L = norm(vcross(n, N));
      let q0 = vsub(fos, vscale(n, vdot(vsub(fos, pos), n)));
      q0 = vsub(q0, vscale(N, vdot(vsub(q0, LM.IAS_P), N)));
      let l = 0;
      for (let t = -4; t <= 4; t += 0.1) {
        const q = vadd(q0, vscale(L, t));
        if (lumenDist(...vsub(q, vscale(N, 0.3)), A0, 'LA') < 0 && lumenDist(...vadd(q, vscale(N, 0.3)), A0, 'RA') < 0) l += 0.1;
      }
      return l;
    };
    const cands = [], fallback = [];
    for (let s = BICAVAL_SLIDE[0]; s <= BICAVAL_SLIDE[1] + 1e-9; s += 0.25) {
      for (let b = 0; b <= BICAVAL_BACK + 1e-9; b += 0.25) for (let lt = -BICAVAL_LAT; lt <= BICAVAL_LAT + 1e-9; lt += 0.5) {
        const pos = vadd(vadd(vadd(ESO, vscale(BODY_AX.S, s)), vscale(BODY_AX.P, b)), vscale(BODY_AX.L, lt));
        const r0 = norm(vsub(fos, pos));
        const e2 = vcross(r0, norm(vcross(r0, BODY_AX.S)));
        for (let ta = -0.15; ta <= 0.151; ta += 0.01) {
          const ax = norm(vadd(r0, vscale(e2, ta)));          // an in-plane ray near the fossa
          const a1 = norm(vcross(ax, BODY_AX.S)), a2 = vcross(ax, a1);
          for (let ph = 0; ph < Math.PI; ph += 0.02) {
            const n = vadd(vscale(a1, Math.cos(ph)), vscale(a2, Math.sin(ph)));
            const fo = Math.abs(vdot(vsub(fos, pos), n));
            if (fo > BICAVAL_FOSSA) continue;
            const cv = Math.min(inPl(A0.svc, pos, n), inPl(A0.ivc, pos, n));
            if (cv >= BICAVAL_CAVA) cands.push({ cv, fo, b, pos, n });
            else if (!fallback.length || cv > fallback[0].cv) fallback[0] = { cv, fo, b, pos, n };
          }
        }
      }
    }
    // among those holding enough of both cavae, the longest septum (then the
    // fossa nearest the plane, the probe nearest the oesophagus)
    let best = null;
    if (!cands.length) cands.push(...fallback);
    for (const c of cands) {
      c.sc = Math.min(c.cv, 2.5) + BICAVAL_SEPT_W * septLen(c.pos, c.n) - 0.3 * c.fo - 0.05 * c.b;
      if (!best || c.sc > best.sc) best = c;
    }
    const { pos, n } = best;
    const onPlane = (q) => vsub(q, vscale(n, vdot(vsub(q, pos), n)));
    const up = norm(vsub(onPlane(A0.svc.a), onPlane(A0.ivc.a)));
    const tgt = vadd(onPlane(fos), vscale(up, BICAVAL_AIM));
    const dir = norm(vsub(tgt, pos));
    let right = vcross(n, dir);
    if (vdot(right, up) < 0) right = vscale(right, -1);                     // SVC screen-right
    _bicaval = teeClear(probeFrom(tgt, dir, right, vdot(vsub(tgt, pos), dir)));
  }
  return { ..._bicaval, pos: _bicaval.pos.slice(), target: _bicaval.target.slice() };
}
// ME RV inflow-outflow (~60-75 deg): RA and tricuspid on the left, the RV
// wrapping round the aortic valve to the RVOT and pulmonary valve on the right.
function meRvioProbe() {
  // the plane through the oesophagus, the tricuspid and the pulmonary valve: the
  // RV wraps from its inflow (screen-left) round the aortic valve to the outflow
  const n = norm(vcross(vsub(T, ESO), vsub(LM.PV, ESO)));
  const tgt = vsub(A, vscale(n, vdot(vsub(A, ESO), n)));
  const dir = norm(vsub(tgt, ESO));
  let right = vsub(LM.PV, T);
  right = vsub(right, vscale(dir, vdot(right, dir)));
  return teeClear(probeFrom(tgt, dir, right, vdot(vsub(tgt, ESO), dir)));
}
// Descending thoracic aorta short axis: the probe turned to face posteriorly;
// the round aorta lies immediately behind the oesophagus.
function descAoSaxProbe() {
  const tgt = LM.DTA_P;
  const raw = vsub(tgt, ESO);
  const dir = norm(vsub(raw, vscale(BODY_AX.S, vdot(raw, BODY_AX.S))));
  const pr = probeFrom(tgt, dir, BODY_AX.L, Math.max(1.0, vdot(raw, dir)));
  // facing the aorta, the LA lies behind and beside the transducer: close in on the
  // aortic wall (>= 0.5 cm stand-off) until no beam of the sector reaches LA blood
  const laInSector = () => {
    for (let a = -0.7; a <= 0.701; a += 0.1) {
      const d = vadd(vscale(pr.dir, Math.cos(a)), vscale(pr.lat, Math.sin(a)));
      for (let r = 0; r <= 3.5; r += 0.1) if (lumenDist(...vadd(pr.pos, vscale(d, r)), A0, 'LA') < 0) return true;
    }
    return false;
  };
  for (let it = 0; it < 30 && laInSector() && vdot(vsub(tgt, pr.pos), pr.dir) > LM.DTA_R + 0.6; it++) {
    pr.pos = vadd(pr.pos, vscale(pr.dir, 0.05));
  }
  return pr;
}
export const TEE_VIEWS = {
  ME4C: { probe: me4cProbe, depth: 14 },
  ME2C: { probe: me2cProbe, depth: 14 },
  MELAX: { probe: melaxProbe, depth: 13 },
  MEAVSAX: { probe: meAvSaxProbe, depth: 10, track: (A) => [0, 0.8 * (A.lv.A[1] - A0.lv.A[1]), 0] },   // follows the root
  MEBICAVAL: { probe: meBicavalProbe, depth: 11 },
  MERVIO: { probe: meRvioProbe, depth: 12 },
  DESCAO: { probe: descAoSaxProbe, depth: 6 },
  TGSAX: { probe: tgsaxProbe, depth: 12 },
  // Mid-oesophageal left-atrial-appendage view: the plane CONTAINS the appendage's
  // long axis, so the narrow ostium, neck and hooked lobes lie in-plane.
  MELAA: { probe: melaaProbe, depth: 10, track: melaaTrack },
};

export const ALL_VIEWS = { ...TTE_VIEWS, ...EXTRA_VIEWS, ...TEE_VIEWS };
