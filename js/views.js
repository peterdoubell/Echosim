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

import { LM, BODY_AX, BODY, anatomyParams, epiDist, lumenDist, bodyClassify, PERI, stomachDist } from './anatomy.js';
import { norm, vadd, vsub, vscale, vdot, vcross, vrot } from './mathutils.js';
import { VALVE_DEFS, geometryAt } from './cardiac-model.js';

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
// ... while its apical end is aimed A4C_APEX_FWD in front of the true apex: the
// mitral annulus lies posterior to the LV long axis, so a plane from the true apex
// to the crux runs behind the mid-cavity and through the posterior papillary
// heads (a slit-like end-systolic cavity crossed by a muscle bar); this one cuts
// the mid and apical LV through its long axis, as the A2C and A3C do.
const A4C_APEX_FWD = 0.5;                        // cm, heart +z (anterior)
const A4C_APEX_LAT = 1.2;                        // cm, toward screen-right
function a4cProbe() {
  const CRUX = vadd(CRUX0, [0, 0, -A4C_TILT]);
  const AP = vadd(APEX, [0, 0, A4C_APEX_FWD]);
  const dir = norm(vsub(CRUX, AP));
  // screen right = LV side (from the tricuspid toward the mitral)
  const pr = probeFrom(CRUX, dir, vsub(M, T), vdot(vsub(CRUX, AP), dir) + 1.4);
  // the transducer sits A4C_APEX_LAT toward the LV lateral wall within that plane
  // (on the apical impulse, which lies lateral to the true apex), still aimed at
  // the crux: the mitral CW line then meets an eccentric MR jet at the same angle
  // as from the true apex
  const pos = vadd(pr.pos, vscale(pr.lat, A4C_APEX_LAT));
  return probeFrom(CRUX, norm(vsub(CRUX, pos)), vsub(M, T), Math.hypot(...vsub(CRUX, pos)));
}
function a2cProbe() {
  // plane through the long axis and the anterior / inferior walls, rotated ~15
  // deg toward lateral from 12 o'clock — as a sonographer rotates from the A4C
  // until the RV and RVOT drop out of the image (anterior wall on the right)
  const tgt = [0, M[1] - 0.2, 0];
  const dir = norm(vsub(tgt, APEX));
  return probeFrom(tgt, dir, [0.26, 0, 0.97], vdot(vsub(tgt, APEX), dir) + 1.4);
}
// A5C: from the A4C the transducer is tilted anteriorly until the LVOT opens up
// with the aortic valve and root: the beam runs from the apex to the aortic
// annulus (LVOT PW / AV CW aligned with it) and the plane, from the A4C's
// tricuspid -> mitral direction, is turned A5C_TURN about the beam toward the
// root, so the root axis lies (within ~3 deg) in plane: the LVOT channel between
// the septum and the anterior mitral leaflet, the cusps and 3.5+ cm of root and
// ascending aorta, with the RV still on the left and the atria in the far field.
// (14 cm deep: the root and the LA behind it, with at most one small section of
// the descending aorta, not the arch and the PA branches beyond them.)
const A5C_TURN = 15;                                 // deg
function a5cProbe() {
  const dir = norm(vsub(A, APEX));
  const right = vrot(inPlane(vsub(M, T), dir), dir, A5C_TURN * Math.PI / 180);
  return probeFrom(A, dir, right, vdot(vsub(A, APEX), dir) + 1.4);
}
// the plane follows the root as it descends in systole (the small tilt that keeps
// the LVOT and the cusps open on screen), on top of the apical window's own track
let _a5cN = null;
function a5cTrack(A) {
  if (!_a5cN) _a5cN = a5cProbe().normal;
  return vadd(APICAL_TRACK(A), vscale(_a5cN, vdot(vsub(A.lv.A, A0.lv.A), _a5cN)));
}
// The beam is angled A3C_STEER within the plane, away from the anteroseptal
// side, from the same apical transducer position: the band of RV in front of the
// anteroseptal wall leaves the right edge (<= 5 % of the sector) while the LVOT
// and the root stay in the sector. (15 cm deep: the LA and the root, with the
// descending aorta behind the LA; the SVC and the right PA, which cross this
// plane beyond the LA roof at 15-17 cm, are not in it.)
const A3C_STEER = -6;                                // deg, in plane
function a3cProbe() {
  // APLAX: the same plane as PLAX (apex + mitral centre + aortic valve centre),
  // imaged from the apex — so the RV drops out; aortic outflow on the right
  const tgt = vscale(vadd(A, M), 0.5);
  const dir = inPlane(vsub(tgt, APEX), PLAX_N);
  const right = inPlane(D_AO, PLAX_N), sd = vdot(vsub(tgt, APEX), dir) + 1.4;
  const pos = vsub(tgt, vscale(dir, sd)), d2 = vrot(dir, PLAX_N, A3C_STEER * Math.PI / 180);
  return probeFrom(vadd(pos, vscale(d2, sd)), d2, right, sd);
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
const PLAX_TIP = 0.35;                                // rad (~20 deg)
const PLAX_SLIDE = 2.0;                              // cm, probe slid toward the apex along the chest wall
function plaxProbe() {
  // aorto-mitral curtain, image centre (projected into the plane)
  const t0 = [-0.35, 0.3, 0.35];
  const tgt = vsub(t0, vscale(PLAX_N, vdot(vsub(t0, M), PLAX_N)));
  // Beam runs posteriorly from the chest wall, tipped PLAX_TIP toward the base —
  // the sonographer's adjustment (one interspace up / probe angulation) that lays
  // the LV long axis near-horizontal on screen with the Ao/LA stacked at right.
  // (~20 deg: the root lies ~47 deg from the LV axis in this plane, so without
  // the tip the root climbs steeply to the near field and the beam through the
  // aortic valve crosses no RVOT. The probe also sits PLAX_SLIDE cm toward the
  // apex (the 3rd-4th interspace rather than over the base), so the beam to the
  // aortic valve comes up through the infundibulum: 2-3 cm of RVOT in front of
  // the root at ED and ~2 cm at ES, the root ~27 deg off the chest wall, and the
  // right edge of the sector stays short of the SVC, which crosses this plane
  // beyond the ascending aorta.)
  const post = vscale(inPlane(BODY_AX.A, PLAX_N), -1);
  const yIn = inPlane(Y, PLAX_N);                   // long axis within the plane
  const perp = inPlane(vsub(post, vscale(yIn, vdot(post, yIn))), PLAX_N); // posterior, perpendicular to it
  const dir = norm(vadd(vscale(perp, Math.cos(PLAX_TIP)), vscale(yIn, Math.sin(PLAX_TIP))));
  // base / aorta to the right of the screen (ASE)
  return chestWallProbe(vsub(tgt, vscale(yIn, PLAX_SLIDE)), dir, yIn, 5.0 + WALL, WALL);
}
// The parasternal beams run posteriorly from the chest wall (as the PLAX's does):
// each short-axis beam is body-posterior taken within its imaging plane. The heart
// is rotated about its long axis in the chest (RV anterior, the inferior wall on
// the diaphragm), so this lays the RV crescent over the front of the LV (~8 to ~1
// o'clock) with the posterior pericardium and lung beyond the inferolateral wall;
// a beam down the heart's own anterior axis would instead run out through the
// diaphragmatic surface into the liver and the gastric fundus (a subcostal
// short axis, not a parasternal one).
const postIn = (n) => inPlane(BODY_AX.P, n);
// psaxProbe(y, sx, rot, tilt, roll): the short-axis plane at level y, its beam
// turned `rot` deg from the heart's anterior axis toward body-posterior (at most
// the full ~45 deg, where it is body-posterior within the plane), tipped `tilt`
// deg cephalad about screen-right (the probe angled up from the same interspace,
// so the lower sector clears the diaphragm), rolled `roll` deg about the beam,
// and the probe slid `sx` cm along the chest wall toward screen-right (the LV
// lateral wall) so the left edge of the sector runs through the RV, not on past
// its diaphragmatic wall into the liver.
function psaxProbe(y, sx, rot, tilt, roll) {
  const a0 = vscale(inPlane(LM.ANT_HEART, Y), -1), p0 = postIn(Y);
  const full = Math.acos(Math.max(-1, Math.min(1, vdot(a0, p0))));
  const ax = norm(vcross(a0, p0));
  let dir = vrot(a0, ax, Math.min(full, rot * Math.PI / 180));
  // viewed from the apex: patient-left (LV lateral wall, +x) on the right
  let right = vcross(Y, dir);
  if (vdot(right, [1, 0, 0]) < 0) right = vscale(right, -1);
  const tl = tilt * Math.PI / 180;
  if (tl) { const d1 = vrot(dir, right, tl), d2 = vrot(dir, right, -tl); dir = vdot(d1, BODY_AX.S) > vdot(d2, BODY_AX.S) ? d1 : d2; }
  if (roll) right = vrot(right, dir, roll * Math.PI / 180);
  return chestWallProbe(vadd([0, y, 0], vscale(right, sx)), dir, right, 4.8 + WALL, WALL);
}
// mid-papillary: 40 deg toward body-posterior, 13 deg cephalad, slid 0.7 cm laterally
const PSAX_MID = [-3.4, 0.7, 40, 13, 0];             // y, sx (cm), rot, tilt, roll (deg)
const PSAX_AV_H = 0.45, PSAX_AV_TILT = [-7.5, 0];    // cm along the root axis; deg, deg
function psaxAvProbe() {
  const tgt = vadd(A, vscale(LM.U_AO, PSAX_AV_H));   // just below the cusp coaptation level
  // the plane is the root's short axis (cusps en face); the beam, body-posterior
  // within it, runs on behind the LA into the posterior mediastinum, not down past
  // the RA into the liver
  let dir = postIn(LM.U_AO);
  let right = vcross(LM.U_AO, dir);
  if (vdot(right, LM.E_SCR) < 0) right = vscale(right, -1);
  // Tipped PSAX_AV_TILT[0] caudally about screen-right: the PA trunk runs back and
  // ~13 deg toward the apex from the root's short-axis plane (U_PA . U_AO = -0.22),
  // so the tipped plane lies along the trunk (PV, trunk and bifurcation in view at
  // every phase, the RVOT PW beam within ~19 deg of the valve axis) and its far
  // field cuts the LA body behind the root (~4 cm wide) rather than its roof. It is
  // not rolled about the beam: the first ~2.5 cm of the LPA (laid in the root's
  // short-axis plane) stays in the plane, so the trunk divides with both branch
  // origins in view ('trousers'). (The review-11 roll of 8 deg, the left edge
  // toward the base, carried the trunk and the infundibulum out of the plane
  // through systole; tipped further than ~8.5 deg the LPA leaves the plane.)
  dir = vrot(dir, right, PSAX_AV_TILT[0] * Math.PI / 180);
  right = vrot(right, dir, PSAX_AV_TILT[1] * Math.PI / 180);
  return chestWallProbe(tgt, dir, right, 4.3 + WALL, WALL);
}
const SC_W = [0.7, 0.5, 0.4];
const SC_TURN = -8;                                  // deg about the beam
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
  // turned SC_TURN about the beam, about the line through the fossa (which stays
  // in plane, so the septum and any ASD are still cut through the fossa): the
  // plane runs from the tricuspid hinge along the RV free wall to the RV apex, so
  // the RV is a full chamber under the liver (~13 cm2, ~half of it at ES)
  right = vrot(right, dir, SC_TURN * Math.PI / 180);
  const nn = vcross(dir, right);
  const tgt2 = vsub(tgt, vscale(nn, vdot(vsub(tgt, fossa), nn)));
  // skin-to-heart depth: ~2 cm of abdominal wall and 3-5 cm of left lobe put
  // the first cardiac structure ~6 cm deep, the whole heart in the mid field
  return chestWallProbe(tgt2, dir, right, 14.0, SC_WALL);
}

// apical windows sit over the apex: a remodelled (longer) LV moves its apex, and
// the apical window with it, out toward the left axilla / caudally
const APICAL_TRACK = (A) => [0, A.lv.apexShift || 0, 0];
// Mid-papillary short axis (PSAX): the papillary heads descend with the annulus
// in systole (the chordae stay taut), so the plane follows `k` of their mean
// descent along its own normal and the muscles stay in plane through systole,
// larger and closer together at ES, at the mid-ventricle the sonographer holds the
// probe on. The probe also slides `lat` of the descent toward screen-right along
// the chest wall, so the far field of the sector does not dip below the heart into
// the liver as the plane moves apically. (Kept small: a larger slide puts the
// window over the lung edge, a larger track shows the liver.)
function papTrack(probe, k, lat = 0) {
  let pr = null;
  return (A) => {
    if (!pr) pr = probe();
    const pap = A.pap || [], drop = pap.length ? pap.reduce((s, p) => s + (p.drop || 0), 0) / pap.length : 0;
    const n = pr.normal;
    return vadd(vscale(n, k * drop * n[1]), vscale(pr.lat, -lat * drop));
  };
}
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
  PSAX: { probe: () => psaxProbe(...PSAX_MID), depth: 13, track: papTrack(() => psaxProbe(...PSAX_MID), 0.3, 0.25) },   // mid-papillary
  // (14.5 cm deep: the PA bifurcation and both branch origins and the LA body
  // behind the root; the sector's far-left corner, low on the patient's right,
  // stops short of the liver dome)
  PSAX_AV: { probe: psaxAvProbe, depth: 14.5, track: (A) => [0, 0.8 * (A.lv.A[1] - A0.lv.A[1]), 0] },   // follows the root's descent
  A4C: { probe: a4cProbe, depth: 17, track: a4cTrack },
  // (15 cm deep: the LA and its roof; the PA bifurcation and the right PA lie
  // just beyond the LA roof in this plane and run along the beam to 17 cm)
  // (review 12: 15.5 cm so the LA roof, 14.6 cm from the apex, lies inside the
  // sector with a 0.9 cm margin for the biplane LA trace)
  A2C: { probe: a2cProbe, depth: 15.5, track: APICAL_TRACK },
  SUBCOSTAL: { probe: subcostalProbe, depth: 18 },
};
// RV inflow (parasternal, tilted from the PLAX toward the right hip): the RV
// inflow in the near field, the tricuspid (anterior and posterior leaflets)
// hinged at its annulus, and the RA in the far field with the coronary sinus and
// the IVC (Eustachian valve) entering it.
// The plane was chosen by search (as for the ME LAA) over rolls about the
// tricuspid inflow axis and tilts off it: rolled RVIT_ROLL deg from the CS
// ostium direction, tilted RVIT_TILT deg off the axis, and set RVIT_CS_W of the
// way from the IVC orifice to the CS ostium. It cuts no septum or LV (the old
// plane through the septal-hugging ostia cut 13-19 % septum) and stays >= 0.15 cm
// off the convex septum through the cycle (rolled 30 / tilted 12 / CS weight 0.7
// it lay tangent to it at end-diastole, and the grazed septum and the anterior
// papillary muscle read as a free 3 cm2 RV mass); it keeps the CS ostium within
// ~0.4 cm, cuts the Eustachian ridge between the CS and the IVC, and leaves no
// great artery in the far field (rolled 35 / tilted 15 it cuts ~9 cm2 of the
// aortic root). The IVC, near vertical, crosses this tilted long-axis plane at
// ~40 deg (as it crosses the PLAX): its mouth is cut, not a length of tube.
// Within it the beam runs posteriorly from the chest wall, swung RVIT_SWING
// toward the inflow axis (the probe angled inferomedially toward the right hip),
// so the beam crosses the valve from the RV into the RA, meeting the tricuspid
// inflow within ~14 deg (the default tricuspid-inflow PW trace).
// RV apex to the left of the screen, base (RA) to the right, as in the PLAX.
// (15.5 cm deep: the whole RA, up to its roof, inside the sector.)
const RVIT_ROLL = 26, RVIT_TILT = 14, RVIT_SWING = 40, RVIT_CS_W = 0.6;    // deg, deg, deg, fraction
const csOstium = (A) => A.cs[A.cs.length - 1].b;
let _rvit = null;
function rvInflowProbe() {
  if (!_rvit) {
    const tv = VALVE_DEFS.tricuspid, ax = tv.n, T = tv.c;
    const cs = csOstium(A0), iv = A0.ivc.a;
    const u0 = inPlane(vsub(cs, T), ax), w0 = vcross(ax, u0), ph = RVIT_ROLL * Math.PI / 180;
    const u = norm(vadd(vscale(u0, Math.cos(ph)), vscale(w0, Math.sin(ph))));
    const n = vrot(norm(vcross(ax, u)), u, RVIT_TILT * Math.PI / 180);
    const tgt = vadd(T, vscale(n, vdot(vsub(cs, T), n) * RVIT_CS_W + vdot(vsub(iv, T), n) * (1 - RVIT_CS_W)));
    const post = inPlane(BODY_AX.P, n), ax2 = inPlane(ax, n), toRA = inPlane(vsub(ax2, vscale(post, vdot(ax2, post))), n);
    const sw = RVIT_SWING * Math.PI / 180;
    const dir = norm(vadd(vscale(post, Math.cos(sw)), vscale(toRA, Math.sin(sw))));
    let right = vcross(n, dir);
    if (right[1] < 0) right = vscale(right, -1);                 // base / RA to the right
    _rvit = chestWallProbe(tgt, dir, right, 8.0, WALL);
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
// In this heart the extension of the IVC axis runs within ~1.2 cm of the aortic
// valve (the IVC is aimed through the fossa at the upper septum, which the
// non-coronary sinus abuts), so a plane holding the whole IVC axis also cuts the
// root, the LA and the RVOT. The view is therefore framed on the hepatic IVC:
// aimed SC_IVC_AIM cm below the RA junction, tipped SC_IVC_CRANIAL and rolled
// SC_IVC_ROLL about the IVC axis, slid SC_IVC_SHIFT toward the patient's right
// and turned SC_IVC_YAW so that, while the hepatic IVC lies along the plane, the
// plane passes to the right of the axis at the junction: the IVC opens into the
// lateral RA body at the right of the sector and the root, the LA, the RVOT/PA
// and the abdominal aorta (all to the left) drop out. (Aimed 3.8 cm below the
// junction and 12.5 cm deep, so the IVC's widening into the RA and the RA
// chamber itself lie inside the sector, ~2.5 cm from its edge, not cut by it.)
const SC_IVC_CRANIAL = -5;                        // deg
const SC_IVC_ROLL = -5;                           // deg
const SC_IVC_AIM = 3.8, SC_IVC_SHIFT = 0.0;       // cm
const SC_IVC_YAW = 6;                             // deg
function subcostalIvcProbe() {
  // transducer at the subxiphoid window, aimed at the hepatic IVC; the plane
  // contains the IVC axis so the vessel runs across the image through the liver
  // into the RA (cranial right)
  const iv = A0.ivc;
  const ax = norm(vsub(iv.a, iv.b));                      // IVC axis, toward the RA
  // sagittal sub-xiphoid probe held back, so the IVC runs ACROSS the image —
  // caudal (liver) screen-left, entering the RA screen-right — rather than
  // diving into the far field, with the abdominal wall and then the liver in the
  // near field
  const tgt = vadd(iv.a, vscale(ax, -SC_IVC_AIM));
  const back = inPlane(vadd(BODY_AX.P, vscale(BODY_AX.S, 0.3)), vcross(ax, BODY_AX.P));
  let perp = norm(vsub(back, vscale(ax, vdot(back, ax))));
  // rolled about the IVC axis toward the patient's right
  const roll = SC_IVC_ROLL * Math.PI / 180;
  perp = norm(vadd(vscale(perp, Math.cos(roll)), vscale(vcross(ax, perp), Math.sin(roll))));
  const th = SC_IVC_CRANIAL * Math.PI / 180;
  const dir = norm(vadd(vscale(perp, Math.cos(th)), vscale(ax, Math.sin(th))));
  let right = norm(vcross(dir, ax));
  if (vdot(right, BODY_AX.L) > 0) right = vscale(right, -1);    // the patient's right
  const yw = SC_IVC_YAW * Math.PI / 180;
  const lat = norm(vadd(vscale(ax, Math.cos(yw)), vscale(right, Math.sin(yw))));
  return chestWallProbe(vadd(tgt, vscale(right, SC_IVC_SHIFT)), dir, lat, 9.0, 6.0);
}
// Suprasternal notch: the aortic arch in long axis, the ascending aorta on the
// left of the screen and the descending aorta on the right, with the right
// pulmonary artery in cross-section beneath the arch. The arch is planar (ARCH
// in anatomy.js), so the plane through the top of the ascending aorta, the arch
// and DTA_TOP is the arch's own plane: the candy cane opens below the summit,
// and beyond the descending limb lies the left lung, which wraps the descending
// aorta posterolaterally (edge-on to the beam, so it shows no pleural line).
function suprasternalProbe() {
  const asc = LM.AO_ASC_TOP, arch = LM.AO_ARCH_MID, desc = LM.DTA_TOP;
  const n = norm(vcross(vsub(asc, arch), vsub(desc, arch)));
  const tgt = vscale(vadd(vadd(asc, desc), vscale(arch, 2)), 0.25);
  const dir = inPlane(vscale(BODY_AX.S, -1), n);
  return probeFrom(tgt, dir, vsub(desc, asc), 7.0);
}
// PSAX-MV: the mitral leaflet level, 1 cm below the annulus at end-diastole.
// The tricuspid annulus lies apical to the mitral (0.85 cm at its centre at ED,
// 1.2 cm at ES as its free-wall side descends) and tilts down toward the RV
// free wall, so a plane exactly perpendicular to the LV long axis at this level
// passes above the tricuspid orifice on the left of the sector: from the
// mid-parasternal position (probe 0.7 cm lateral, no roll) it shows 2-4 cm2 of
// RA floor beside the inferoseptal wall. The plane is therefore rolled 8 deg
// about the beam (its right-heart edge toward the apex, its lateral edge still
// below the LA: the small angulation a sonographer makes to keep the RV, not
// the RA, beside the septum), with the beam turned 40 deg toward body-posterior
// as in the papillary-level short axis. It follows PSAX_MV_TRACK of the base's
// systolic descent, so at ES it sits ~0.55 cm below the annulus and still cuts
// the closed leaflets.
// At end-systole those two constraints leave a window of only ~3 mm (a plane
// high enough for the closed mitral leaflets crosses the posterior tricuspid
// annulus where it reaches 6 cm to the right of the LV centre; more roll brings
// the lateral LA in), so the probe is still slid 1.6 cm along the chest wall
// toward the LV lateral wall (it was 2.5 cm): the LV centre lies ~12 deg left
// of the beam axis (it was 19 deg), no RA or LA lumen is in the plane at any
// phase, and of the tricuspid valve only the leaflets hanging into the RV
// inflow show (<= 0.45 cm2).
const PSAX_MV = [-1.0, 1.6, 40, 0, -8];              // y, sx (cm), rot, tilt, roll (deg)
const PSAX_MV_TRACK = 0.7;
const psaxMvProbe = () => psaxProbe(...PSAX_MV);
function psaxMvTrack(A) {
  if (!A.msTip) return [0, PSAX_MV_TRACK * (A.lv.base - A0.lv.base), 0];
  // rheumatic MS: the plane slid along the long axis until it passes just above the
  // funnel tip (A.msTip), where the orifice is planimetered
  const p = psaxMvProbe(), q = vsub(vadd(A.msTip, [0, 0.04, 0]), p.pos);
  return [0, vdot(q, p.normal) / p.normal[1], 0];
}
// additional apical / parasternal / subcostal / suprasternal views
export const EXTRA_VIEWS = {
  A5C: { probe: a5cProbe, depth: 14, track: a5cTrack },
  A3C: { probe: a3cProbe, depth: 15, track: APICAL_TRACK },
  // mitral leaflet level, 1 cm below the annulus (closed coaptation line at ED,
  // fish-mouth in diastole): the probe follows the posterolateral annular
  // excursion (the small tilt a sonographer makes), so the plane stays at the
  // leaflets and never climbs into the posterior AV groove (CS / RA) in systole
  // In rheumatic MS the plane sits at the funnel tip (A.msTip, the open fish-mouth
  // orifice, ~6 deg off perpendicular to the inflow axis) — where the orifice is
  // planimetered — instead of at the leaflet bodies, which cut the dome obliquely.
  PSAX_MV: { probe: psaxMvProbe, depth: 13, track: psaxMvTrack },
  RVIT: { probe: rvInflowProbe, depth: 15.5, track: rvInflowTrack },
  SC_IVC: { probe: subcostalIvcProbe, depth: 12.5 },
  SSN: { probe: suprasternalProbe, depth: 16 },
};

// ---- TEE -------------------------------------------------------------------
// Multiplane TEE probes image a 90 deg sector (the TTE phased arrays here +/-0.66 rad)
export const TEE_HALF = Math.PI / 4, TTE_HALF = 0.66;
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
// The TEE probe travels one path. TEE_PATH_KNOTS: [s, r] — s cm along the
// body's cranial axis from the LA-level reference ESO, r cm toward the patient's
// right (taken perpendicular to the cranial axis), joined by a smooth monotone
// cubic (no corners, so the beam carried along it on advance turns gradually).
// Behind the atria the oesophagus lies against the posterior LA wall, behind the
// right half of the LA at the level of the mitral annulus and the roof, by the
// right pulmonary-vein ostia (r 2.2: ~1.4 cm right of the body midline, LA blood
// ~1.1 cm in front; the mid-oesophageal station below), right of and in front of
// the descending aorta, which runs posterolateral to it on the left; at the level
// of the arch it lies to the right of the aortic isthmus (never inside the aorta),
// and below the atria it swings left and forward to the hiatus. (The oesophagus
// runs anywhere from the left to the right pulmonary veins behind the LA; this
// course is where the LV long axis, continued from the apex through the mitral
// valve, meets the posterior LA wall (~1.9 cm right of the LA centre), and the one
// from which the four-chamber sector holds the interatrial septum and a secundum
// defect with some margin: with the station at r 1.5 / 1.8 / 2.0 / 2.2 the fossa
// lies 40 / 38 / 37 / 35 deg left of the beam, against the 45 deg sector edge,
// and the defect's LA-RA lumen contact in the ME4C reads 0.2 / 0.4 / 0.6 / 0.8 cm.)
// Below the gastro-oesophageal junction (s = TEE_GEJ) the path
// enters the stomach: the fundus contact under the inferior wall (the transgastric
// views) and, further on toward the apex, the deep transgastric position.
const TEE_PATH_KNOTS = [[-6, 0], [-4, 0.6], [-2, 1.7], [-0.5, 2.2], [2.5, 2.2], [4.5, 1.6], [7.0, 1.4], [9.0, 0.6], [11, 0]];
export const TEE_GEJ = -6;
const _S = BODY_AX.S, _ER = norm(vsub(BODY_AX.R, vscale(_S, vdot(BODY_AX.R, _S))));
// Fritsch-Butland slopes of the monotone cubic through the knots (flat at both ends)
const _KM = (() => {
  const K = TEE_PATH_KNOTS, n = K.length, d = [], m = new Array(n).fill(0);
  for (let i = 0; i + 1 < n; i++) d.push((K[i + 1][1] - K[i][1]) / (K[i + 1][0] - K[i][0]));
  for (let i = 1; i + 1 < n; i++) {
    const h0 = K[i][0] - K[i - 1][0], h1 = K[i + 1][0] - K[i][0];
    if (d[i - 1] * d[i] > 0) m[i] = 3 * (h0 + h1) / ((2 * h1 + h0) / d[i - 1] + (h1 + 2 * h0) / d[i]);
  }
  return m;
})();
function teePathR(s) {
  const K = TEE_PATH_KNOTS;
  if (s <= K[0][0]) return K[0][1];
  for (let i = 0; i + 1 < K.length; i++) {
    if (s > K[i + 1][0]) continue;
    const h = K[i + 1][0] - K[i][0], t = (s - K[i][0]) / h, t2 = t * t, t3 = t2 * t;
    return (2 * t3 - 3 * t2 + 1) * K[i][1] + (t3 - 2 * t2 + t) * h * _KM[i] + (-2 * t3 + 3 * t2) * K[i + 1][1] + (t3 - t2) * h * _KM[i + 1];
  }
  return K[K.length - 1][1];
}
export function teePathPoint(s) {
  return vadd(vadd(ESO, vscale(_S, s)), vscale(_ER, teePathR(s)));
}
// The mid-oesophageal station: ONE probe position and beam from which the ME
// four-chamber, commissural, two-chamber and long-axis views are reached by
// turning the omniplane alone, as on a real multiplane probe. The transducer sits
// on the oesophageal path ME_S cm above the LA-level reference, behind the right
// half of the LA where the path runs ~1.4 cm right of the body midline, and
// the beam is aimed ME_AIM of the way from the mitral-valve centre to the
// epicardial apex, ME_POST cm posterior (heart -z, off the LVOT: the transducer any
// higher, or the aim any more anterior, lays the aortic root wall in the
// four-chamber plane), so it runs within a few degrees of the LV long axis: the
// apex then stays within ~0.5 cm of the imaging plane at every omniplane angle
// (rotation turns the plane about the beam, and the beam runs down the LV) — the
// mitral valve centred and the LV apex in the far field in each view. From this
// station (rather than ~2.2 cm to the left, behind the LA centre) the
// interatrial septum, the fossa and a secundum defect lie inside the 90 deg
// four-chamber sector (~35 deg left of the beam), the 2C plane at +90 deg holds
// the appendage beside the anterior mitral annulus, and the ME LAX holds the
// LVOT, cusps and root ~29 deg right of the beam.
const ME_S = 1.0, ME_AIM = 0.65, ME_POST = 0.4;
function meBaseProbe() {
  const pos = teePathPoint(ME_S), tgt = vadd(vadd(M, vscale(vsub(APEX, M), ME_AIM)), [0, 0, -ME_POST]);
  const dir = norm(vsub(tgt, pos));
  const { lft } = teeFrame(dir);                     // 0 deg: the patient's left on the screen-right
  return teeClear(probeFrom(tgt, dir, lft, Math.hypot(...vsub(tgt, pos))));
}
// omniplane angle (deg, searched in [lo, hi]) of the ME station that lays point q in the plane
function meOmniThrough(q, lo, hi) {
  const b = meBaseProbe();
  let best = null;
  for (let o = lo; o <= hi + 1e-9; o += 0.25) {
    const p = teeManipulate(b, { omni: o }), d = Math.abs(vdot(vsub(q, p.pos), p.normal));
    if (!best || d < best.d) best = { o, d };
  }
  return best.o;
}
// ME4C (~0-10 deg): ME4C_TRIM beyond the omniplane angle that brings the
// tricuspid-valve centre into the plane (the mitral valve is on the beam): both AV
// valves, the crux and the atrial septum in profile, both ventricles, the LV apex
// in the far field (the trim keeps the RV apex short of the LV's and the plane off
// the aortic root through systole: from the station behind the right half of the
// LA, 1.75 deg leaves the RV lumen ending ~0.9 cm short of the LV apex at end-
// systole, where 0.75 deg left it 0.82 cm short, against the audit's 0.8).
// MECOMM = ME4C + 65 with a small turn (P3-A2-P1 and the papillary muscles; below), ME2C = ME4C + 90
// (the inferior wall and CS on the left, the anterior wall and LAA on the right),
// MELAX = the angle that brings the aortic-valve centre into the plane (~150 deg:
// the LVOT, cusps and root to the screen-right, the inferolateral wall on the left).
const ME4C_TRIM = 1.75;                                // deg
// The commissural preset is the omniplane at ME4C + MECOMM_OMNI with a small turn
// of the shaft, MECOMM_TURN (as the view is found at the bedside: omniplane to
// ~60-70 deg, then the probe turned a few degrees until the muscles sit on the
// walls). The transducer does not move. With the beam left on the ME4C line the
// + 60 deg plane holds both papillary heads (0.02 / 0.06 cm off) while their
// bases lie 0.5-0.8 cm to either side of it, so it runs up the free flank of
// both muscles: two free columns with a channel behind each at end-diastole and
// a slit cavity with a detached pocket at end-systole. Turned, the plane still
// crosses both commissures (it moves ~0.4 cm at the annulus), passes through
// the wall-attached base of the posteromedial muscle and just clear of the
// anterolateral head, and the cavity is one piece at both phases; the turn
// carries the plane toward the inferior RV, which the extra 5 deg of omniplane
// takes back out. Scan (audit rows: RV area cm2 / lower-70 % cleft / double
// barrel / island cm2 / one piece ED / ES):
//   + 60, no turn      0.49 / 0.85 / 2.15 / 2.74 / 1.00 / 0.83
//   + 44..56, no turn  ... / 0.50-0.60 / 1.50-2.05 / 0 / 1.00 / 0.92-0.93
//   + 60, 4.5 deg      1.53 / 0.65 / 0.85 / 0 / 1.00 / 0.97
//   + 63, 4.5 deg      1.18 / 0.60 / 0.95 / 0 / 1.00 / 0.970
//   + 65, 4.5 deg      0.95 / 0.55 / 1.00 / 0 / 1.00 / 0.973   (chosen)
//   + 65, 4 or 5 deg   0.89, 1.07 / 0.55, 0.65 / 1.15, 0.85 / 0 / 0.95-0.96 / 0.97
//   + 66, 4.5 deg      0.87 / 0.55 / 1.05 / 0 / 0.957 / 0.973
// (6 deg of turn loses the anterolateral muscle from the view altogether.)
const MECOMM_OMNI = 65, MECOMM_TURN = 4.5;            // deg
let _meOm = null;
const meOmni = () => _meOm || (_meOm = (() => {
  const o4 = meOmniThrough(T, -15, 40) + ME4C_TRIM;
  return { ME4C: o4, MECOMM: o4 + MECOMM_OMNI, ME2C: o4 + 90, MELAX: meOmniThrough(A, 110, 175) };
})());
const meView = (key) => { let c = null; return () => { if (!c) c = teeManipulate(key === 'MECOMM' ? teeManipulate(meBaseProbe(), { rotate: MECOMM_TURN }) : meBaseProbe(), { omni: meOmni()[key] }); return { ...c, pos: c.pos.slice(), dir: c.dir.slice(), lat: c.lat.slice(), normal: c.normal.slice(), target: c.target.slice() }; }; };
const me4cProbe = meView('ME4C');
const melaxProbe = meView('MELAX');
// The transducer is anteflexed and slid TG_ROLL further round the LV long axis
// toward the diaphragmatic face than the reference gastric contact, along the
// greater curvature, so it looks up through the INFERIOR wall: the posteromedial
// papillary muscle then sits near the top of the sector (near field, ~11 o'clock)
// and the anterolateral one at the right (~2.5 o'clock), the ASE transgastric
// orientation used to name the segments and their coronary territories. It rests
// where that beam line meets the gastric lumen. (The slide is no more than the
// fundus allows: a few degrees further and the beam line misses the lumen, leaving
// the transducer in the tissue beyond it. Steering the beam off the LV centre from
// the contact instead would push the RV off the side of the sector.)
const TG_ROLL = 6;                                   // deg
function tgsaxProbe() {
  // transgastric: probe in the fundus below the heart, looking up through the
  // inferior wall; mid-papillary short axis, lateral wall on the right
  // (the transducer rests on the gastric wall; anatomy.js draws that wall in the
  // first ~0.45 cm, then the diaphragm and the inferior wall)
  const G = LM.GASTRIC, c = vadd(G.p, vscale(G.n, 4.6));
  const n = vrot(G.n, Y, TG_ROLL * Math.PI / 180);
  let t = 4.6;
  while (t < 8 && stomachDist(...vsub(c, vscale(n, t))) > 0) t += 0.02;
  return probeFrom(c, n, [1, 0, 0], t);
}
// ME LAA (~60-90 deg, turned left): the appendage along its length as a finger
// on the screen-right of the LA, the LSPV beside it across the warfarin ridge.
// The plane holds the oesophagus and the LAA ostium; its roll about that line is
// the one that lays the most appendage lumen in plane while keeping the LSPV
// ostium within 0.4 cm of it and the least ventricular myocardium in the sector
// at end-diastole and end-systole (the LV wall swings into a plane grazing the
// AV groove as the base descends).
const A_ES = anatomyParams(1, 0, {}, null, 0.42);
let _melaa = null;
// the appendage rides down with the AV plane in systole: the operator follows it
// (a slight probe flex), so the plane shifts along its normal with the LAA
function melaaTrack(A) {
  const p = melaaProbe(), n = p.normal;
  const c = vscale(vadd(A.la.aa[0].b, A.pv.veins[0].w), 0.5), c0 = vscale(vadd(A0.la.aa[0].b, A0.pv.veins[0].w), 0.5);
  return vscale(n, vdot(vsub(c, c0), n));
}
// The plane is then chosen among beams turned up to MELAA_TURN / MELAA_TILT about
// the line to the appendage neck and rolled about the beam: the one that lays the
// most appendage (beyond the LA body), LA and LSPV lumen in the sector at
// end-diastole and end-systole (the plane following the appendage, melaaTrack),
// with the least RV, aortic root, myocardium (LV wall above 1.5 cm2 heavily) and
// empty space, while the LSPV ostium and the warfarin ridge stay within 0.4 cm
// of the plane: the LUPV and the windsock appendage then both run to the right
// of the LA, the ridge between them.
// (The plane is chosen from the LA-level reference behind the appendage, and the
// transducer then slides within that plane to where it crosses the oesophageal
// path, ~1 cm lower and ~2 cm to the right: the same plane and beam direction,
// seen from the probe's place in the oesophagus.)
const MELAA_TURN = [-8, 0, 8], MELAA_TILT = [-16, -8, 0, 8, 16], MELAA_DEPTH = 8;
function melaaScore(p, A) {
  const c = { laa: 0, la: 0, pv: 0, rv: 0, ao: 0, myo: 0, lvm: 0, out: 0 }, h = 0.3;
  for (let d = h; d < MELAA_DEPTH; d += h) for (let x = -d * 0.75; x <= d * 0.75; x += h) {
    if (Math.hypot(x, d) > MELAA_DEPTH) continue;
    const q = vadd(p.pos, vadd(vscale(p.dir, d), vscale(p.lat, x)));
    const k = bodyClassify(q[0], q[1], q[2], A).code;
    if (k === BODY.LA) {
      if (lumenDist(...q, A, 'LAA') < 0 && lumenDist(...q, A, 'LA_BODY') > 0) c.laa++;
      else if (lumenDist(...q, A, 'PV') < 0) c.pv++;
      else c.la++;
    } else if (k === BODY.RV) c.rv++;
    else if (k === BODY.AO) c.ao++;
    else if (k === BODY.MYO) { c.myo++; if (lumenDist(...q, A, 'LV_EPI') < 0) c.lvm++; }
    else if (k === BODY.OUT) c.out++;
  }
  for (const key in c) c[key] *= h * h;
  return c;
}
function melaaProbe() {
  if (_melaa) return { ..._melaa, pos: _melaa.pos.slice(), target: _melaa.target.slice() };
  const aa = A0.la.aa, v = A0.pv.veins[0], tip = aa[2].b;
  const rg = vscale(vadd(A0.pv.ridge.a, A0.pv.ridge.b), 0.5);
  const tgt0 = vscale(vadd(aa[0].b, v.w), 0.5);
  const d0 = norm(vsub(tgt0, ESO));
  const e1 = norm(vcross(d0, Y)), e2 = vcross(e1, d0);
  const lead = (A) => vscale(vadd(A.la.aa[0].b, A.pv.veins[0].w), 0.5);
  let pick = null;
  for (const ta of MELAA_TURN) for (const tb of MELAA_TILT) {
    const dir = vrot(vrot(d0, e1, ta * Math.PI / 180), e2, tb * Math.PI / 180);
    const f1 = norm(vcross(dir, Y)), f2 = vcross(f1, dir);
    for (let roll = 0; roll < 180; roll += 10) {
      let lat = vadd(vscale(f1, Math.cos(roll * Math.PI / 180)), vscale(f2, Math.sin(roll * Math.PI / 180)));
      const n = vcross(dir, lat);
      if (Math.abs(vdot(vsub(v.a, ESO), n)) > 0.4 || Math.abs(vdot(vsub(rg, ESO), n)) > 0.4) continue;
      if (vdot(lat, vsub(tip, aa[0].a)) < 0) lat = vscale(lat, -1);       // appendage tip on screen-right
      const p = probeFrom(vadd(ESO, vscale(dir, 3)), dir, lat, 3);
      const pES = { ...p, pos: vadd(p.pos, vscale(n, vdot(vsub(lead(A_ES), lead(A0)), n))) };
      const s = [melaaScore(p, A0), melaaScore(pES, A_ES)];
      const lo = (k) => Math.min(s[0][k], s[1][k]), hi = (k) => Math.max(s[0][k], s[1][k]);
      const score = 3 * Math.min(4, lo('laa')) + 0.5 * Math.min(6, lo('la')) + 1.5 * Math.min(1.5, lo('pv'))
        - hi('rv') - 0.6 * hi('ao') - 0.2 * hi('myo') - 2 * Math.max(0, hi('lvm') - 1.5) - 0.15 * hi('out');
      if (!pick || score > pick.score) pick = { score, dir, lat };
    }
  }
  const pos = pathInPlane(ESO, vcross(pick.dir, pick.lat), -4, 3), sd = vdot(vsub(tgt0, pos), pick.dir);
  _melaa = teeClear(probeFrom(vadd(pos, vscale(pick.dir, sd)), pick.dir, pick.lat, sd));
  return melaaProbe();
}
// The point where the oesophageal path (s in [s0, s1]) crosses the plane through
// q with normal n (projected into it): a plane chosen from the anatomy is imaged
// from the probe's place in the oesophagus.
function pathInPlane(q, n, s0, s1) {
  let best = null;
  for (let s = s0; s <= s1 + 1e-9; s += 0.01) {
    const p = teePathPoint(s), d = vdot(vsub(p, q), n);
    if (!best || Math.abs(d) < Math.abs(best.d)) best = { p, d };
  }
  return vsub(best.p, vscale(n, best.d));
}
// ME two-chamber (ME4C + 90 deg): LA near field, LV below, the LAA and anterior
// wall on the right, the inferior wall and coronary sinus on the left.
const me2cProbe = meView('ME2C');
// ME aortic-valve short axis (~30-60 deg): the three cusps en face, the LA in
// the near field, the RVOT in the far field; NCC beside the interatrial septum.
// The transducer stays on the oesophageal path, MEAVSAX_S cm above the LA-level
// reference (0.25 cm below the ME station: the probe barely moves between the ME
// views and the AV SAX), and the plane holds the transducer and the cusps (the
// beam aimed at the cusp level) while turned as near to perpendicular to the root
// axis as that allows: from the path behind the right half of the LA, within ~3
// deg of a true root short axis (from behind the LA centre it was ~17 deg off). The
// tricuspid valve then lies within ~0.3-0.6 cm of the plane, so TR is recorded by
// CW here (aligned, as from the ME RV inflow-outflow), and the secundum ASD PW
// meets the shunt at ~12 deg.
const MEAVSAX_S = 0.75, MEAVSAX_FOSSA = 0.7;
function meAvSaxProbe() {
  const tgt = vadd(A, vscale(LM.U_AO, 0.75)), pos = teePathPoint(MEAVSAX_S);
  const dir = norm(vsub(tgt, pos));
  const n0 = norm(vsub(LM.U_AO, vscale(dir, vdot(LM.U_AO, dir))));
  let l0 = vcross(n0, dir);
  if (vdot(l0, LM.E_SCR) < 0) l0 = vscale(l0, -1);
  // rolled about the beam no further than needed to keep the centre of the fossa ovalis (the
  // interatrial septum beside the non-coronary cusp) within MEAVSAX_FOSSA of the plane
  let lat = l0;
  for (let r = 0; r <= 20; r += 0.5) {
    const cand = [r, -r].map((d) => vrot(l0, dir, d * Math.PI / 180))
      .filter((l) => Math.abs(vdot(vsub(LM.IAS.asdC, pos), vcross(dir, l))) <= MEAVSAX_FOSSA);
    if (cand.length) { lat = cand[0]; break; }
  }
  return teeClear(probeFrom(tgt, dir, lat, Math.hypot(...vsub(tgt, pos))));
}
// ME bicaval (~90-110 deg): the RA with the SVC entering on the right of the
// screen and the IVC on the left, the interatrial septum and fossa in profile
// between the LA (near field) and the RA — the view for PFO / sinus venosus ASD.
// The transducer is withdrawn BICAVAL_UP cm from the LA-level reference and offset
// BICAVAL_OFF (cm, posterior and toward the patient's left), where it was before
// the oesophageal path moved to the right half of the LA (review 12): ~2.4 cm left
// of the moved path. It is kept there with the aim below (the septal / caval
// group's plane): the plane through the SVC and IVC axis points and the fossa is
// near-coronal and never reaches the moved oesophagus, so from any point on it
// (GEJ to 4 cm above the LA reference) the plane through the two caval points
// passes 0.9-1.1 cm from the fossa, and planes refitted through the path point,
// the fossa and the caval axes cut the cavae for < 1.5 cm or lose the crista, the
// Eustachian valve or the LA near field (see the TEE path rows). The plane holds
// the transducer, the SVC axis BICAVAL_T[0] cm above its orifice and the IVC axis
// BICAVAL_T[1] cm below its orifice: both cavae then run in from the sides as
// tubes (cut along ~2 cm each), the fossa lies within 0.45 cm of the plane, the
// septum crosses the image for 3+ cm and the beam passes through the middle of
// the LA body (not its septal edge), so the LA is a broad, echo-free near field.
// The beam is aimed BICAVAL_AIM of the way from the IVC point to the SVC point.
const BICAVAL_UP = 2.5;                          // cm withdrawn along the oesophagus
const BICAVAL_OFF = [-0.3, -0.3];                 // cm: body-posterior, patient's-right share (within the lumen)
const BICAVAL_T = [1.5, 0.75];                   // cm along the SVC / IVC axes from their orifices
const BICAVAL_AIM = 0.5;
let _bicaval = null;
function meBicavalProbe() {
  if (!_bicaval) {
    const S = BODY_AX.S, e1 = norm(vsub(BODY_AX.P, vscale(S, vdot(BODY_AX.P, S)))), e2 = vcross(S, e1);
    const pos = vadd(vadd(vadd(ESO, vscale(S, BICAVAL_UP)), vscale(e1, BICAVAL_OFF[0])), vscale(e2, BICAVAL_OFF[1]));
    const sv = A0.svc, iv = A0.ivc;
    const P1 = vadd(sv.a, vscale(norm(vsub(sv.b, sv.a)), BICAVAL_T[0]));
    const P2 = vadd(iv.a, vscale(norm(vsub(iv.b, iv.a)), BICAVAL_T[1]));
    const tgt = vadd(P2, vscale(vsub(P1, P2), BICAVAL_AIM));
    const dir = norm(vsub(tgt, pos));
    // SVC (cranial) on the screen-right
    _bicaval = teeClear(probeFrom(tgt, dir, vsub(P1, P2), Math.hypot(...vsub(tgt, pos))));
  }
  return { ..._bicaval, pos: _bicaval.pos.slice(), target: _bicaval.target.slice() };
}
// ME RV inflow-outflow (~60-75 deg): RA and tricuspid on the left, the RV
// wrapping round the aortic valve to the RVOT and pulmonary valve on the right.
const MERVIO_S = 0;                                  // cm up the oesophageal path
const MERVIO_AV_H = 0.3;                             // cm above the aortic annulus: under the cusps' coaptation at ED
function meRvioProbe() {
  // the plane from the oesophagus through the aortic valve at cusp level that
  // passes closest to the tricuspid and the pulmonary valve centres: the RV wraps
  // from its inflow (screen-left) round the aortic valve to the outflow, and the
  // shut aortic cusps make their 'Y' in the middle of the sector. (The plane
  // through the oesophagus and the two right-sided valve centres alone cuts the
  // root ~0.9 cm above the annulus from the mid-oesophageal station, above the
  // coaptation of one cusp pair (that arm of the 'Y' is missing at 0.9 cm, faint
  // at 0.45-0.6 cm); holding the cusp level costs 0.7-0.8 cm of off-centre cut
  // on each of the two valves, inside their annuli.)
  const ESO = teePathPoint(MERVIO_S);
  const n0 = norm(vcross(vsub(T, ESO), vsub(LM.PV, ESO)));
  const u = norm(vsub(vadd(A, vscale(LM.U_AO, MERVIO_AV_H)), ESO));
  const e1 = norm(vcross(u, n0)), e2 = vcross(u, e1);
  // n = e1 sin(th) + e2 cos(th), perpendicular to u, minimising the squared
  // distances of the two valve centres from the plane
  const a = [vdot(vsub(T, ESO), e1), vdot(vsub(T, ESO), e2)], b = [vdot(vsub(LM.PV, ESO), e1), vdot(vsub(LM.PV, ESO), e2)];
  const mxx = a[0] * a[0] + b[0] * b[0], myy = a[1] * a[1] + b[1] * b[1], mxy = a[0] * a[1] + b[0] * b[1];
  const th = 0.5 * Math.atan2(2 * mxy, mxx - myy) + Math.PI / 2;      // minor axis of the 2 x 2 form
  let n = norm(vadd(vscale(e1, Math.cos(th)), vscale(e2, Math.sin(th))));
  if (vdot(n, n0) < 0) n = vscale(n, -1);
  const tgt = vsub(A, vscale(n, vdot(vsub(A, ESO), n)));
  const dir = norm(vsub(tgt, ESO));
  let right = vsub(LM.PV, T);
  right = vsub(right, vscale(dir, vdot(right, dir)));
  right = vsub(right, vscale(n, vdot(right, n)));
  return teeClear(probeFrom(tgt, dir, right, vdot(vsub(tgt, ESO), dir)));
}
// the tricuspid annulus descends ~1-2 cm in systole and this plane's normal runs
// largely along the LV axis, so the plane follows the cardiac base along its
// normal (the small advance / flex that holds the TV and RV inflow in view): the
// tricuspid annulus mostly, the aortic root (which descends about half as far)
// for the rest — so at end-systole the root is still cut near its annulus, not
// through the LVOT below it
let _rvioN = null;
const RVIO_TV_W = 0.7;
function meRvioTrack(A) {
  if (!_rvioN) _rvioN = meRvioProbe().normal;
  const d = vadd(vscale(vsub(A.valves.tricuspid.c, A0.valves.tricuspid.c), RVIO_TV_W), vscale(vsub(A.valves.aortic.c, A0.valves.aortic.c), 1 - RVIO_TV_W));
  return vscale(_rvioN, vdot(d, _rvioN));
}
// Descending thoracic aorta short axis: the probe turned to face posteriorly;
// the round aorta lies immediately behind the oesophagus. The probe sits on the
// path DESCAO_S cm below the LA-level reference, in the lower oesophagus where it
// runs in front of and just right of the aorta on its way to the hiatus (behind
// the LA the path lies further right, behind the right half of the atrium),
// facing the aortic axis at the same level.
const DESCAO_GAP = 1.05;                             // probe face -> aortic lumen (cm)
const DESCAO_S = -4.5;
let _descao = null;
function descAoSaxProbe() {
  if (!_descao) _descao = descAoSaxBuild();
  return { ..._descao, pos: _descao.pos.slice(), dir: _descao.dir.slice(), lat: _descao.lat.slice(), normal: _descao.normal.slice(), target: _descao.target.slice() };
}
function descAoSaxBuild() {
  const ESO = teePathPoint(DESCAO_S);
  // (the descending aorta runs straight down the body axis below DTA_P)
  const tgt = vadd(LM.DTA_P, vscale(BODY_AX.S, vdot(vsub(ESO, LM.DTA_P), BODY_AX.S)));
  const raw = vsub(tgt, ESO);
  const dir = norm(vsub(raw, vscale(BODY_AX.S, vdot(raw, BODY_AX.S))));
  const pr = probeFrom(tgt, dir, BODY_AX.L, Math.max(1.0, vdot(raw, dir)));
  // facing the aorta, the LA lies behind and beside the transducer: close in on the
  // aortic wall (>= 0.5 cm stand-off) until no beam of the sector reaches LA blood
  const laInSector = () => {
    for (let a = -TEE_HALF; a <= TEE_HALF + 1e-3; a += TEE_HALF / 8) {
      const d = vadd(vscale(pr.dir, Math.cos(a)), vscale(pr.lat, Math.sin(a)));
      for (let r = 0; r <= 3.5; r += 0.1) if (lumenDist(...vadd(pr.pos, vscale(d, r)), A0, 'LA') < 0) return true;
    }
    return false;
  };
  // but never nearer than DESCAO_GAP: the oesophageal wall and periaortic tissue lie
  // between the probe face and the aorta, whose whole round section (wall
  // included) then fits inside the sector instead of being clipped by its apex
  const gap = () => vdot(vsub(tgt, pr.pos), pr.dir) - LM.DTA_R;
  for (let it = 0; it < 30 && laInSector() && gap() > DESCAO_GAP + 0.05; it++) {
    pr.pos = vadd(pr.pos, vscale(pr.dir, 0.05));
  }
  if (gap() < DESCAO_GAP) pr.pos = vsub(pr.pos, vscale(pr.dir, DESCAO_GAP - gap()));
  // the shaft turned DESCAO_TURN (to the patient's left) so the edge of the 90 deg
  // sector clears the left lower pulmonary vein, the aorta still whole in the sector
  return teeManipulate(pr, { rotate: DESCAO_TURN });
}
const DESCAO_TURN = -4;                              // deg
export const TEE_VIEWS = {
  ME4C: { probe: me4cProbe, depth: 16 },
  ME2C: { probe: me2cProbe, depth: 16 },
  MELAX: { probe: melaxProbe, depth: 15 },
  MEAVSAX: { probe: meAvSaxProbe, depth: 10, track: (A) => [0, 0.8 * (A.lv.A[1] - A0.lv.A[1]), 0] },   // follows the root
  MEBICAVAL: { probe: meBicavalProbe, depth: 11 },
  MERVIO: { probe: meRvioProbe, depth: 12, track: meRvioTrack },              // follows the tricuspid annulus
  DESCAO: { probe: descAoSaxProbe, depth: 8 },
  TGSAX: { probe: tgsaxProbe, depth: 12 },
  // Mid-oesophageal left-atrial-appendage view: the plane CONTAINS the appendage's
  // long axis, so the ostium, neck and windsock body lie in-plane beside the LUPV.
  MELAA: { probe: melaaProbe, depth: MELAA_DEPTH, track: melaaTrack },
};

// ---- one TEE probe: omniplane angle, advance / withdraw, rotation, flexion ----
// Every TEE view above is a preset of one probe on one path (teePathPoint, then
// the stomach). The probe can then be manipulated as on a real TEE: the
// omniplane angle turns the imaging plane about the beam (0 deg transverse with
// the patient's left on the screen-right, 90 deg longitudinal with cranial on
// the screen-right, 180 deg the mirror of 0), advance (+) / withdraw (-) moves
// the transducer along the one path (teePathAt: the oesophagus, then through the
// gastro-oesophageal junction along the gastric polyline to the transgastric and
// deep transgastric stations), the beam carried with the path tangent; rotation
// turns the shaft about the local path tangent (+ to the patient's right,
// clockwise) and flexion tilts the tip (+ anteflexion: the beam swings
// cranially — in the oesophagus from the ME4C up into the LVOT and aortic valve
// (ME 5C), in the stomach up into the heart; - retroflexion: caudally, from the
// ME4C down to the coronary sinus).
function teeFrame(dir) {
  const up = norm(vsub(_S, vscale(dir, vdot(_S, dir))));
  let lft = vcross(dir, up);
  if (vdot(lft, BODY_AX.L) < 0) lft = vscale(lft, -1);
  return { up, lft };
}
// omniplane angle (deg) of a probe's imaging plane
export function teeOmni(p) {
  const { up, lft } = teeFrame(p.dir);
  return Math.atan2(vdot(p.lat, up), vdot(p.lat, lft)) * 180 / Math.PI;
}
// The whole probe path as one coordinate u (cm): u >= TEE_GEJ is the oesophagus
// (u = s of teePathPoint, + toward the mouth); below the junction u runs on along
// the gastric polyline (GEJ -> TG SAX contact -> deep TG contact), u = TEE_GEJ
// minus the distance travelled into the stomach. Advancing the probe lowers u.
let _gastric = null;
function gastricPolyline() {
  if (!_gastric) _gastric = [teePathPoint(TEE_GEJ), tgsaxProbe().pos, dtg5cProbe().pos];
  return _gastric;
}
export const TEE_PATH_TOP = 12;                      // cm: the upper oesophagus, above the arch
let _pathMin = null;
export function teePathMin() {
  if (_pathMin == null) {
    const g = gastricPolyline();
    let L = 0;
    for (let i = 0; i + 1 < g.length; i++) L += Math.hypot(...vsub(g[i + 1], g[i]));
    _pathMin = TEE_GEJ - L;
  }
  return _pathMin;
}
export function teePathAt(u) {
  if (u >= TEE_GEJ) return teePathPoint(Math.min(u, TEE_PATH_TOP));
  const g = gastricPolyline();
  let L = TEE_GEJ - u;
  for (let i = 0; i + 1 < g.length; i++) {
    const seg = vsub(g[i + 1], g[i]), l = Math.hypot(...seg);
    if (L <= l || i + 2 === g.length) return vadd(g[i], vscale(seg, Math.min(L, l) / l));
    L -= l;
  }
  return g[g.length - 1].slice();
}
// unit tangent of the path at u, pointing toward the mouth (withdrawal); taken
// over +/-0.6 cm so the frame turns through the junction rather than jumping
export function teePathTangent(u) {
  return norm(vsub(teePathAt(Math.min(TEE_PATH_TOP, u + 0.6)), teePathAt(Math.max(teePathMin(), u - 0.6))));
}
// path coordinate of the path point nearest q
export function teePathParam(q) {
  const near = (u0, u1, du) => {
    let best = null;
    for (let u = Math.max(teePathMin(), u0); u <= Math.min(TEE_PATH_TOP, u1) + 1e-9; u += du) {
      const d = Math.hypot(...vsub(q, teePathAt(u)));
      if (!best || d < best.d) best = { u, d };
    }
    return best.u;
  };
  const u = near(teePathMin(), TEE_PATH_TOP, 0.5);
  return near(u - 0.5, u + 0.5, 0.02);
}
// points along the path from u1 down to u0 (the 3D oesophagus and probe shaft)
export function teePathPolyline(u0 = teePathMin(), u1 = TEE_PATH_TOP, step = 0.5) {
  const out = [];
  for (let u = u1; u > u0 - 1e-9; u -= step) out.push(teePathAt(u));
  return out;
}
const ROT_SIGN = vdot(vrot(BODY_AX.A, _S, 0.1), BODY_AX.R) > vdot(BODY_AX.A, BODY_AX.R) ? 1 : -1;
// the minimal rotation taking unit a onto unit b, applied to v
function turnWith(v, a, b) {
  const k = vcross(a, b), s = Math.hypot(...k), c = vdot(a, b);
  if (s < 1e-9) return v.slice();
  return vrot(v, k, Math.atan2(s, c));
}
export function teeManipulate(p, { omni = null, advance = 0, rotate = 0, flex = 0 } = {}) {
  const om = omni == null ? teeOmni(p) : omni;
  let pos = p.pos.slice(), dir = p.dir.slice();
  let target = (p.target || vadd(p.pos, vscale(p.dir, 5))).slice();
  if (!rotate && !flex && !advance && omni == null) return { ...p, pos, dir, lat: p.lat.slice(), normal: p.normal.slice(), target };
  let tan = _S;
  if (advance || rotate) {
    const u0 = teePathParam(pos);
    let u1 = u0;
    if (advance) {
      // the transducer moves along the path; the beam (and any offset of the
      // transducer from the path centreline) turns with the path tangent
      u1 = Math.max(teePathMin(), Math.min(TEE_PATH_TOP, u0 - advance));
      const t0 = teePathTangent(u0), t1 = teePathTangent(u1), q0 = teePathAt(u0);
      const rel = vsub(target, pos);
      pos = vadd(teePathAt(u1), turnWith(vsub(pos, q0), t0, t1));
      dir = norm(turnWith(dir, t0, t1)); target = vadd(pos, turnWith(rel, t0, t1));
    }
    tan = teePathTangent(u1);
  }
  if (rotate) {
    const a = ROT_SIGN * rotate * Math.PI / 180;
    dir = vrot(dir, tan, a); target = vadd(pos, vrot(vsub(target, pos), tan, a));
  }
  if (flex) {
    // about the axis that tips the beam toward the body's cranial axis: +flex
    // (anteflexion) raises dir.S, -flex (retroflexion) lowers it
    let k = vcross(dir, _S);
    if (Math.hypot(...k) < 0.15) k = p.normal;
    const a = flex * Math.PI / 180;
    dir = vrot(dir, k, a); target = vadd(pos, vrot(vsub(target, pos), k, a));
  }
  const { up, lft } = teeFrame(dir), t = om * Math.PI / 180;
  const lat = norm(vadd(vscale(lft, Math.cos(t)), vscale(up, Math.sin(t))));
  return { pos, dir, lat, normal: vcross(dir, lat), target };
}
// a TEE view as an omniplane preset of another view's probe position and beam
const omniOf = (probe, omni) => { let c = null; return () => { if (!c) c = teeManipulate(probe(), { omni }); return { ...c, pos: c.pos.slice(), target: c.target.slice() }; }; };

// ME aortic-valve long axis (~120-140 deg): the beam on the aortic valve, the
// plane holding the root axis — LVOT, the cusps, sinuses and the proximal
// ascending aorta running to the screen-right.
function meAvLaxProbe() {
  const tgt = vadd(A, vscale(LM.U_AO, 0.5)), pos = teePathPoint(MEAVLAX_UP);
  const dir = norm(vsub(tgt, pos));
  return teeClear(probeFrom(tgt, dir, vsub(LM.U_AO, vscale(dir, vdot(LM.U_AO, dir))), Math.hypot(...vsub(tgt, pos))));
}
const MEAVLAX_UP = 1.0;                              // cm withdrawn: the beam meets the root side-on
// ME mitral commissural (ME4C + MECOMM_OMNI and a MECOMM_TURN turn of the shaft): the
// plane through both commissures — P3 on the left, A2 in the middle, P1 on the
// right, the papillary muscles on the walls below.
const meCommProbe = meView('MECOMM');
// Transgastric views are manipulations of the one TG SAX probe position (the
// fundal contact on the gastric path), so they lie on the path by construction
// and follow the anatomy: a small search picks each preset's controls.
const cachedProbe = (build) => { let c = null; return () => { if (!c) c = build(); return { ...c, pos: c.pos.slice(), dir: c.dir.slice(), lat: c.lat.slice(), normal: c.normal.slice(), target: c.target.slice() }; }; };
const distToPlane = (p, q) => Math.abs(vdot(vsub(q, p.pos), p.normal));
const inPlaneAng = (p, q) => { const w = vsub(q, p.pos); return Math.atan2(vdot(w, p.lat), vdot(w, p.dir)); };
// TG two-chamber (~90 deg): the plane holding both papillary muscles (end-diastole
// and end-systole) with the LV long axis: the inferior wall in the near field,
// the anterior wall beyond, both papillary muscles with their chordae running to
// the mitral leaflets. Searched over the omniplane angle, a little flex, a
// shaft turn and up to TG2C_ADV cm of advance / withdrawal along the gastric path,
// for the heads (the worst of the two muscles), the mid-bodies (TG2C_W[0]), the
// mitral chordae fanning from the heads to the leaflets (TG2C_W[1] x the share of
// chordal samples > 0.3 cm off the plane), the apex (TG2C_W[2] beyond 0.5 cm), and
// the heads and the mitral valve inside the sector (TG2C_W[3]: within 38 deg of the
// beam, 1 cm to depth - 1 cm; a plane can hold them all behind or beside the
// transducer): the muscles then read as columns on the wall with their chordal
// fans, not as free head islands. A head more than TG2C_HEAD cm off the plane
// costs TG2C_W[4] times more. (The four head positions (two muscles, ED and ES)
// are not coplanar: the heads move ~0.3 cm inward in systole. With the valve
// inside the sector the best common plane holds them within ~0.36 cm; the search
// settles at 0.37 / 0.38 cm with about half the chordae in plane at end-diastole.)
// (From the TG mid station, ~3 cm from the apex,
// the LV spans more than the 90 deg sector: the apex lies off its left edge.)
const TG2C_ADV = [-2, 3], TG2C_W = [0.5, 3.0, 0.1, 1.0, 8.0], TG2C_HEAD = 0.3, TG2C_DEPTH = 13;
let _tg2cG = null;
const tg2cGeo = () => _tg2cG || (_tg2cG = [geometryAt(0, {}), geometryAt(0.42, {})]);   // end-diastole, end-systole (as A0, A_ES)
const tg2cProbe = cachedProbe(() => {
  const tg = tgsaxProbe(), Gs = tg2cGeo();
  const heads = Gs.flatMap((G) => G.A.pap.map((pp) => pp.b));
  const mids = Gs.flatMap((G) => G.A.pap.map((pp) => vscale(vadd(pp.a, pp.b), 0.5)));
  const chords = Gs.flatMap((G) => G.chordae.filter((c) => !c.tv).flatMap((c) => [0, 0.25, 0.5, 0.75, 1].map((t) => vadd(c.a, vscale(vsub(c.b, c.a), t)))));
  let best = null;
  const tryIt = (adv, th, f, rot) => {
    if (adv < TG2C_ADV[0] || adv > TG2C_ADV[1]) return;
    const p = teeManipulate(tg, { omni: th, flex: f, advance: adv, rotate: rot });
    const m = Math.max(...heads.map((h) => distToPlane(p, h))), mm = Math.max(...mids.map((h) => distToPlane(p, h)));
    const off = chords.reduce((n, q) => n + (distToPlane(p, q) > 0.3 ? 1 : 0), 0) / chords.length;
    const outside = (q) => { const a = Math.abs(inPlaneAng(p, q)) * 180 / Math.PI, z = vdot(vsub(q, p.pos), p.dir); return Math.max(0, a - 38) / 10 + Math.max(0, 1 - z) + Math.max(0, z - (TG2C_DEPTH - 1)); };
    const s = m + TG2C_W[4] * Math.max(0, m - TG2C_HEAD) + TG2C_W[0] * mm + TG2C_W[1] * off + TG2C_W[2] * Math.max(0, distToPlane(p, APEX) - 0.5)
      + TG2C_W[3] * [...heads, M].reduce((a, q) => a + outside(q), 0)
      + 0.01 * (Math.abs(adv) + Math.abs(f) / 10 + Math.abs(rot) / 5);
    if (!best || s < best.s) best = { s, p, k: [adv, th, f, rot] };
  };
  // coarse grid, then a fine one about the best
  for (let adv = TG2C_ADV[0]; adv <= TG2C_ADV[1]; adv += 0.5) for (let th = 40; th <= 140; th += 4) for (let f = -30; f <= 30; f += 6) for (const rot of [-30, -20, -10, 0, 10, 20, 30]) tryIt(adv, th, f, rot);
  const [a0, t0, f0, r0] = best.k;
  for (let adv = a0 - 1; adv <= a0 + 1; adv += 0.5) for (let th = t0 - 4; th <= t0 + 4; th += 1) for (let f = f0 - 6; f <= f0 + 6; f += 2) for (let rot = r0 - 5; rot <= r0 + 5; rot += 2.5) tryIt(adv, th, f, rot);
  return best.p;
});
// TG long axis (~120-140 deg): from the TG mid station (the TG SAX contact, the
// probe advanced no more than TGLAX_ADV[1] cm along the greater curvature), with
// the flex and turn of the probe controls (+/-30 deg), the omniplane turned until
// the plane holds the LVOT, the cusps in profile and the proximal ascending aorta
// (the root axis within a few degrees of the plane), the aortic valve in the far
// field to the screen-right. From this station the beam meets the LVOT flow at
// ~40 deg (the model's root / LVOT geometry); the deep TG 5C, advanced to the
// apex, is the window aligned with ejection for the AV CW / LVOT PW pair.
const TGLAX_ADV = [0, 1.5], TGLAX_OMNI = [90, 160], TGLAX_OMNI_REF = 125;   // (the readout kept near the textbook 120-140 deg)
const tgLaxProbe = cachedProbe(() => {
  const tg = tgsaxProbe(), ax = norm(vsub(A, LM.LVOT0));
  let best = null;
  for (let adv = TGLAX_ADV[0]; adv <= TGLAX_ADV[1]; adv += 0.5) for (let f = -30; f <= 30; f += 3) for (let rot = -30; rot <= 30; rot += 5) {
    const b = teeManipulate(tg, { advance: adv, flex: f, rotate: rot });
    for (let th = 0; th < 180; th += 2) {
      const p = teeManipulate(b, { omni: th }), om = teeOmni(p);
      if (om < TGLAX_OMNI[0] || om > TGLAX_OMNI[1]) continue;
      const avOff = distToPlane(p, A), avAng = inPlaneAng(p, A) * 180 / Math.PI;
      if (avOff > 0.3 || avAng < 3 || avAng > 38) continue;
      const rootOut = Math.asin(Math.min(1, Math.abs(vdot(LM.U_AO, p.normal)))) * 180 / Math.PI;
      const beam = Math.acos(Math.min(1, Math.abs(vdot(norm(vsub(A, p.pos)), ax)))) * 180 / Math.PI;
      const s = Math.max(0, rootOut - 3) + 0.15 * beam + avOff + 0.02 * (Math.abs(f) + Math.abs(rot)) + 0.03 * Math.max(0, Math.abs(om - TGLAX_OMNI_REF) - 10);
      if (!best || s < best.s) best = { s, p };
    }
  }
  return best.p;
});
// body-classifier areas (cm2) of a probe's sector at a coarse step (preset searches)
function sectorAreas(p, A, depth, h = 0.3) {
  const c = {};
  for (let d = h; d < depth; d += h) for (let x = -d; x <= d; x += h) {
    if (Math.hypot(x, d) > depth || Math.abs(Math.atan2(x, d)) > TEE_HALF) continue;
    const q = vadd(p.pos, vadd(vscale(p.dir, d), vscale(p.lat, x)));
    const k = bodyClassify(q[0], q[1], q[2], A).code;
    c[k] = (c[k] || 0) + h * h;
  }
  return (k) => c[k] || 0;
}
// TG basal short axis (~0 deg): the TG SAX probe withdrawn toward the junction,
// anteflexed and turned until the plane cuts the mitral leaflets (within 0.5 cm
// of the mitral centre, the fish-mouth) across the LV long axis, the most LV and
// the least atrial / aortic lumen in the sector.
const tgBasalSaxProbe = cachedProbe(() => {
  const tg = tgsaxProbe(), axl = norm(vsub(M, APEX));
  let best = null;
  for (let adv = -3; adv <= 1; adv += 0.5) for (let f = -10; f <= 50; f += 4) for (const rot of [-10, 0, 10]) {
    const p = teeManipulate(tg, { advance: adv, flex: f, rotate: rot });
    const mo = distToPlane(p, M);
    if (mo > 0.5) continue;
    const tilt = Math.acos(Math.min(1, Math.abs(vdot(axl, p.normal)))) * 180 / Math.PI;
    const c = sectorAreas(p, A0, 12, 0.5);
    const s = -c(BODY.LV) + 3 * mo + 0.2 * Math.max(0, tilt - 15) + 0.3 * (c(BODY.AO) + c(BODY.RA) + c(BODY.LA));
    if (!best || s < best.s) best = { s, p };
  }
  return best.p;
});
// TG RV inflow (~100-160 deg): the probe at the TG mid station (the TG SAX
// contact) at the two-chamber angle, turned to the patient's right onto the
// tricuspid valve (within 0.8 cm of the plane), the omniplane and flex adjusted
// for the most RV (end-diastole and end-systole) and RA: the RA, the tricuspid
// leaflets and their chordae, the RV inflow and apex.
const tgRvInProbe = cachedProbe(() => {
  const b = teeManipulate(tgsaxProbe(), { omni: 90 }), om0 = teeOmni(b);
  let best = null;
  for (let rot = 20; rot <= 80; rot += 10) for (let d = 0; d <= 90; d += 15) for (const f of [-25, -10, 5]) {
    const p = teeManipulate(b, { rotate: rot, omni: om0 + d, flex: f });
    const tv = distToPlane(p, T);
    if (tv > 0.8) continue;
    const c0 = sectorAreas(p, A0, 14, 0.5), c1 = sectorAreas(p, A_ES, 14, 0.5);
    const s = -Math.min(c0(BODY.RV), 2.5 * c1(BODY.RV)) - 0.5 * Math.min(c0(BODY.RA), 6) + 0.3 * (c0(BODY.LV) + 0.3 * c0(BODY.LUNG)) + 2 * tv;
    if (!best || s < best.s) best = { s, p };
  }
  return best.p;
});
// ME ascending-aortic short axis (~0-40 deg): the probe withdrawn up the
// oesophagus toward the level of the right pulmonary artery, the beam on the
// ascending aorta (between the STJ and the arch landmark) and the plane turned to
// cut it square (the aorta round) with the RPA in plane behind it, the main PA
// beside it. ME ascending-aortic long axis: the same probe at +90 deg, the
// ascending aorta along its length with the RPA in cross-section behind it.
// (From the arch landmarks, so they follow the arch.)
function ascAxis() { return norm(vsub(LM.AO_ASC_TOP, LM.AO_STJ)); }
const meAscSaxProbe = cachedProbe(() => {
  const U = ascAxis(), rpa = vadd(LM.PA_BIF, vscale(LM.U_RPA, 1.5));
  let best = null;
  for (let s = 2; s <= 9; s += 0.5) {
    const pos = teePathPoint(s);
    for (let h = 0; h <= 3; h += 0.5) {
      const Pa = vadd(LM.AO_STJ, vscale(U, h)), dir = norm(vsub(Pa, pos));
      const base = probeFrom(Pa, dir, teeFrame(dir).lft, Math.hypot(...vsub(Pa, pos)));
      for (let o = 0; o < 180; o += 3) {
        const p = teeManipulate(base, { omni: o });
        const tilt = Math.acos(Math.min(1, Math.abs(vdot(p.normal, U)))) * 180 / Math.PI;
        const sc = tilt / 10 + distToPlane(p, rpa);
        if (!best || sc < best.sc) best = { sc, p };
      }
    }
  }
  return best.p;
});
const meAscLaxProbe = cachedProbe(() => { const p = meAscSaxProbe(); return teeManipulate(p, { omni: teeOmni(p) + 90 }); });
// Upper-oesophageal aortic arch long axis (~0 deg): the probe withdrawn behind the
// arch (UE_ARCH_DROP cm below the level of the arch summit landmark), the beam on
// the summit and the plane near transverse, at the omniplane angle within
// UE_ARCH_OMNI that lays the most aortic lumen in the sector (the arch rises
// steeply from the isthmus to the summit here, so the plane tips a little toward
// it); UE arch short axis: the same probe at +90 deg, the arch in cross-section.
// (From the arch landmarks, so they follow the arch.)
const UE_ARCH_DROP = 1.0, UE_ARCH_OMNI = [0, 20], UE_ARCH_DEPTH = 8;
const ueArchLaxProbe = cachedProbe(() => {
  const s = vdot(vsub(LM.AO_ARCH_MID, ESO), _S) - UE_ARCH_DROP;
  const pos = teePathPoint(Math.min(TEE_PATH_TOP, s)), tgt = LM.AO_ARCH_MID;
  const dir = norm(vsub(tgt, pos));
  const base = probeFrom(tgt, dir, teeFrame(dir).lft, Math.hypot(...vsub(tgt, pos)));
  const lumen = (p, h = 0.2) => {
    let n = 0;
    for (let d = h; d < UE_ARCH_DEPTH; d += h) for (let x = -d; x <= d; x += h) {
      if (Math.hypot(x, d) > UE_ARCH_DEPTH || Math.abs(Math.atan2(x, d)) > TEE_HALF) continue;
      if (lumenDist(...vadd(p.pos, vadd(vscale(p.dir, d), vscale(p.lat, x))), A0, 'AO') < 0) n++;
    }
    return n * h * h;
  };
  let best = null;
  for (let o = UE_ARCH_OMNI[0]; o <= UE_ARCH_OMNI[1]; o += 1) {
    const p = teeManipulate(base, { omni: o }), a = lumen(p);
    if (!best || a > best.a) best = { a, p };
  }
  return best.p;
});
const ueArchSaxProbe = cachedProbe(() => { const p = ueArchLaxProbe(); return teeManipulate(p, { omni: teeOmni(p) + 90 }); });
// Deep transgastric five-chamber (0 deg): the probe advanced from the fundus
// into the gastric body beneath the apex and anteflexed, so the beam runs up the
// LV long axis through the LVOT to the aortic valve — aligned with ejection,
// the TEE window for the aortic-valve gradient (CW). The probe rests DTG_BACK cm
// beyond the apex on the LVOT axis.
const DTG_BACK = 2.2;
function dtg5cProbe() {
  const ax = norm(vsub(A, APEX)), pos = vsub(APEX, vscale(ax, DTG_BACK));
  const dir = norm(vsub(A, pos));
  return probeFrom(A, dir, vsub(vsub(M, T), vscale(dir, vdot(vsub(M, T), dir))), Math.hypot(...vsub(A, pos)));
}
export const TEE_EXTRA = {
  MEAVLAX: { probe: meAvLaxProbe, depth: 10, track: (A) => [0, 0.8 * (A.lv.A[1] - A0.lv.A[1]), 0] },
  MECOMM: { probe: meCommProbe, depth: 15 },
  TG2C: { probe: tg2cProbe, depth: TG2C_DEPTH },
  TGLAX: { probe: tgLaxProbe, depth: 15 },
  DTG5C: { probe: dtg5cProbe, depth: 16 },
  DESCAOLAX: { probe: omniOf(descAoSaxProbe, 90), depth: 8 },
  // review 12: the rest of the ASE/SCA set
  MEASCSAX: { probe: meAscSaxProbe, depth: 10 },
  MEASCLAX: { probe: meAscLaxProbe, depth: 10 },
  // (the plane follows the mitral annulus down in systole along its normal, as the
  // operator follows the base with slight flexion: PSAX-MV does the same)
  TGBSAX: { probe: tgBasalSaxProbe, depth: 12, track: (A) => { const n = tgBasalSaxProbe().normal; return vscale(n, vdot(vsub(A.valves.mitral.c, A0.valves.mitral.c), n)); } },
  TGRVIN: { probe: tgRvInProbe, depth: 14 },
  UEARCHLAX: { probe: ueArchLaxProbe, depth: UE_ARCH_DEPTH },
  UEARCHSAX: { probe: ueArchSaxProbe, depth: UE_ARCH_DEPTH },
};
Object.assign(TEE_VIEWS, TEE_EXTRA);
// Every TEE view carries its half-angle, which the renderer, the Doppler site
// geometry and the audits read through sectorHalfOf.
for (const v of Object.values(TEE_VIEWS)) v.half = TEE_HALF;
export const sectorHalfOf = (view) => (view && view.half) || TTE_HALF;
// the oesophageal / gastric path every TEE probe position lies on: the
// oesophagus from the gastro-oesophageal junction up past the arch, then the
// stomach from the junction to the fundal contact (TG views) and on to the
// deep transgastric position beneath the apex
export function teePathDist(q) {
  let best = 1e9;
  for (let s = TEE_GEJ; s <= TEE_PATH_TOP; s += 0.1) best = Math.min(best, Math.hypot(...vsub(q, teePathPoint(s))));
  const g = gastricPolyline();
  for (let i = 0; i + 1 < g.length; i++) {
    const ab = vsub(g[i + 1], g[i]), t = Math.max(0, Math.min(1, vdot(vsub(q, g[i]), ab) / vdot(ab, ab)));
    best = Math.min(best, Math.hypot(...vsub(q, vadd(g[i], vscale(ab, t)))));
  }
  return best;
}
// The TEE preset nearest a (manipulated) probe: transducer distance (cm) and the
// angles (deg) between the beams and between the imaging planes (normals, either
// sense). `match` when the probe is at that preset's station and plane — the live
// view is then named after it (e.g. ME4C turned to 60 deg reads ME commissural).
export function nearestTeePreset(p) {
  let best = null;
  for (const [name, v] of Object.entries(TEE_VIEWS)) {
    const q = v.probe();
    const dPos = Math.hypot(...vsub(p.pos, q.pos));
    const dDir = Math.acos(Math.max(-1, Math.min(1, vdot(p.dir, q.dir)))) * 180 / Math.PI;
    const dPl = Math.acos(Math.min(1, Math.abs(vdot(p.normal, q.normal)))) * 180 / Math.PI;
    const score = dPos * 10 + dDir + dPl;
    if (!best || score < best.score) best = { name, dPos, dDir, dPl, score };
  }
  best.match = best.dPos < 0.6 && best.dDir < 8 && best.dPl < 8;
  return best;
}
export const ALL_VIEWS = { ...TTE_VIEWS, ...EXTRA_VIEWS, ...TEE_VIEWS };
