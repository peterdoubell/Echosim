// thorax3d.js — the thoracic context of the 3D panel: lungs, diaphragm, liver,
// oesophagus and stomach, the pericardial sac and the chest wall, as fields for
// the surface-nets mesher (marching.js), plus the structure labels and the
// camera presets. Kept free of Three.js so the node audits can mesh, count and
// project them; heart3d.js turns the meshes into translucent context.
//
// The lungs and the liver come from the same body classifier the echo uses
// (anatomy.js bodyClassify), so the 3D context is the tissue the imaging plane
// cuts through; the stomach and the diaphragm come from their shared SDFs.

import { bodyClassify, BODY, BODY_AX, LM, diaphragmBelow, stomachDist, epiDist } from './anatomy.js';
import { TEE_VIEWS } from './views.js';

const { L, P, S, A: AN } = BODY_AX;
const MIDLINE = -2.0;                                   // p . L of the mid-sagittal plane (anatomy.js)
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const unit = (a) => { const l = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0] / l, a[1] / l, a[2] / l]; };
const mad = (a, d, k) => [a[0] + d[0] * k, a[1] + d[1] * k, a[2] + d[2] * k];
// body-frame coordinates of a point: l (+ patient left of the midline), a (+ anterior), s (+ superior)
const bodyOf = (x, y, z) => [x * L[0] + y * L[1] + z * L[2] - MIDLINE, -(x * P[0] + y * P[1] + z * P[2]), x * S[0] + y * S[1] + z * S[2]];
const world = (l, a, s) => [0, 1, 2].map((i) => (l + MIDLINE) * L[i] + a * AN[i] + s * S[i]);

// Thoracic cavity: a rounded (superelliptic) cross-section ~26 cm wide and
// ~16 cm deep, the anterior chest wall ~1.6 cm under the parasternal window,
// tapering to the lung apices above the clavicles; the diaphragm closes it below.
const THX = { ac: -1.5, rl: 13, ra: 8, p: 2.4, sTop: 16, sTaper: 4, sLow: -14 };
// approximate signed distance (cm, negative inside) of the cavity grown by `grow` cm
function cavity(l, a, s, grow = 0) {
  const k = s > THX.sTaper ? Math.sqrt(Math.max(0.02, 1 - ((s - THX.sTaper) / (THX.sTop - THX.sTaper)) ** 2)) : 1;
  const rl = THX.rl * k + grow, ra = THX.ra * (0.55 + 0.45 * k) + grow;
  const e = Math.pow(Math.pow(Math.abs(l) / rl, THX.p) + Math.pow(Math.abs(a - THX.ac) / ra, THX.p), 1 / THX.p);
  return Math.max((e - 1) * Math.min(rl, ra), THX.sLow - grow - s, s - THX.sTop - grow);
}

// The world box that holds the whole context (the body-frame box l -15..15,
// a -12..10, s -19..17 turned into world axes) and its voxel size.
export const CONTEXT_CELL = 0.7;
export const CONTEXT_BOX = (() => {
  const mn = [1e9, 1e9, 1e9], mx = [-1e9, -1e9, -1e9];
  for (const l of [-15.5, 15.5]) for (const a of [-12, 10.5]) for (const s of [-19, 17]) {
    const w = world(l, a, s); for (let i = 0; i < 3; i++) { mn[i] = Math.min(mn[i], w[i]); mx[i] = Math.max(mx[i], w[i]); }
  }
  return { min: mn, max: mx };
})();

// Oesophagus: the tube the TEE probe travels in, through the mid-oesophageal
// and descending-aorta probe stations, down through the hiatus to the cardia
// (and the transgastric probe station in the stomach); radius ~0.8 cm.
function oesophagusPath() {
  const me = TEE_VIEWS.ME4C.probe().pos, da = TEE_VIEWS.DESCAO.probe().pos, tg = TEE_VIEWS.TGSAX.probe().pos;
  const top = mad(me, S, 13);
  const hiatus = mad(mad(da, S, -4.2), L, 0.6);
  return [top, me, da, hiatus, tg];
}
function capsule(p, a, b, r) {
  const ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]], ap = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
  const t = Math.max(0, Math.min(1, dot(ap, ab) / dot(ab, ab)));
  return Math.hypot(ap[0] - ab[0] * t, ap[1] - ab[1] * t, ap[2] - ab[2] * t) - r;
}

// Context organs, drawn back to front (renderOrder) under the heart. Colours are
// the panel's: grey-blue aerated lung, brown liver, muscle-red diaphragm,
// pink-tan gut, pale fibrous pericardium, and a faint chest wall.
export const CONTEXT = [
  { key: 'chestWall', color: 0x5d7088, opacity: 0.07, back: true },
  { key: 'lungs', color: 0x9fb0c8, opacity: 0.13 },
  { key: 'diaphragm', color: 0xa8584a, opacity: 0.28 },
  { key: 'liver', color: 0x7a3a22, opacity: 0.38 },
  { key: 'stomach', color: 0xc99a7c, opacity: 0.32 },
  { key: 'oesophagus', color: 0xd09080, opacity: 0.5 },
  { key: 'pericardium', color: 0xeadfc6, opacity: 0.09 },
];

// Build the context meshes from the end-diastolic anatomy A with the mesher
// `buildMesh(field, box, cell)`. The lungs and liver are classifier labels, so
// they are sampled once on the grid as indicators, box-blurred (3x3x3) and meshed
// at the half level; the rest are fields. Returns { key: mesh | null }.
export function buildContext(A, buildMesh, heartBox) {
  const box = CONTEXT_BOX, c = CONTEXT_CELL;
  const nx = Math.ceil((box.max[0] - box.min[0]) / c) + 1, ny = Math.ceil((box.max[1] - box.min[1]) / c) + 1, nz = Math.ceil((box.max[2] - box.min[2]) / c) + 1;
  const lung = new Float32Array(nx * ny * nz), liver = new Float32Array(nx * ny * nz);
  const below = new Float32Array(nx * ny * nz), cav = new Float32Array(nx * ny * nz);
  let q = 0;
  for (let k = 0; k < nz; k++) for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++, q++) {
    const x = box.min[0] + i * c, y = box.min[1] + j * c, z = box.min[2] + k * c;
    const [l, a, s] = bodyOf(x, y, z);
    const cv = cavity(l, a, s); cav[q] = cv;
    const db = diaphragmBelow(x, y, z); below[q] = db;
    if (cv > 0.5 || s < -19) continue;                    // outside the thorax and upper abdomen
    const code = bodyClassify(x, y, z, A).code;
    if (code === BODY.LUNG && db < -0.3) lung[q] = 1;     // (gastric gas is labelled LUNG below the dome)
    else if (code === BODY.LIVER) liver[q] = 1;
  }
  const blur = (f) => {
    const t = new Float32Array(f.length), o = new Float32Array(f.length);
    const pass = (src, dst, st, n, dim) => {
      for (let p = 0; p < src.length; p++) {
        const idx = dim === 0 ? p % nx : dim === 1 ? ((p / nx) | 0) % ny : (p / (nx * ny)) | 0;
        let s = src[p] * 2, w = 2;
        if (idx > 0) { s += src[p - st]; w++; }
        if (idx < n - 1) { s += src[p + st]; w++; }
        dst[p] = s / w;
      }
    };
    pass(f, t, 1, nx, 0); pass(t, o, nx, ny, 1); pass(o, t, nx * ny, nz, 2);
    return t;
  };
  const lungB = blur(lung), liverB = blur(liver);
  const at = (x, y, z) => {
    const i = Math.round((x - box.min[0]) / c), j = Math.round((y - box.min[1]) / c), k = Math.round((z - box.min[2]) / c);
    return i < 0 || j < 0 || k < 0 || i >= nx || j >= ny || k >= nz ? -1 : (k * ny + j) * nx + i;
  };
  const grid = (arr, outside) => (x, y, z) => { const p = at(x, y, z); return p < 0 ? outside : arr[p]; };
  const lungF = grid(lungB, 0), liverF = grid(liverB, 0), belowF = grid(below, 0), cavF = grid(cav, 1);
  // distance to the dome sheet: the height below it over its gradient (the height
  // field climbs steeply where the sheet drapes round the sac)
  const sheetF = (x, y, z) => {
    const p = at(x, y, z); if (p < 0) return 1e3;
    const i = p % nx, j = ((p / nx) | 0) % ny, k = (p / (nx * ny)) | 0;
    const g = (d, n, st) => (d > 0 && d < n - 1 ? (below[p + st] - below[p - st]) / (2 * c) : 0);
    const gm = Math.max(0.35, Math.hypot(g(i, nx, 1), g(j, ny, nx), g(k, nz, nx * ny)));
    return Math.abs(below[p] + 0.2) / gm;
  };
  const oes = oesophagusPath();
  const fields = {
    // the chest wall: the outer skin ~1.8 cm outside the cavity (drawn from the inside)
    chestWall: (x, y, z) => { const [l, a, s] = bodyOf(x, y, z); return Math.max(cavity(l, a, s, 1.8), -15 - s, s - 12); },
    lungs: (x, y, z) => 0.5 - lungF(x, y, z),
    // the muscular sheet about the dome surface, inside the cavity's outline (drawn
    // ~1 cm thick: a thinner sheet falls between the 7 mm voxels and perforates)
    diaphragm: (x, y, z) => Math.max(sheetF(x, y, z) - 0.5, cavF(x, y, z)),
    liver: (x, y, z) => 0.5 - liverF(x, y, z),
    stomach: (x, y, z) => Math.max(stomachDist(x, y, z) - 0.45, -belowF(x, y, z)),
    oesophagus: (x, y, z) => {
      let d = 1e9; for (let n = 0; n + 1 < oes.length; n++) d = Math.min(d, capsule([x, y, z], oes[n], oes[n + 1], 0.8));
      return d;
    },
  };
  const out = {};
  for (const [k, f] of Object.entries(fields)) out[k] = buildMesh(f, box, c);
  // the fibrous pericardium: a sac ~3 mm outside the epicardium of the chambers
  // and the roots of the great vessels, meshed in the (finer) heart box
  out.pericardium = heartBox ? buildMesh((x, y, z) => epiDist(x, y, z, A) - 0.3, heartBox, 0.45) : null;
  return out;
}

// Structure labels of the 3D panel: [text, world anchor] from the end-diastolic
// anatomy. The left atrial appendage at its tip and the four pulmonary veins at
// their hilar ends (both hard to pick out of the LA otherwise), the cavae and
// great arteries, and the context organs.
export function structureLabels(A) {
  const out = [];
  const aa = A.la.aa; if (aa && aa.length) out.push(['LAA', aa[aa.length - 1].b]);
  for (const v of (A.pv && A.pv.veins) || []) out.push([v.name, v.b]);
  out.push(['SVC', A.svc.b]);
  const ivcTop = dot(A.ivc.a, S) > dot(A.ivc.b, S) ? A.ivc.a : A.ivc.b;      // its RA end
  out.push(['IVC', mad(ivcTop, S, -3.5)]);
  out.push(['Ao', mad(A.valves.aortic.c, LM.U_AO, 4.5)]);
  out.push(['PA', mad(A.valves.pulmonic.c, LM.U_PA, 2.5)]);
  const ctx = [
    ['L lung', world(10, 1, 3)], ['R lung', world(-10, 1, 3)], ['Liver', world(-7, 2, -10)],
    ['Diaphragm', world(-11, -2, -7)], ['Stomach', world(6, 0, -9)], ['Oesophagus', mad(TEE_VIEWS.ME4C.probe().pos, S, 11)],
  ];
  return { heart: out, context: ctx };
}

// 3D camera presets per view, in body axes (l, a, s): the camera sits on that
// side of the heart with patient-superior up, in the anterior hemisphere and
// never far left of the sternal plane, so the SVC always projects to the
// patient's right of the ascending aorta (screen left, as on a CT or an anterior
// anatomical view) — the left-lateral cameras used before showed it on the
// aortic arch's concave side, as if left of the aorta.
export const CAMERA_PRESETS = {
  PLAX: [-0.3, 1, 0.35], PSAX: [0.1, 1, 0.05], PSAX_AV: [0, 1, 0.3], A4C: [0.12, 0.85, -0.5], A2C: [0.15, 0.9, -0.3],
  SUBCOSTAL: [-0.3, 0.6, -0.75], default: [0.05, 1, 0.25],
};
export const CAMERA_UP = S;
export function cameraDir(name) {
  const b = CAMERA_PRESETS[name] || CAMERA_PRESETS.default;
  return unit([0, 1, 2].map((i) => b[0] * L[i] + b[1] * AN[i] + b[2] * S[i]));
}
