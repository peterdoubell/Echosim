// anatomy.js — an anatomically-structured signed-distance model of the heart.
// Replaces the old "stack of ellipsoids" with a topologically correct body:
// a continuous myocardium, a crescentic RV wrapping a shared interventricular
// septum, atria with appendages, an aortic root with sinuses of Valsalva + arch,
// the main pulmonary artery crossing anteriorly, papillary muscles and the RV
// moderator band. Both the 2D echo cross-section and the 3D marching-cubes
// surface sample this same field, so the two views are always consistent.
//
// Heart space (cm): +y base (top), -y apex (bottom), +x patient-left (LV side),
// +z anterior (toward the transducer).

import {
  smin, ssub, sdSphere, sdEllipsoid, sdCapsule, sdRoundCone,
} from './sdf.js';

// ---- BSE/ASE normal reference ranges -------------------------------------
// Adult echocardiographic normal ranges, exported so the Measurements UI can
// draw the "normal band" next to each value. Units are centimetres except EF
// (%) and areas (cm^2). Source: BSE/ASE chamber-quantification guidance.
export const REF = {
  lviddNormal:  [4.2, 5.9],  // LV internal diameter, end-diastole (men)
  lvidsNormal:  [2.5, 4.0],  // LV internal diameter, end-systole
  lvWallNormal: [0.6, 1.2],  // IVSd / PWd (septal / posterior wall thickness)
  efNormal:     [55, 70],    // ejection fraction (%)
  laNormal:     [3.0, 4.0],  // LA antero-posterior diameter
  raAreaNormal: [10, 18],    // RA area (<= 18 cm^2 normal)
  aoRootNormal: [2.9, 3.7],  // aortic sinus-of-Valsalva diameter
  paNormal:     [1.5, 2.5],  // main pulmonary-artery diameter
  rvd1Normal:   [2.5, 4.1],  // RV basal diameter (RVD1)
};

// ---- anatomical configuration --------------------------------------------
// End-diastolic ("resting") constants, calibrated toward the BSE/ASE normal
// ranges above. All lengths in centimetres, in heart space (+y base, -y apex,
// +x patient-left / LV side, +z anterior). This is the SINGLE SOURCE OF TRUTH
// for chamber sizes: anatomyParams() scales it per cardiac phase for the 2D/3D
// views, and cardiac-model.js derives its coarse BASE flow/measurement proxies
// from these same numbers (imported) so the panel matches the visual anatomy.
//
// Packing note: the LV minor semi-axis is 2.3 cm -> LVIDd 4.6 cm (mid-normal).
// It is not pushed to the population mean (~4.8) because the fixed 3D voxel box
// (heart3d BBOX) and the fixed 2D scan-plane targets (main.js) frame this
// compact schematic layout; a larger LV would clip the lateral wall / crowd the
// atria. The long axis is likewise kept schematic (base-apex ~6 cm) so the LV
// base does not overrun the atria — EF and LVIDd (the panel values) do not
// depend on the long axis, only on the minor axis and the systolic ratios.
const CFG = {
  lv:  { c: [1.2, -1.15, 0.15], r: [2.3, 3.05, 2.3], wall: 0.9 },
  // RV = two crescent lobes hugging the LV. Centres are pulled toward the shared
  // septum (x) and well anterior (+z, toward the transducer) and the lobes are
  // widened + splayed so PSAX reads as a BROAD crescent wrapping the septum (a
  // sweeping arc, not a compact lateral lobe). The long-axis semi-axis r[1] is
  // kept short so in A4C the RV stops well short of the LV apex (a normal RV is
  // shorter than the LV). The large rvBlend fuses the two lobes into one arc and
  // the ssub carve (rvCarve) keeps that arc hugging the interventricular septum.
  rvA: { c: [-0.55, -0.5, 1.4],  r: [2.05, 2.5, 1.8] },
  rvB: { c: [-1.6, 0.05, 1.1],   r: [1.8, 2.05, 1.55] },
  rvWall:  0.38,   // RV free-wall thickness (normal <= 0.5 cm)
  rvBlend: 1.4,    // smin radius fusing the two RV lobes into one broad crescent
  rvCarve: 0.3,    // ssub radius carving the LV epicardium out of the RV crescent
  la:  { c: [1.15, 2.7, -1.0], r: 1.75,   // radius 1.75 -> LA diameter 3.5 cm
         // Left atrial appendage. NOT a straight tube: a narrow oval ostium
         // (~1.5-2.0 cm long-axis clinically) opens into a neck, then a multi-lobed
         // body that hooks back on itself — the "chicken-wing" morphology, the
         // commonest of the four described variants (~48%). It sits anterolaterally
         // in the left AV groove, overlying the circumflex. The narrow neck is why
         // it is THE site of thrombus in atrial fibrillation, and why it looks so
         // different from the wide-mouthed RA appendage. Three tapering segments:
         // ostium -> neck, neck -> body lobe, then the hooked apical lobe.
         aa: [
           { a: [1.72, 2.80, -0.70], b: [2.25, 2.68, -0.15], r1: 0.50, r2: 0.38 },
           { a: [2.25, 2.68, -0.15], b: [2.62, 2.34,  0.34], r1: 0.38, r2: 0.30 },
           { a: [2.62, 2.34,  0.34], b: [2.36, 2.86,  0.70], r1: 0.30, r2: 0.19 },
         ] },
  ra:  { c: [-2.0, 2.55, 0.2], r: 1.7,
         // Right atrial appendage: BROAD-BASED and triangular, with a wide mouth
         // continuous with the atrium — the key teaching discriminator from the
         // narrow-necked LAA. One wide cone tapering to a blunt tip.
         aa: [{ a: [-2.50, 2.74, 0.48], b: [-3.28, 2.08, 1.02], r1: 0.85, r2: 0.34 }] },
  // Crista terminalis waypoints, given as UNIT directions from the RA centre and
  // scaled by the live RA radius so the ridge keeps hugging the wall as the
  // atrium fills and empties. The arc runs from the SVC orifice (superior),
  // around the lateral wall, to the IVC (inferior) — the C-shaped muscular ridge
  // dividing the smooth-walled sinus venarum behind from the trabeculated
  // pectinate appendage in front. A prominent crista is a classic echo
  // pseudo-mass mistaken for RA thrombus.
  cristaU: [
    [-0.18,  0.96, 0.07],
    [-0.79,  0.55, 0.10],
    [-0.99,  0.00, 0.07],
    [-0.76, -0.53, 0.03],
    [-0.20, -0.94, 0.00],
  ],
  cristaFrac: 0.86,   // fraction of the RA radius: sits just inside the wall
  cristaR: 0.19,      // ridge thickness ~3.8 mm (crista is typically 3-6 mm)
  // Eustachian valve (valve of the IVC): the crescentic flap at the IVC-RA
  // junction, directed toward the fossa ovalis — in fetal life it streams IVC
  // blood across the foramen ovale. Endpoints are RA-relative like the crista.
  eustachU: { a: [-0.30, -0.90, 0.05], b: [0.42, -0.72, -0.20], fa: 0.88, fb: 0.80, r: 0.10 },
  // Systemic venous returns entering the RA. Short stubs only — enough to give
  // the crista its two anatomical anchors and the Eustachian valve its orifice,
  // and to make the subcostal IVC assessment (a real bedside skill) possible.
  svc: { a: [-2.20, 3.85, 0.18], b: [-2.35, 5.30, 0.12], r: 0.62 },  // ~1.2 cm calibre
  ivc: { a: [-2.15, 1.05, 0.12], b: [-2.30, -0.70, 0.02], r: 0.82 }, // ~1.6 cm calibre
  // Interatrial septum. Modelled explicitly because without it the LA and RA
  // blood pools abut directly — i.e. every heart would have an open ASD. The
  // septum is a disc normal to the LA->RA axis whose thickness ramps from a thin
  // central FOSSA OVALIS membrane out to the thicker muscular limbus. That
  // thin fossa is genuine teaching anatomy: it is where PFO/secundum ASD occur
  // and where transseptal puncture is aimed, and its echo dropout is the classic
  // false-positive for an ASD.
  ias: { r: 1.30, tLimbus: 0.17, tFossa: 0.055, rFossa: 0.60 },
  // Pulmonary veins draining into the POSTERIOR left atrium — left and right,
  // superior and inferior. Ostial calibre here ~0.9-1.0 cm (clinically 1.0-2.0 cm).
  // `u` is a unit direction from the LA centre, so the ostia stay seated on the
  // atrial wall as it fills and empties; each vein then runs outward posteriorly.
  // The right-sided veins enter close to the interatrial septum, which is why the
  // right superior vein is the one seen entering the LA on a standard apical
  // four-chamber view.
  pv: [
    { name: 'LSPV', u: [ 0.45,  0.55, -0.70], r: 0.53, len: 1.5 },
    { name: 'LIPV', u: [ 0.58, -0.42, -0.70], r: 0.48, len: 1.4 },
    { name: 'RSPV', u: [-0.50,  0.50, -0.70], r: 0.55, len: 1.5 },
    { name: 'RIPV', u: [-0.55, -0.45, -0.70], r: 0.50, len: 1.4 },
  ],
  // The muscular ridge between the left superior pulmonary vein and the LAA
  // ostium — the "warfarin" (coumadin) ridge. It is a normal structure with a
  // bulbous tip that is regularly mistaken for a thrombus or mass on TEE, so it
  // is worth showing rather than smoothing away.
  pvRidgeR: 0.13,
  aoValve:  [0.42, 1.5, 0.14], // aortic-valve centre; the three sinuses ring it
  // Sinus radius enlarged from the old 0.66, but capped here: in this compressed
  // schematic the aortic root sits at the crux, so a full BSE-size root (~1.6 cm
  // radius) would balloon into the interventricular septum / RV inflow. 0.75 ->
  // an effective root diameter ~2.5 cm, the largest that stays clear of the crux.
  aoSinusR: 0.75,
  aoWall:   0.2,               // aortic-root wall thickness
  paWall:   0.16,              // pulmonary-artery wall thickness
  dilatedScale: 1.4,           // LV dilatation factor -> dilated LVIDd ~6.4 cm (severe)
  laDilation:   1.16,          // LA enlargement factor in dilated cardiomyopathy
};

// Exported so cardiac-model.js can build its coarse BASE ellipsoid proxies from
// the same numbers (single source of truth; no silent geometry drift).
export { CFG };

// Interatrial-septum frame, derived once from the atrial centres: the plate sits
// midway between the atria with its normal along the LA->RA axis, so it is the
// true dividing surface rather than an axis-aligned approximation.
const IAS_C = [
  (CFG.la.c[0] + CFG.ra.c[0]) / 2,
  (CFG.la.c[1] + CFG.ra.c[1]) / 2,
  (CFG.la.c[2] + CFG.ra.c[2]) / 2,
];
const IAS_N = (() => {
  const d = [CFG.ra.c[0] - CFG.la.c[0], CFG.ra.c[1] - CFG.la.c[1], CFG.ra.c[2] - CFG.la.c[2]];
  const l = Math.hypot(d[0], d[1], d[2]) || 1;
  return [d[0] / l, d[1] / l, d[2] / l];
})();

// Place the crista-terminalis waypoints against the live RA wall and return the
// consecutive segments that make up the ridge.
function buildCrista(c, r) {
  const f = CFG.cristaFrac * r, pts = CFG.cristaU;
  const segs = [];
  for (let i = 0; i < pts.length - 1; i++) {
    const p = pts[i], q = pts[i + 1];
    segs.push({
      a: [c[0] + p[0] * f, c[1] + p[1] * f, c[2] + p[2] * f],
      b: [c[0] + q[0] * f, c[1] + q[1] * f, c[2] + q[2] * f],
      r: CFG.cristaR,
    });
  }
  return segs;
}

// Seat the pulmonary-vein ostia on the live LA wall and run each vein outward
// posteriorly. Returns segments plus the LSPV/LAA ridge that sits between the
// left superior vein and the appendage mouth.
function buildPV(c, r) {
  const veins = CFG.pv.map((v) => {
    const ul = Math.hypot(v.u[0], v.u[1], v.u[2]) || 1;
    return {
      name: v.name, r: v.r,
      a: [c[0] + v.u[0] * r * 0.80, c[1] + v.u[1] * r * 0.80, c[2] + v.u[2] * r * 0.80],
      b: [c[0] + v.u[0] * (r + v.len), c[1] + v.u[1] * (r + v.len), c[2] + v.u[2] * (r + v.len)],
      // unit INFLOW axis (distal vein -> atrium), precomputed: the Doppler field
      // reads this per sample, so normalising here keeps it out of the hot loop
      ax: [-v.u[0] / ul, -v.u[1] / ul, -v.u[2] / ul],
    };
  });
  // ridge: spans the LA wall between the LSPV ostium and the LAA ostium
  const lspv = veins[0].a, laa = CFG.la.aa[0].a;
  const ridge = {
    a: [(lspv[0] + c[0] * 0.12) / 1.12, (lspv[1] + c[1] * 0.12) / 1.12, (lspv[2] + c[2] * 0.12) / 1.12],
    b: [laa[0], laa[1], laa[2]],
    r: CFG.pvRidgeR,
  };
  return { veins, ridge };
}

// The Eustachian valve, likewise anchored to the live RA wall at the IVC orifice.
function buildEustachian(c, r) {
  const e = CFG.eustachU, fa = e.fa * r, fb = e.fb * r;
  return {
    a: [c[0] + e.a[0] * fa, c[1] + e.a[1] * fa, c[2] + e.a[2] * fa],
    b: [c[0] + e.b[0] * fb, c[1] + e.b[1] * fb, c[2] + e.b[2] * fb],
    r: e.r,
  };
}

// Build the per-phase parameter bundle. k = contraction 0..1, kick = atrial
// kick 0..1, path = pathology flags. `mech` (optional) carries the volume-exact
// LV scaling from the lumped-parameter circulation (sShort/sLong/lvWall) so the
// SDF cavity tracks the modelled PV loop; when omitted the legacy kinematic
// scaling is used (keeps anatomyParams callable standalone / backward-compatible).
export function anatomyParams(k, kick, path = {}, mech = null) {
  let contract = 1.0, lvScale = 1.0, lvWallMul = 1.0, laScale = 1.0;
  if (path.dilated) { contract = 0.26; lvScale = CFG.dilatedScale; lvWallMul = 0.82; laScale = CFG.laDilation; }
  if (path.lvh || path.aorticStenosis) lvWallMul = 1.7;

  // RV pressure/volume overload (pulmonary hypertension): the RV dilates and the
  // interventricular septum flattens/bows toward the LV, so the LV reads D-shaped in
  // PSAX (the "D-sign"). rvpoScale enlarges the RV; septFlat clips the LV septal side.
  const rvpoScale = path.rvpo ? 1.45 : 1;
  const septFlat = path.rvpo ? 0.42 : 0;   // pronounced septal flattening → clear D-sign

  // RV / atrial cavity scale: shrinks in systole (short axis more than long)
  const sx = 1 - k * 0.30 * contract, sy = 1 - k * 0.14 * contract, sz = 1 - k * 0.30 * contract;
  // LV cavity scale: volume-exact from the circulation when coupled, else legacy
  const lsx = mech ? mech.sShort : sx, lsy = mech ? mech.sLong : sy, lsz = mech ? mech.sShort : sz;
  const lvWall = mech ? mech.lvWall : CFG.lv.wall * lvWallMul * (1 + k * 0.55 * contract);
  const rvWall = CFG.rvWall * (1 + k * 0.45);

  // Apex-anchored longitudinal contraction: real ventricles contract toward a
  // near-fixed apex, so the mitral annulus (LV base) descends by the FULL long-axis
  // shortening (MAPSE ~12-15 mm) rather than half of it. The apex stays put and the
  // centre shifts apically as the long axis shortens; this is volume-preserving
  // (only the centre moves, the radii are unchanged) so EF/LVIDd/LVIDs are intact.
  const r1d = CFG.lv.r[1] * lvScale;         // end-diastolic LV long semi-axis
  const apexY = CFG.lv.c[1] - r1d;           // fixed apex (−y end)
  const dyLV = r1d * (lsy - 1);              // centre apical shift (≤0 in systole)
  const lvCy = CFG.lv.c[1] + dyLV;           // apex-anchored LV centre y
  // material map: a point at end-diastolic height y0 moves with the apex-anchored
  // wall to apexY + (y0−apexY)·lsy. Used for the papillary muscles + valve plane.
  const axialMap = (y0) => apexY + (y0 - apexY) * lsy;

  // atria: reservoir fills through systole, empties in diastole, kicks at end
  const aFill = (0.86 + 0.32 * (1 - k) - 0.34 * kick);

  // aortic root: the ascending column is angled slightly rightward (-x) and
  // anterior (+z) so its septal (anterior) wall reads continuous with the
  // interventricular septum — aortomitral / septal continuity.
  const av = CFG.aoValve;
  const ascB = [0.05, 3.7, 0.45];   // top of the ascending column (arch start)

  return {
    k, kick, contract,
    // LV cavity: `r` are the LIVE (per-phase) semi-axes; `rDia`/`wallDia` are the
    // end-diastolic (unshortened, un-thickened) reference used by the regional
    // wall-motion (RWMA) blend, and `rwma` describes the hypokinetic/akinetic
    // segment (null when normal). See rwmaBlend() in the LV distance fields.
    lv: {
      c: [CFG.lv.c[0], lvCy, CFG.lv.c[2]],   // apex-anchored centre (base descends → MAPSE)
      r: [CFG.lv.r[0] * lvScale * lsx, CFG.lv.r[1] * lvScale * lsy, CFG.lv.r[2] * lvScale * lsz],
      wall: lvWall,
      rDia: [CFG.lv.r[0] * lvScale, CFG.lv.r[1] * lvScale, CFG.lv.r[2] * lvScale],
      wallDia: CFG.lv.wall * lvWallMul,
      rwma: buildRwma(path),
      septFlat,   // >0 flattens the LV septal (−x) side → D-shaped LV (RV overload)
    },
    // RV built from two blended blobs (carved by the LV epicardium at classify).
    // rvpoScale dilates the RV in pressure/volume overload.
    rvA: { c: CFG.rvA.c, r: [CFG.rvA.r[0] * sx * rvpoScale, CFG.rvA.r[1] * sy * rvpoScale, CFG.rvA.r[2] * sz * rvpoScale] },
    rvB: { c: CFG.rvB.c, r: [CFG.rvB.r[0] * sx * rvpoScale, CFG.rvB.r[1] * sy * rvpoScale, CFG.rvB.r[2] * sz * rvpoScale] },
    rvWall, rvBlend: CFG.rvBlend, rvCarve: CFG.rvCarve,
    rvot: { a: [-1.25, 0.55, 0.92], b: [-0.95, 1.95, 1.0], r1: 0.9 * (1 - k * 0.2), r2: 0.8 },
    la: { c: CFG.la.c, r: CFG.la.r * laScale * aFill, aa: CFG.la.aa },
    // pulmonary veins + the LSPV/LAA ridge, seated on the live atrial wall
    pv: buildPV(CFG.la.c, CFG.la.r * laScale * aFill),
    ra: { c: CFG.ra.c, r: CFG.ra.r * aFill, aa: CFG.ra.aa },
    // right-atrial internal anatomy, rebuilt per phase so the crista keeps
    // hugging the (breathing) atrial wall rather than floating in the lumen
    crista: buildCrista(CFG.ra.c, CFG.ra.r * aFill),
    eustachian: buildEustachian(CFG.ra.c, CFG.ra.r * aFill),
    svc: CFG.svc, ivc: CFG.ivc,
    // interatrial septum: a disc normal to the LA->RA axis, thin over the fossa
    // ovalis. A secundum ASD is modelled by driving the fossa thickness NEGATIVE,
    // which perforates the membrane exactly where such defects actually occur.
    ias: {
      c: IAS_C, n: IAS_N, r: CFG.ias.r,
      tLimbus: CFG.ias.tLimbus, rFossa: CFG.ias.rFossa,
      tFossa: path.asd ? -0.03 : CFG.ias.tFossa,
    },
    // aortic root: three sinuses of Valsalva around the aortic valve, then a
    // tapered ascending column and a short arch.
    ao: {
      valve: av,
      // sinusR enlarged toward the BSE sinus-of-Valsalva diameter; the three
      // overlapping sinus spheres give an effective root diameter ~2.8 cm.
      sinus: [[av[0] + 0.6, av[1] + 0.07, av[2]], [av[0] - 0.35, av[1] + 0.07, av[2] + 0.52], [av[0] - 0.35, av[1] + 0.07, av[2] - 0.52]],
      sinusR: CFG.aoSinusR,
      // ascending column ~1.9 cm across; slightly narrower than the sinuses
      // (sino-tubular junction) then tapering up into the arch.
      asc: { a: [av[0] - 0.02, 1.95, av[2] + 0.04], b: ascB, r1: 0.95, r2: 0.88 },
      arch: [[ascB, [-0.3, 4.2, 0.0], 0.84], [[-0.3, 4.2, 0.0], [-1.35, 4.0, -0.55], 0.8]],
      wall: CFG.aoWall,
    },
    // main pulmonary artery: crosses anterior-leftward over the aorta, bifurcates.
    // Main-trunk radius ~1.1 -> ~2.2 cm diameter (BSE normal <= 2.5 cm).
    pa: {
      main: { a: [-0.95, 1.95, 1.0], b: [-0.15, 3.35, 0.55], r1: 1.1 * (1 - k * 0.15), r2: 1.0 },
      branch: [[[-0.15, 3.35, 0.55], [-1.2, 3.5, 0.2], 0.62], [[-0.15, 3.35, 0.55], [0.7, 3.4, 0.7], 0.57]],
      wall: CFG.paWall,
    },
    // papillary muscles: thick bellies (r~0.45) seated squarely in the PSAX
    // mid-cavity short-axis plane (y ~ -1.2), tapering up toward the mitral
    // valve. Their bases are spread WIDE in x and one is offset anteriorly (+z)
    // vs the other (posterior, -z) so in the short-axis cut they sit at roughly
    // 4 and 8 o'clock about the LV, not bunched together inferiorly.
    // papillary muscles follow the apex-anchored wall (their apical base barely
    // moves; the belly descends with the annulus) via the material map.
    pap: [
      { a: [0.05, axialMap(-2.6), 0.2],  b: [0.5, axialMap(0.2), -0.05], r: 0.45 },
      { a: [2.4, axialMap(-2.6), -0.65], b: [1.8, axialMap(0.2), -0.5],  r: 0.45 },
    ],
    // moderator band: crosses the RV cavity near the apex. Endpoints are nudged
    // toward the A4C plane's centre-z and thickened slightly so the band reads
    // clearly crossing the RV apex within the (anterior-tipped) A4C scan plane.
    mod: { a: [-1.05, -1.8, 1.2], b: [-2.1, -1.7, 1.1], r: 0.34 },
    // RV anterior papillary muscle — anatomically CONTINUOUS with the moderator
    // band (its base is the band's free-wall insertion) and rising toward the
    // tricuspid, whose anterior-leaflet chordae it anchors. Completes the
    // septum → moderator band → papillary → tricuspid apparatus.
    rvPap: { a: [-2.1, -1.7, 1.1], b: [-2.02, -0.28, 0.78], r: 0.3 },
    // apex-anchor parameters so the valve planes (cardiac-model) descend with the
    // mitral/tricuspid annulus by the same material map: y → apexY + (y−apexY)·lsy.
    axial: { apexY, lsy },
  };
}

// ---- regional wall-motion abnormality (RWMA) -----------------------------
// A per-segment contractility defect (ischemia / infarct). Each named segment
// maps to an angular sector about the LV long axis (short-axis x–z plane), or to
// the apical cap. `severity` 0..1: 0 = normal, 1 = akinetic (the wall keeps its
// end-diastolic position and does not thicken). path.rwma may be `true`, a region
// name string, or { region, severity }.
const RWMA_REGIONS = {
  septal:       { ang: Math.PI,          half: 1.15 },
  anteroseptal: { ang: Math.PI * 0.75,   half: 1.20 },
  anterior:     { ang: Math.PI / 2,      half: 1.15 },
  lateral:      { ang: 0,                half: 1.15 },
  inferior:     { ang: -Math.PI / 2,     half: 1.15 },
  posterior:    { ang: -Math.PI / 2,     half: 1.15 },
  apical:       { apical: true },
};
function buildRwma(path) {
  if (!path.rwma) return null;
  const spec = typeof path.rwma === 'object' ? path.rwma : {};
  const name = spec.region || (typeof path.rwma === 'string' ? path.rwma : 'septal');
  const base = RWMA_REGIONS[name] || RWMA_REGIONS.septal;
  const sev = spec.severity != null ? Math.max(0, Math.min(1, spec.severity)) : 0.85;
  return Object.assign({ sev, region: name }, base);
}
// Blend weight g in [0,1]: how "diastolic" (non-contracting) this point's wall is.
function rwmaBlend(x, y, z, A) {
  const rw = A.lv.rwma;
  if (!rw) return 0;
  const c = A.lv.c;
  if (rw.apical) {
    // apex is toward -y: the defect runs from mid-cavity to the apical cap
    return rw.sev * Math.max(0, Math.min(1, (c[1] - y) / 2.2));
  }
  let d = Math.atan2(z - c[2], x - c[0]) - rw.ang;
  while (d > Math.PI) d -= 2 * Math.PI;
  while (d < -Math.PI) d += 2 * Math.PI;
  const t = d / rw.half;
  return rw.sev * Math.exp(-t * t);
}

// ---- component distance fields -------------------------------------------
// The LV is bullet-shaped, not a symmetric ellipsoid: the base (+y) is broad and
// the apex (−y) tapers. sdLVbody narrows the short-axis (x,z) radii toward the
// apex so the classified cavity and the 3D surface read as a real ventricle. The
// numeric EF/LVIDd come from the ellipsoid proxy in cardiac-model (unchanged), so
// this is a morphology-only refinement. APEX_TAPER = short-axis fraction at the tip.
const APEX_TAPER = 0.28;
function lvTaper(y, cy, ry) {
  const u = (y - cy) / ry;                 // −1 ≈ apex, +1 ≈ base
  const s = u < 0 ? Math.min(1, -u / 0.85) : 0;
  return 1 - APEX_TAPER * s * s;           // short-axis → ~0.72× near the apex
}
function sdLVbody(x, y, z, cx, cy, cz, rx, ry, rz) {
  const t = lvTaper(y, cy, ry);
  return sdEllipsoid(x, y, z, cx, cy, cz, rx * t, ry, rz * t);
}
// Flatten the LV septal (−x) side by intersecting with a half-space at
// xp = cx − rx(1 − sf); sf>0 pushes the septum toward the LV → D-shaped cavity.
function septClip(d, x, cx, rx, sf) {
  if (!sf) return d;
  const xp = cx - rx * (1 - sf);
  return Math.max(d, xp - x);
}
function dLVlumen(x, y, z, A) {
  const c = A.lv.c, r = A.lv.r;
  const g = A.lv.rwma ? rwmaBlend(x, y, z, A) : 0;
  if (g > 0) {
    // in the akinetic segment the endocardium stays at its end-diastolic radius
    // (it does not move inward in systole) -> a hypokinetic/akinetic wall.
    const rd = A.lv.rDia;
    return sdLVbody(x, y, z, c[0], c[1], c[2],
      r[0] + (rd[0] - r[0]) * g, r[1] + (rd[1] - r[1]) * g, r[2] + (rd[2] - r[2]) * g);
  }
  return septClip(sdLVbody(x, y, z, c[0], c[1], c[2], r[0], r[1], r[2]), x, c[0], r[0], A.lv.septFlat);
}
function lvEpi(x, y, z, A) {
  const c = A.lv.c, r = A.lv.r, w = A.lv.wall;
  const g = A.lv.rwma ? rwmaBlend(x, y, z, A) : 0;
  if (g > 0) {
    // segment wall also stops thickening: endocardium + wall both revert to the
    // end-diastolic state, so the epicardium barely moves (no systolic excursion).
    const rd = A.lv.rDia, we = w + (A.lv.wallDia - w) * g;
    return sdLVbody(x, y, z, c[0], c[1], c[2],
      r[0] + (rd[0] - r[0]) * g + we, r[1] + (rd[1] - r[1]) * g + we, r[2] + (rd[2] - r[2]) * g + we);
  }
  return septClip(sdLVbody(x, y, z, c[0], c[1], c[2], r[0] + w, r[1] + w, r[2] + w), x, c[0], r[0] + w, A.lv.septFlat);
}
// gentler apex taper for the RV so it converges to a triangular apex in A4C
// (the LV forms the true apex; the RV apex sits just short of it) while the
// short-axis crescent in PSAX — a cut near the base — stays broad.
function sdRVbody(x, y, z, c, r) {
  const u = (y - c[1]) / r[1];               // −1 ≈ apex, +1 ≈ base
  const s = u < 0 ? Math.min(1, -u / 0.9) : 0;
  const t = 1 - 0.20 * s * s;                // short-axis → ~0.80× near the apex
  return sdEllipsoid(x, y, z, c[0], c[1], c[2], r[0] * t, r[1], r[2] * t);
}
function dRVlumen(x, y, z, A) {
  const a = sdRVbody(x, y, z, A.rvA.c, A.rvA.r);
  const b = sdRVbody(x, y, z, A.rvB.c, A.rvB.r);
  let d = smin(a, b, A.rvBlend);
  const rvot = sdRoundCone(x, y, z, A.rvot.a[0], A.rvot.a[1], A.rvot.a[2], A.rvot.b[0], A.rvot.b[1], A.rvot.b[2], A.rvot.r1, A.rvot.r2);
  d = smin(d, rvot, 0.5);
  // carve the RV around the LV so it hugs the interventricular septum; a small
  // carve radius keeps the RV a continuous thin crescent wrapping the septum
  return ssub(d, lvEpi(x, y, z, A), A.rvCarve);
}
// An appendage lumen: a chain of tapering lobes. The lobes are joined tightly
// (small smin) so the individual lobes stay legible instead of melting into one
// sausage — the lobulation is the anatomy worth teaching.
function dAppendage(x, y, z, segs) {
  let d = 1e9;
  for (let i = 0; i < segs.length; i++) {
    const s = segs[i];
    d = smin(d, sdRoundCone(x, y, z, s.a[0], s.a[1], s.a[2], s.b[0], s.b[1], s.b[2], s.r1, s.r2), 0.22);
  }
  return d;
}
// Interatrial septum: a disc normal to the LA->RA axis whose half-thickness ramps
// from the thin central fossa-ovalis membrane out to the muscular limbus. When
// the fossa thickness is negative (ASD) the membrane is perforated, leaving a
// hole at exactly the site secundum defects occur.
function dIAS(x, y, z, A) {
  const s = A.ias;
  const dx = x - s.c[0], dy = y - s.c[1], dz = z - s.c[2];
  const u = dx * s.n[0] + dy * s.n[1] + dz * s.n[2];        // along the septal normal
  const px = dx - u * s.n[0], py = dy - u * s.n[1], pz = dz - u * s.n[2];
  const rad = Math.sqrt(px * px + py * py + pz * pz);        // in-plane radius
  const span = s.r - s.rFossa;
  const f = rad <= s.rFossa ? 0 : Math.min(1, (rad - s.rFossa) / (span || 1));
  const t = s.tFossa + (s.tLimbus - s.tFossa) * f * f * (3 - 2 * f);
  return Math.max(Math.abs(u) - t, rad - s.r);               // slab ∩ disc
}
// Pulmonary-vein lumens (continuous with the LA). The generous blend radius
// gives each ostium the funnel/flare a real pulmonary vein has where it meets
// the atrium, rather than a tube butted onto a sphere.
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
  const s = sdSphere(x, y, z, A.la.c[0], A.la.c[1], A.la.c[2], A.la.r);
  // Cheap far-field rejection. The appendage, the veins and the septum all lie
  // within ~1.7 cm of the atrial sphere's surface, so a sample further out than
  // that cannot be affected by any of them and needs only the sphere. This runs
  // for every pixel of every frame, so the early-out matters: it keeps the LA's
  // sub-anatomy off the hot path for the (many) samples nowhere near it.
  if (s > 2.2) return s;
  // small blend radius => the appendage joins the LA through a NARROW NECK
  let d = smin(s, dAppendage(x, y, z, A.la.aa), 0.16);
  d = smin(d, dPV(x, y, z, A), 0.38);           // veins flare into the atrium
  return ssub(d, dIAS(x, y, z, A), 0.06);       // the septum is muscle, not blood
}
function dRAlumen(x, y, z, A) {
  const s = sdSphere(x, y, z, A.ra.c[0], A.ra.c[1], A.ra.c[2], A.ra.r);
  // far-field rejection as for the LA, but a wider band: the caval stubs run
  // further from the atrial sphere than the left-sided sub-anatomy does.
  if (s > 3.0) return s;
  // large blend radius => a BROAD-BASED mouth continuous with the atrium
  let d = smin(s, dAppendage(x, y, z, A.ra.aa), 0.7);
  // the systemic veins are continuous with the atrium, not separate pools
  const v = A.svc, w = A.ivc;
  d = smin(d, sdCapsule(x, y, z, v.a[0], v.a[1], v.a[2], v.b[0], v.b[1], v.b[2], v.r), 0.45);
  d = smin(d, sdCapsule(x, y, z, w.a[0], w.a[1], w.a[2], w.b[0], w.b[1], w.b[2], w.r), 0.45);
  return ssub(d, dIAS(x, y, z, A), 0.06);
}
// Crista terminalis + Eustachian valve — muscular/fibrous ridges standing proud
// of the RA wall into the blood pool (same treatment as the papillary muscles).
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
function dAOlumen(x, y, z, A) {
  const ao = A.ao;
  let d = 1e9;
  for (const s of ao.sinus) d = smin(d, sdSphere(x, y, z, s[0], s[1], s[2], ao.sinusR), 0.35);
  d = smin(d, sdRoundCone(x, y, z, ao.asc.a[0], ao.asc.a[1], ao.asc.a[2], ao.asc.b[0], ao.asc.b[1], ao.asc.b[2], ao.asc.r1, ao.asc.r2), 0.4);
  for (const seg of ao.arch) d = smin(d, sdCapsule(x, y, z, seg[0][0], seg[0][1], seg[0][2], seg[1][0], seg[1][1], seg[1][2], seg[2]), 0.4);
  return d;
}
function dPAlumen(x, y, z, A) {
  const pa = A.pa;
  let d = sdRoundCone(x, y, z, pa.main.a[0], pa.main.a[1], pa.main.a[2], pa.main.b[0], pa.main.b[1], pa.main.b[2], pa.main.r1, pa.main.r2);
  for (const seg of pa.branch) d = smin(d, sdCapsule(x, y, z, seg[0][0], seg[0][1], seg[0][2], seg[1][0], seg[1][1], seg[1][2], seg[2]), 0.4);
  return d;
}
function dPap(x, y, z, A) {
  let d = 1e9;
  for (const p of A.pap) d = smin(d, sdCapsule(x, y, z, p.a[0], p.a[1], p.a[2], p.b[0], p.b[1], p.b[2], p.r), 0.3);
  return d;
}
function dMod(x, y, z, A) {
  const m = A.mod;
  let d = sdCapsule(x, y, z, m.a[0], m.a[1], m.a[2], m.b[0], m.b[1], m.b[2], m.r);
  const p = A.rvPap;                       // RV anterior papillary, fused to the band
  if (p) d = smin(d, sdCapsule(x, y, z, p.a[0], p.a[1], p.a[2], p.b[0], p.b[1], p.b[2], p.r), 0.35);
  return d;
}

// Myocardial + atrial shell only (no great vessels). Used as the 3D muscle
// surface so the aorta/PA can render as distinct vessel meshes.
export function myoDist(x, y, z, A) {
  let d = lvEpi(x, y, z, A);
  // RV free wall = RV lumen inflated, still carved around LV so the septum is shared
  const rvEpi = ssub(smin(
    sdEllipsoid(x, y, z, A.rvA.c[0], A.rvA.c[1], A.rvA.c[2], A.rvA.r[0] + A.rvWall, A.rvA.r[1] + A.rvWall, A.rvA.r[2] + A.rvWall),
    sdEllipsoid(x, y, z, A.rvB.c[0], A.rvB.c[1], A.rvB.c[2], A.rvB.r[0] + A.rvWall, A.rvB.r[1] + A.rvWall, A.rvB.r[2] + A.rvWall), A.rvBlend),
    lvEpi(x, y, z, A) - A.lv.wall * 0.4, 0.4);
  d = smin(d, rvEpi, 0.5);
  // atrial walls (thin) — including the appendages, veins and septum, so those
  // are enclosed by wall rather than floating as bare lumen outside the
  // epicardial boundary (which would let pericardial fluid pool over them).
  // Each atrium's sub-anatomy sits behind a cheap bounding test on its sphere:
  // this is the per-sample hot path, and most samples are nowhere near it.
  const sLA = sdSphere(x, y, z, A.la.c[0], A.la.c[1], A.la.c[2], A.la.r + 0.24);
  d = smin(d, sLA, 0.5);
  if (sLA < 2.2) {
    d = smin(d, dAppendage(x, y, z, A.la.aa) - 0.16, 0.4);
    d = smin(d, dPV(x, y, z, A) - 0.14, 0.4);       // thin-walled pulmonary veins
    d = smin(d, dIAS(x, y, z, A), 0.25);            // interatrial septum
  }
  const sRA = sdSphere(x, y, z, A.ra.c[0], A.ra.c[1], A.ra.c[2], A.ra.r + 0.24);
  d = smin(d, sRA, 0.5);
  if (sRA < 3.0) {
    d = smin(d, dAppendage(x, y, z, A.ra.aa) - 0.16, 0.4);
    // thin walls on the systemic venous stubs entering the RA
    const v = A.svc, w = A.ivc;
    d = smin(d, sdCapsule(x, y, z, v.a[0], v.a[1], v.a[2], v.b[0], v.b[1], v.b[2], v.r + 0.14), 0.4);
    d = smin(d, sdCapsule(x, y, z, w.a[0], w.a[1], w.a[2], w.b[0], w.b[1], w.b[2], w.r + 0.14), 0.4);
  }
  // AV / interatrial septum: a thin muscle bridge at the crux where the four
  // chambers meet (mitral & tricuspid annuli, inter-atrial & inter-ventricular
  // septa). Gives the chambers a shared muscular crux instead of a blood gap.
  const crux = sdCapsule(x, y, z, -0.8, 0.6, 0.0, -0.2, 1.6, -0.1, 0.35);
  d = smin(d, crux, 0.35);
  return d;
}

// Whole-heart epicardial surface including the great-vessel walls (used for
// pericardial-effusion classification, which pools outside this boundary).
export function epiDist(x, y, z, A) {
  let d = myoDist(x, y, z, A);
  d = smin(d, dAOlumen(x, y, z, A) - A.ao.wall, 0.4);
  d = smin(d, dPAlumen(x, y, z, A) - A.pa.wall, 0.4);
  return d;
}

// Distance to a named blood pool (for 3D chamber surfaces).
export function lumenDist(x, y, z, A, which) {
  switch (which) {
    case 'LV': return dLVlumen(x, y, z, A);
    case 'RV': return dRVlumen(x, y, z, A);
    case 'LA': return dLAlumen(x, y, z, A);
    case 'RA': return dRAlumen(x, y, z, A);
    case 'AO': return dAOlumen(x, y, z, A);
    case 'PA': return dPAlumen(x, y, z, A);
  }
  return 1e9;
}

// Tissue codes must match cardiac-model TISSUE (imported there). We return a
// small code + echogenicity; cardiac-model maps codes to its TISSUE enum.
export const BODY = { OUT: 0, MYO: 1, LV: 2, RV: 3, LA: 4, RA: 5, AO: 6, PAP: 10 };

// Classify a point as blood / myocardium / vessel (valves & effusion handled by
// the caller). Returns { code, echo }.
export function bodyClassify(x, y, z, A) {
  const dlv = dLVlumen(x, y, z, A);
  const drv = dRVlumen(x, y, z, A);

  // papillary muscles & moderator band read as muscle inside the blood pool
  if (dlv < 0.2 && dPap(x, y, z, A) < 0) return { code: BODY.PAP, echo: 0.7 };
  if (drv < 0.2 && dMod(x, y, z, A) < 0) return { code: BODY.PAP, echo: 0.66 };

  const dla = dLAlumen(x, y, z, A);
  const dra = dRAlumen(x, y, z, A);

  // right-atrial ridges standing proud of the wall: the crista terminalis reads
  // as muscle (a classic pseudo-mass), the Eustachian valve as brighter fibrous
  // tissue. Gated on being in/near the RA so the tests cost nothing elsewhere.
  if (dra < 0.25) {
    if (dCrista(x, y, z, A) < 0) return { code: BODY.PAP, echo: 0.68 };
    if (dEustachian(x, y, z, A) < 0) return { code: BODY.PAP, echo: 0.82 };
  }
  // the LSPV/LAA ("warfarin") ridge stands into the left atrium — a normal
  // structure routinely mistaken for a mass, so it reads as bright tissue
  if (dla < 0.25) {
    const g = A.pv.ridge;
    if (sdCapsule(x, y, z, g.a[0], g.a[1], g.a[2], g.b[0], g.b[1], g.b[2], g.r) < 0) {
      return { code: BODY.PAP, echo: 0.78 };
    }
  }
  const dao = dAOlumen(x, y, z, A);
  const dpa = dPAlumen(x, y, z, A);

  // nearest blood pool that we are inside
  let best = 0.0, code = BODY.OUT;
  const pick = (d, c) => { if (d < best) { best = d; code = c; } };
  pick(dlv, BODY.LV); pick(drv, BODY.RV); pick(dla, BODY.LA); pick(dra, BODY.RA);
  pick(dao, BODY.AO); pick(dpa, BODY.AO);
  if (code !== BODY.OUT) return { code, echo: 0.03 };

  // myocardium: inside the epicardial body but outside every lumen
  if (epiDist(x, y, z, A) < 0) {
    return { code: BODY.MYO, echo: 0.55 };
  }
  return { code: BODY.OUT, echo: 0 };
}
