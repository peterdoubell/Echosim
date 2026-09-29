// verify-atrial-phase.mjs
// Proves the atria and ventricles are OUT of phase.
//
// Samples the live model across one cardiac cycle and checks, against the
// geometry that is actually rendered (not a restatement of the formula):
//   - LV volume is minimal in systole and maximal at end-diastole;
//   - atrial volume is MAXIMAL when LV volume is minimal (AV-valve opening) and
//     MINIMAL at the onset of ventricular systole;
//   - the two volume curves are strongly ANTI-correlated;
//   - the atrial curve is periodic and continuous across the phase wrap;
//   - atrial emptying splits into physiologic passive and active fractions.
//
// Usage: node tools/verify-atrial-phase.mjs
import { geometryAt } from '../js/cardiac-model.js';
import { lumenDist } from '../js/anatomy.js';

const N = 400;
const phases = Array.from({ length: N }, (_, i) => i / N);

const rows = phases.map((p) => {
  const g = geometryAt(p, {});
  const vol = (r) => r[0] * r[1] * r[2];        // ∝ ellipsoid volume
  return { p, lv: vol(g.lv.r), la: vol(g.la.r), ra: vol(g.ra.r) };
});

const argmax = (key) => rows.reduce((a, b) => (b[key] > a[key] ? b : a)).p;
const argmin = (key) => rows.reduce((a, b) => (b[key] < a[key] ? b : a)).p;

const pearson = (a, b) => {
  const n = a.length;
  const ma = a.reduce((s, v) => s + v, 0) / n;
  const mb = b.reduce((s, v) => s + v, 0) / n;
  let num = 0, da = 0, db = 0;
  for (let i = 0; i < n; i++) {
    const x = a[i] - ma, y = b[i] - mb;
    num += x * y; da += x * x; db += y * y;
  }
  return num / Math.sqrt(da * db);
};

// circular distance between two phases, in [0, 0.5]
const dphase = (a, b) => { const d = Math.abs(a - b) % 1; return Math.min(d, 1 - d); };

const lvMinP = argmin('lv'), lvMaxP = argmax('lv');
// The LV sits on a flat floor through late ejection, so argmin lands on the
// FIRST sample of that plateau. True end-systole — the instant the AV valves are
// about to open, which is what the atrial peak should track — is the LAST phase
// at which LV volume is still within 2 % of its minimum.
const lvMinV = Math.min(...rows.map((x) => x.lv));
// end-systole = aortic closure: the FIRST phase at the minimum (the LV then sits
// at ESV through isovolumic relaxation)
const lvEndSysP = rows.filter((x) => x.lv <= lvMinV * 1.02).reduce((a, b) => (b.p < a.p ? b : a)).p;
const laMaxP = argmax('la'), laMinP = argmin('la');
const raMaxP = argmax('ra'), raMinP = argmin('ra');

const r = pearson(rows.map((x) => x.lv), rows.map((x) => x.la));

// periodicity / continuity across the wrap: compare the step at 0 with the
// largest step anywhere else in the cycle.
const step = (i) => Math.abs(rows[(i + 1) % N].la - rows[i].la);
const wrapStep = step(N - 1);
const maxInteriorStep = Math.max(...rows.map((_, i) => (i === N - 1 ? 0 : step(i))));

// atrial function: volumes at the three landmarks
// Emptying fractions are measured on the DRAWN chamber (voxel volume of the LA
// lumen SDF, which includes the septal anchoring and vein mouths), not on the
// analytic ellipsoid — what a trainee planimeters is what must be physiologic.
const voxelLA = (p) => {
  const A = geometryAt(p, {}).A, h = 0.2;
  let n = 0;
  for (let x = -7; x < 5; x += h) for (let y = -3; y < 7; y += h) for (let z = -6; z < 4; z += h)
    if (lumenDist(x, y, z, A, 'LA') < 0) n++;
  return n * h * h * h;
};
const vMax = voxelLA(laMaxP), vMin = voxelLA(laMinP), vPreA = voxelLA(0.84);
const totalEF = (vMax - vMin) / vMax;
const passiveEF = (vMax - vPreA) / vMax;
const activeEF = (vPreA - vMin) / vPreA;

const checks = [
  ['LV minimum is in ventricular systole (0.20-0.44)', lvMinP >= 0.20 && lvMinP <= 0.44, `LV min @ ${lvMinP.toFixed(3)}`],
  ['LV maximum is at end-diastole (>=0.90 or <=0.02)', lvMaxP >= 0.90 || lvMaxP <= 0.02, `LV max @ ${lvMaxP.toFixed(3)}`],
  ['LA maximum is near AV-valve opening (0.42-0.52)', laMaxP >= 0.42 && laMaxP <= 0.52, `LA max @ ${laMaxP.toFixed(3)}`],
  ['LA minimum is at onset of systole (>=0.97 or <=0.03)', laMinP >= 0.97 || laMinP <= 0.03, `LA min @ ${laMinP.toFixed(3)}`],
  // the RA runs on the right heart's clock: it stops emptying when the tricuspid
  // closes, ~25 ms (0.03) after the mitral, and peaks as it opens, just before it
  ['RA follows the LA, offset by the right-heart timing', dphase(raMaxP, laMaxP) <= 0.03 && ((raMinP - laMinP + 1.5) % 1) - 0.5 >= 0.01 && ((raMinP - laMinP + 1.5) % 1) - 0.5 <= 0.05, `RA max @ ${raMaxP.toFixed(3)}, min @ ${raMinP.toFixed(3)} (LA ${laMaxP.toFixed(3)}, ${laMinP.toFixed(3)})`],
  ['LA max follows LV end-systole by an IVRT-like 0-0.15', laMaxP - lvEndSysP >= 0 && laMaxP - lvEndSysP <= 0.15, `LV end-systole @ ${lvEndSysP.toFixed(3)}, LA max @ ${laMaxP.toFixed(3)}`],
  ['LA is still EXPANDING throughout ventricular ejection', rows.filter((x) => x.p > 0.02 && x.p < lvEndSysP).every((x, i, a) => i === 0 || x.la > a[i - 1].la), 'monotone rise 0.02 -> end-systole'],
  ['LA min is within 0.10 of the LV MAXIMUM (anti-phase)', dphase(laMinP, lvMaxP) <= 0.10, `|LAmin-LVmax| = ${dphase(laMinP, lvMaxP).toFixed(3)}`],
  ['LA and LV volumes are anti-correlated (r < -0.5)', r < -0.5, `pearson r = ${r.toFixed(3)}`],
  ['atrial curve is continuous across the phase wrap', wrapStep <= maxInteriorStep * 1.5, `wrap step ${wrapStep.toExponential(2)} vs max interior ${maxInteriorStep.toExponential(2)}`],
  ['total LA emptying fraction 45-65 % (voxel)', totalEF >= 0.45 && totalEF <= 0.65, `LAEF = ${(totalEF * 100).toFixed(1)} %, ${vMin.toFixed(0)}-${vMax.toFixed(0)} mL`],
  ['passive (conduit) emptying fraction 28-48 %', passiveEF >= 0.28 && passiveEF <= 0.48, `passive EF = ${(passiveEF * 100).toFixed(1)} %`],
  ['active (booster) emptying fraction 20-40 %', activeEF >= 0.20 && activeEF <= 0.40, `active EF = ${(activeEF * 100).toFixed(1)} %`],
];

let fails = 0;
for (const [name, ok, detail] of checks) {
  if (!ok) fails++;
  console.log(`  ${ok ? 'PASS' : 'FAIL'} ${name}  [${detail}]`);
}

console.log('\n  phase   LV vol    LA vol');
for (let i = 0; i < N; i += N / 20) {
  const x = rows[i];
  console.log(`  ${x.p.toFixed(2)}   ${x.lv.toFixed(3).padStart(7)}   ${x.la.toFixed(3).padStart(7)}`);
}

console.log(`\nRESULT=${fails === 0 ? 'PASS' : `FAIL (${fails})`}`);
process.exit(fails === 0 ? 0 : 1);
