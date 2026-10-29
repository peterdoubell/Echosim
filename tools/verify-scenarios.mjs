// Selected pathology subcases must change the physiology, not just their labels.
// Run with: node tools/verify-scenarios.mjs
import assert from 'node:assert/strict';
import { Heartbeat, vcgToLeads } from '../js/electrophysiology.js';
import { hemoSummary, geometryAt, regionalStrain } from '../js/cardiac-model.js';
import { rwmaBlend } from '../js/anatomy.js';
import { DCM_SEVERITIES } from '../js/scenarios.js';

let checks = 0;
function check(name, test) { assert.ok(test, name); checks++; console.log(`PASS ${name}`); }
const normal = new Heartbeat(); normal.advance(0.12);
const baseline = { ...normal.leads() };
// Expected clinical lead groups are deliberately independent of the UI config.
const territories = {
  septal: { leads: ['V1', 'V2'], wall: 'mid anteroseptal', opposite: 'mid inferolateral', angle: 5 * Math.PI / 6 },
  anteroseptal: { leads: ['V1', 'V2', 'V3', 'V4'], wall: 'mid anteroseptal', opposite: 'mid inferior', angle: 3 * Math.PI / 4 },
  anterior: { leads: ['V3', 'V4'], wall: 'mid anterior', opposite: 'mid inferior', angle: Math.PI / 2 },
  lateral: { leads: ['I', 'aVL', 'V5', 'V6'], reciprocal: ['III', 'aVF'], wall: 'mid anterolateral', opposite: 'mid anteroseptal', angle: 0 },
  inferior: { leads: ['II', 'III', 'aVF'], reciprocal: ['I', 'aVL'], wall: 'mid inferior', opposite: 'mid anterior', angle: -Math.PI / 2 },
  inferolateral: { leads: ['II', 'III', 'aVF', 'I', 'V5', 'V6'], reciprocal: ['V1', 'V2'], wall: 'mid inferolateral', opposite: 'mid anteroseptal', angle: -Math.PI / 4 },
};
for (const [region, expected] of Object.entries(territories)) {
  const path = { rwma: region };
  const heart = new Heartbeat(); heart.advance(0.12);
  const phase = heart.mechanical().phase;
  heart.setPathology(path); // must update an already-running beat
  const leads = heart.leads();
  check(`${region}: selected lead group shows ST injury`, expected.leads.every((l) => leads[l] - baseline[l] > 0.1));
  if (expected.reciprocal) check(`${region}: reciprocal ST depression`, expected.reciprocal.every((l) => leads[l] - baseline[l] < -0.1));
  check(`${region}: changing case keeps cardiac time/phase`, heart._t === 0.12 && heart.mechanical().phase === phase);
  check(`${region}: all leads retain electrical identities`, Math.abs(leads.II - leads.I - leads.III) < 1e-12
    && Math.abs(leads.aVL - (leads.I - leads.II / 2)) < 1e-12);
  const strain = regionalStrain(path).segments;
  const wall = strain.find((s) => s.name === expected.wall).strain;
  const opposite = strain.find((s) => s.name === expected.opposite).strain;
  check(`${region}: affected wall has less shortening`, wall > -8 && opposite < -16);
  const A = geometryAt(hemoSummary(path).tMinVol, path).A;
  const a = expected.angle, y = A.lv.apexY + 5;
  const affected = rwmaBlend(2 * Math.cos(a), y, 2 * Math.sin(a), A);
  const remote = rwmaBlend(-2 * Math.cos(a), y, -2 * Math.sin(a), A);
  check(`${region}: rendered myocardial deformation follows selected wall`, affected > 0.65 && remote < 0.05);
  heart.setPathology({});
  check(`${region}: clearing case removes injury immediately`, Math.abs(heart.leads().V2 - baseline.V2) < 1e-12);
}

// Explicit STEMI rhythm must not silently replace an inferior case with an
// anterior ECG. It remains an anterior example when there is no selected RWMA.
{
  const heart = new Heartbeat(); heart.setRhythm('stemi'); heart.setPathology({ rwma: 'inferior' }); heart.advance(0.12);
  check('regional injury takes precedence over generic STEMI example', heart.leads().III - baseline.III > 0.4 && heart.leads().aVL - baseline.aVL < -0.3);
}

let lastEF = Infinity, lastEDV = 0, lastStrain = -Infinity;
const results = {};
for (const severity of ['borderline', 'impaired', 'severe']) {
  const path = { dilated: true, dcmSeverity: severity };
  const h = hemoSummary(path), range = DCM_SEVERITIES[severity];
  check(`DCM ${severity}: measured EF in BSE band`, h.EF >= range.minEF && h.EF <= range.maxEF && Math.abs(h.EF - range.targetEF) < 1);
  const ed = geometryAt(0, path), es = geometryAt(h.tMinVol, path);
  const volume = (g) => g.lv.r.reduce((a, v) => a * v, 1);
  const visualEF = 100 * (1 - volume(es) / volume(ed));
  check(`DCM ${severity}: geometry shortening matches measured EF`, Math.abs(visualEF - h.EF) < 1);
  const strains = regionalStrain(path).segments.map((s) => s.strain);
  check(`DCM ${severity}: global rather than regional dysfunction`, Math.max(...strains) - Math.min(...strains) < 1e-12);
  check(`DCM ${severity}: progressive dilatation and loss of function`, h.EF < lastEF && h.EDV > lastEDV && strains[0] > lastStrain);
  results[severity] = h.EF;
  lastEF = h.EF; lastEDV = h.EDV; lastStrain = strains[0];
  console.log(`     EF ${h.EF.toFixed(1)}%, EDV ${h.EDV.toFixed(0)} mL, GLS ${strains[0].toFixed(1)}%`);
}
check('default DCM retains the original severe teaching case', hemoSummary({ dilated: true }).EF === results.severe);
for (const severity of ['impaired', 'borderline', 'severe']) {
  check(`DCM ${severity}: switching back returns the correct cached model`, hemoSummary({ dilated: true, dcmSeverity: severity }).EF === results[severity]);
}

// At 150 bpm the six-second display contains ~15 beats, so retaining just six
// beats would erase most of a paused strip when a region is selected.
{
  const heart = new Heartbeat(), control = new Heartbeat();
  heart.setHR(150); control.setHR(150);
  heart.advance(8); control.advance(8);
  heart.leads();
  const time = heart.time, phase = heart.mechanical().phase;
  const cachedVcg = [...heart.vcg];
  const earlyBeat = heart._beats.find((b) => b.t >= time - 6);
  check('tachycardic history covers the full six-second display', earlyBeat.t < time - 5.5 && heart._beats.length > 12);
  const before = vcgToLeads(heart.vcgAt(earlyBeat.t + 0.12));
  const out = [0, 0, 0];
  check('retrospective sample supports a reusable output vector', heart.vcgAt(earlyBeat.t + 0.12, out) === out);
  heart.vcgAt(time + 0.03);
  check('past/future sampling leaves current caches and phase unchanged', cachedVcg.every((v, i) => v === heart.vcg[i])
    && heart.time === time && heart.mechanical().phase === phase);
  heart.setPathology({ rwma: 'inferior' });
  const after = vcgToLeads(heart.vcgAt(earlyBeat.t + 0.12));
  check('paused history changes immediately to the selected region', after.III - before.III > 0.4 && after.aVL - before.aVL < -0.3);
  check('paused history repaint does not advance cardiac time', heart.time === time && heart.mechanical().phase === phase);
  heart.setPathology({}); heart.advance(1); control.advance(1);
  const actual = heart.leads(), expected = control.leads();
  check('retrospective reads do not affect later ECG or rhythm', heart.time === control.time
    && heart.mechanical().phase === control.mechanical().phase && Object.keys(actual).every((l) => Math.abs(actual[l] - expected[l]) < 1e-12));
}
console.log(`\n${checks} scenario checks passed.`);
