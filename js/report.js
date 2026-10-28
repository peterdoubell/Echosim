// Structured echocardiography report.
//
// Builds a BSE-style transthoracic report from the simulated heart: the
// measurements a sonographer would take (made on the model with the same
// calipers as tools/verify-anatomy.mjs, so the report and the anatomy audit
// agree), indexed to the patient's body surface area or height, graded against
// the BSE 2020 reference intervals in js/references.js, and written up as
// findings per structure plus a conclusion. Per-pathology templates give the
// learner a scaffold for writing the report themselves.
//
// Everything here is pure: measureModel() samples the model (a few hundred ms,
// cached per pathology), buildReport() grades and words it, and reportText()
// renders plain text for copying. The DOM side lives in js/main.js.

import { geometryAt, classify, TISSUE, hemoSummary, hemodynamics, eaRatio, longitudinalStrain, gradeOf, shuntLabel } from './cardiac-model.js';
import { lumenDist, LM, rapEstimate } from './anatomy.js';
import { TTE_VIEWS } from './views.js';
import {
  REF_SOURCE, bsa, gradeGraded, gradeNormal, gradeLimit, gradeLvef, gradeLavi, lvGeometry,
  gradeAs, gradeMs, gradeRegurg, phProbability, effusionSize, GLS_NORMAL, EA_NORMAL, IVC_NORMAL_MAX,
} from './references.js';

// ---- small vector helpers -----------------------------------------------------
const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const mul = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const unit = (a) => { const l = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0] / l, a[1] / l, a[2] / l]; };
const inPlane = (v, n) => unit(sub(v, mul(n, dot(v, n))));
const Y = [0, 1, 0];

// Caliper: [entry, exit] distances of the first run along unit d whose samples
// satisfy `is` (step 0.01 cm). `gap` lets a caliper cross a short interruption.
function run(G, p0, d, is, maxLen = 16, start = 0, gap = 0, path = {}) {
  let inRun = false, a = 0, lastIn = 0;
  for (let t = start; t <= maxLen; t += 0.01) {
    const p = add(p0, mul(d, t));
    const c = classify(p[0], p[1], p[2], G, path);
    const k = is(c.tissue, p, c.echo);
    if (k && !inRun) { inRun = true; a = t; }
    if (k) lastIn = t;
    if (!k && inRun && t - lastIn > gap) return [a, lastIn + 0.01];
  }
  return inRun ? [a, lastIn] : null;
}
const isT = (...ks) => (t) => ks.includes(t);
const len = (r) => (r ? r[1] - r[0] : null);
// inside a caval tube beyond its RA orifice (excluded from the RA area trace)
function inCapsule(q, v) {
  const ab = sub(v.b, v.a), t = dot(sub(q, v.a), ab) / dot(ab, ab);
  if (t <= 0 || t >= 1) return false;
  return Math.hypot(...sub(q, add(v.a, mul(ab, t)))) < v.r;
}

// ---- model measurements ----------------------------------------------------------
// Returns raw values in clinical units (mm, mL, cm2, m/s, mmHg, %), unindexed.
const _cache = new Map();
export function measureModel(path = {}) {
  const key = JSON.stringify(path);
  if (_cache.has(key)) return _cache.get(key);
  const sum = hemoSummary(path);
  const ES = sum.tMinVol, ED = 0;
  const GED = geometryAt(ED, path), GES = geometryAt(ES, path);
  // the RV's own end-systole (maximum RV contraction) for TAPSE and FAC
  let rvEs = ES, kMax = -1;
  for (let ph = 0.38; ph <= 0.56; ph += 0.01) { const k = geometryAt(ph, path).kRV; if (k > kMax + 1e-6) { kMax = k; rvEs = ph; } }
  const GRES = geometryAt(rvEs, path);
  const m = {};

  // LV linear dimensions: PLAX, at the leaflet tips, perpendicular to the long axis
  const plax = TTE_VIEWS.PLAX.probe();
  const perpAx = inPlane(plax.dir, Y);
  const lvLine = (G) => {
    const p0 = sub([0, G.A.lv.M[1] - 1.3, 0], mul(perpAx, 7));
    const lvRun = run(G, p0, perpAx, isT(TISSUE.LV, TISSUE.VALVE), 16, 0, 0, path);
    if (!lvRun) return null;
    let s0 = lvRun[0];
    while (s0 > 0) { const p = add(p0, mul(perpAx, s0 - 0.01)); if (classify(p[0], p[1], p[2], G, path).tissue !== TISSUE.MYO) break; s0 -= 0.01; }
    const pw = run(G, p0, perpAx, isT(TISSUE.MYO), 16, lvRun[1], 0, path);
    return { lvid: len(lvRun), ivs: lvRun[0] - s0, pw: pw && Math.abs(pw[0] - lvRun[1]) < 0.05 ? len(pw) : null };
  };
  const ed = lvLine(GED), es = lvLine(GES);
  if (ed) { m.lvidd = ed.lvid * 10; m.ivsd = ed.ivs * 10; m.pwd = ed.pw != null ? ed.pw * 10 : null; }
  if (es) m.lvids = es.lvid * 10;
  // LV mass (ASE linear cube formula, g) and relative wall thickness
  if (ed && ed.pw != null) {
    const D = ed.lvid, S = ed.ivs, P = ed.pw;
    m.lvmass = 0.8 * (1.04 * ((D + S + P) ** 3 - D ** 3)) + 0.6;
    m.rwt = (2 * P) / D;
  }

  // LV volumes and EF: biplane-equivalent volumes from the circulation the
  // cavity geometry is driven by (the model's ground truth)
  m.lvedv = sum.EDV; m.lvesv = sum.ESV; m.lvef = sum.EF;
  m.gls = longitudinalStrain(path);
  m.ea = eaRatio(path);
  m.mapse = (GED.A.valves.mitral.c[1] - GES.A.valves.mitral.c[1]) * 10;

  // aortic root (PLAX, perpendicular to the root axis, inner edge, end-diastole)
  const U = LM.U_AO, perpAo = inPlane(plax.dir, U);
  const aoRun = (G, h) => {
    const c = add(G.A.valves.aortic.c, mul(U, h));
    return len(run(G, sub(c, mul(perpAo, 3)), perpAo, (t, p) => lumenDist(p[0], p[1], p[2], G.A, 'AOROOT') < 0, 6));
  };
  m.sov = Math.max(aoRun(GED, 0.7) || 0, aoRun(GED, 0.85) || 0, aoRun(GED, 1.0) || 0) * 10 || null;
  m.stj = (aoRun(GED, 2.1) || 0) * 10 || null;
  {
    const asc = GED.A.ao.asc, ua = unit(sub(asc.b, asc.a));
    const c = add(asc.a, mul(ua, 1.4));
    const pa = unit(sub(perpAo, mul(ua, dot(perpAo, ua))));
    const r = run(GED, sub(c, mul(pa, 3)), pa, (t, p) => lumenDist(p[0], p[1], p[2], GED.A, 'AOROOT') < 0, 6);
    m.asc = r ? len(r) * 10 : null;
  }

  // LA volume at end-systole (voxel integration of the LA blood pool)
  {
    let n = 0; const h = 0.2;
    for (let x = -3; x <= 5; x += h) for (let y = -3; y <= 7; y += h) for (let z = -5; z <= 4; z += h) {
      if (lumenDist(x, y, z, GES.A, 'LA') < 0) n++;
    }
    m.lavol = n * h * h * h;
  }

  // A4C plane: RA area (ES), RV basal diameter and areas (ED / RV end-systole)
  const a4c = TTE_VIEWS.A4C.probe();
  const onA4c = (G, track) => {
    const o = track ? TTE_VIEWS.A4C.track(G.A) : [0, 0, 0];
    return { pos: add(a4c.pos, o), lat: a4c.lat, dir: a4c.dir, normal: a4c.normal };
  };
  const areaIn = (G, code, c, excl) => {
    const P = onA4c(G, true); let n = 0; const h = 0.1;
    const c0 = sub(c, mul(P.normal, dot(sub(c, P.pos), P.normal)));
    for (let i = -6; i <= 6; i += h) for (let j = -6; j <= 6; j += h) {
      const p = add(add(c0, mul(P.lat, i)), mul(P.dir, j));
      if (excl && (inCapsule(p, G.A.svc) || inCapsule(p, G.A.ivc))) continue;
      if (classify(p[0], p[1], p[2], G, path).tissue === code) n++;
    }
    return n * h * h;
  };
  m.raArea = areaIn(GES, TISSUE.RA, GES.A.ra.c, true);
  m.rvedArea = areaIn(GED, TISSUE.RV, GED.A.rv.c);
  const rvesArea = areaIn(GRES, TISSUE.RV, GRES.A.rv.c);
  m.rvfac = m.rvedArea > 0 ? (1 - rvesArea / m.rvedArea) * 100 : null;
  {
    const P = onA4c(GED, true);
    const w = inPlane(sub(GED.A.valves.mitral.c, GED.A.valves.tricuspid.c), P.normal);
    const onP = (p) => sub(p, mul(P.normal, dot(sub(p, P.pos), P.normal)));
    const c = onP([-3.5, GED.A.valves.tricuspid.c[1] - 1.0, 0]);
    const r = run(GED, sub(c, mul(w, 5)), w, isT(TISSUE.RV, TISSUE.VALVE), 9, 0, 0, path);
    m.rvd1 = r ? len(r) * 10 : null;
  }
  // TAPSE: lateral tricuspid annular excursion at the RV's end-systole
  { const a = GED.A.valves.tricuspid.lat, b = GRES.A.valves.tricuspid.lat; m.tapse = Math.hypot(...sub(a, b)) * 10; }

  // Doppler and haemodynamics
  const H = GED.hemo;
  m.avVmax = H.vAoPeak;
  m.avPeakGrad = 4 * H.vAoPeak ** 2;
  // mean gradient: the mean LV-aortic pressure difference over ejection (what a
  // traced CW envelope's mean 4v^2 estimates); AVA by continuity (SV / VTI)
  {
    let dp = 0, s1 = 0, n = 0;
    for (let i = 0; i < 240; i++) {
      const h = hemodynamics(i / 240, path);
      if (h.Qao > 1) { dp += Math.max(0, h.Plv - h.Pao); s1 += h.Qao / sum.QaoMax; n++; }
    }
    if (n) {
      m.avMeanGrad = dp / n;
      const vti = H.vAoPeak * 100 * (s1 / n) * sum.ejectTime;   // cm
      m.ava = vti > 0 ? sum.forwardSV / vti : null;
    }
  }
  m.mvEVel = H.vMitPeak;
  if (path.mitralStenosis) { m.mva = sum.mvaHakki; m.mvMeanGrad = sum.meanGradMV; }
  if (path.mr || sum.regurgVol > 1) { m.mrRf = sum.regurgFraction * 100; m.mrRvol = sum.regurgVol; }
  m.trVmax = H.vTRPeak || null;
  m.ivc = 2 * GED.A.ivc.r * 10;
  m.ivcCollapse = GED.A.ivc.collapse != null ? GED.A.ivc.collapse * 100 : null;
  m.rap = rapEstimate(GED.A.ivc);
  if (m.trVmax) m.pasp = 4 * m.trVmax ** 2 + m.rap;

  // pericardial effusion: echo-free depth behind the LV posterior wall (PLAX, ED)
  if (path.effusion) {
    const c = [0, GED.A.lv.M[1] - 1.3, 0];
    const r = run(GED, c, perpAx, isT(TISSUE.PERICARDIUM), 8, 0, 0, path);
    m.effusion = r ? len(r) : null;
  }
  m.shunt = shuntLabel(path);
  m.rwma = path.rwma || null;
  m.grade = gradeOf(path);
  m.esPhase = ES; m.rvEsPhase = rvEs;
  // linear measurements are reported (and graded) in whole millimetres, as a
  // sonographer reads them off the caliper
  for (const k of ['lvidd', 'lvids', 'ivsd', 'pwd', 'sov', 'stj', 'asc', 'rvd1', 'tapse', 'mapse', 'ivc']) if (m[k] != null) m[k] = Math.round(m[k]);
  _cache.set(key, m);
  return m;
}

// ---- report model ------------------------------------------------------------------
const f0 = (v) => (v == null || !Number.isFinite(v) ? '—' : v.toFixed(0));
const f1 = (v) => (v == null || !Number.isFinite(v) ? '—' : v.toFixed(1));
const f2 = (v) => (v == null || !Number.isFinite(v) ? '—' : v.toFixed(2));
const RWMA_WORDS = {
  septal: 'septal and apical hypokinesis (LAD territory)', anteroseptal: 'anteroseptal hypokinesis (LAD territory)',
  anterior: 'anterior hypokinesis (LAD territory)', lateral: 'lateral hypokinesis (circumflex territory)',
  inferior: 'inferior hypokinesis (RCA territory)', inferolateral: 'inferolateral hypokinesis (circumflex territory)',
  posterior: 'inferolateral hypokinesis (circumflex territory)', apical: 'apical hypokinesis (LAD territory)',
};

// Patient defaults: the model is one idealised adult, 175 cm / 75 kg.
export const DEFAULT_PATIENT = { sex: 'male', age: 45, height: 175, weight: 75 };

// Grades every measurement and returns { patient, rows (grouped tables),
// findings (sentences per section), conclusion (list), learner } for the UI.
export function buildReport(path, patient = DEFAULT_PATIENT, learner = {}) {
  const m = measureModel(path);
  const sx = patient.sex === 'female' ? 'female' : 'male';
  const B = bsa(patient.height, patient.weight), Hm = patient.height / 100;
  const ix = (v) => (v == null ? null : v / B);
  const row = (label, value, unit, g, idx) => ({ label, value, unit, idx, ref: g ? g.range : '', level: g ? g.level : 'na', word: g ? g.word : '' });

  const lvedvi = ix(m.lvedv), lvesvi = ix(m.lvesv), lvmi = ix(m.lvmass), lavi = ix(m.lavol);
  const g = {
    lvidd: gradeGraded('lvidd', m.lvidd, sx), lvids: gradeGraded('lvids', m.lvids, sx),
    lvedvi: gradeGraded('lvedvi', lvedvi, sx), lvesvi: gradeGraded('lvesvi', lvesvi, sx),
    lvmi: gradeGraded('lvmi', lvmi, sx), ivsd: gradeNormal('ivsd', m.ivsd, sx), pwd: gradeNormal('pwd', m.pwd, sx),
    lvef: gradeLvef(m.lvef), lavi: gradeLavi(lavi),
    sovh: gradeNormal('sovh', m.sov / Hm, sx), stjh: gradeNormal('stjh', m.stj / Hm, sx), asch: gradeNormal('asch', m.asc / Hm, sx),
    rvd1: gradeNormal('rvd1', m.rvd1, sx), tapse: gradeLimit('tapse', m.tapse, sx), rvfac: gradeLimit('rvfac', m.rvfac, sx),
    rvedai: gradeLimit('rvedai', ix(m.rvedArea), sx), raArea: gradeLimit('raArea', m.raArea, sx), raAreai: gradeLimit('raAreai', ix(m.raArea), sx),
  };
  const geom = lvGeometry(m.rwt, lvmi, sx);
  const as = path.aorticStenosis ? gradeAs(m.avVmax, m.avMeanGrad, m.ava) : null;
  const ms = path.mitralStenosis ? gradeMs(m.mva) : null;
  const mr = m.mrRf != null ? gradeRegurg(m.mrRf, m.mrRvol) : null;
  const ph = m.trVmax ? phProbability(m.trVmax) : null;
  const eff = m.effusion != null ? effusionSize(m.effusion) : null;
  const glsOk = m.gls != null && m.gls <= GLS_NORMAL;
  const eaOk = m.ea != null && m.ea >= EA_NORMAL[0] && m.ea <= EA_NORMAL[1];

  const rows = [
    { title: 'Left ventricle', rows: [
      row('IVSd', f0(m.ivsd), 'mm', g.ivsd), row('LVIDd', f0(m.lvidd), 'mm', g.lvidd),
      row('PWd', f0(m.pwd), 'mm', g.pwd), row('LVIDs', f0(m.lvids), 'mm', g.lvids),
      row('LVEDV', f0(m.lvedv), 'mL', g.lvedvi, `${f0(lvedvi)} mL/m²`), row('LVESV', f0(m.lvesv), 'mL', g.lvesvi, `${f0(lvesvi)} mL/m²`),
      row('LVEF', f0(m.lvef), '%', g.lvef),
      row('LV mass', f0(m.lvmass), 'g', g.lvmi, `${f0(lvmi)} g/m²`),
      row('RWT', f2(m.rwt), '', { range: '≤ 0.42', level: m.rwt > 0.42 ? 'mild' : 'normal', word: geom || '' }),
      row('GLS', f1(m.gls), '%', { range: `≤ ${GLS_NORMAL}`, level: glsOk ? 'normal' : 'abnormal', word: glsOk ? 'normal' : 'reduced' }),
      row('MAPSE', f0(m.mapse), 'mm', { range: '≥ 10', level: m.mapse >= 10 ? 'normal' : 'abnormal', word: m.mapse >= 10 ? 'normal' : 'reduced' }),
      row('E/A', f1(m.ea), '', { range: `${EA_NORMAL[0]}–${EA_NORMAL[1]}`, level: eaOk ? 'normal' : 'abnormal', word: eaOk ? 'normal' : (m.ea < EA_NORMAL[0] ? 'low' : 'high') }),
    ] },
    { title: 'Atria', rows: [
      row('LA volume', f0(m.lavol), 'mL', g.lavi, `${f0(lavi)} mL/m²`),
      row('RA area', f1(m.raArea), 'cm²', g.raArea, `${f1(ix(m.raArea))} cm²/m²`),
    ] },
    { title: 'Right ventricle', rows: [
      row('RVD1 (base)', f0(m.rvd1), 'mm', g.rvd1),
      row('RVED area', f1(m.rvedArea), 'cm²', g.rvedai, `${f1(ix(m.rvedArea))} cm²/m²`),
      row('RV FAC', f0(m.rvfac), '%', g.rvfac), row('TAPSE', f0(m.tapse), 'mm', g.tapse),
    ] },
    { title: 'Aorta', rows: [
      row('Sinus of Valsalva', f0(m.sov), 'mm', g.sovh, `${f1(m.sov / Hm)} mm/m`),
      row('Sino-tubular junction', f0(m.stj), 'mm', g.stjh, `${f1(m.stj / Hm)} mm/m`),
      row('Proximal ascending', f0(m.asc), 'mm', g.asch, `${f1(m.asc / Hm)} mm/m`),
    ] },
    { title: 'Doppler', rows: [
      row('AV Vmax', f1(m.avVmax), 'm/s', { range: '< 2.5', level: as ? as.level : (m.avVmax < 2.5 ? 'normal' : 'mild'), word: as ? as.word : '' }),
      row('AV mean gradient', f0(m.avMeanGrad), 'mmHg', null),
      row('AVA (continuity)', f2(m.ava), 'cm²', null),
      row('MV E velocity', f1(m.mvEVel), 'm/s', null),
      ...(m.mva != null ? [row('MVA (pressure half-time)', f2(m.mva), 'cm²', ms ? { range: '> 1.5', level: ms.level, word: ms.word } : null),
        row('MV mean gradient', f1(m.mvMeanGrad), 'mmHg', null)] : []),
      ...(m.mrRf != null ? [row('MR regurgitant fraction', f0(m.mrRf), '%', { range: '< 30', level: mr.level, word: mr.word }),
        row('MR regurgitant volume', f0(m.mrRvol), 'mL', null)] : []),
      ...(m.trVmax ? [row('TR Vmax', f1(m.trVmax), 'm/s', { range: '≤ 2.8', level: ph.level, word: `${ph.word} probability of PH` }),
        row('Est. PASP', f0(m.pasp), 'mmHg', null)] : []),
      row('IVC', f0(m.ivc), 'mm', { range: `≤ ${IVC_NORMAL_MAX * 10}`, level: m.ivc <= IVC_NORMAL_MAX * 10 ? 'normal' : 'abnormal', word: `RAP ~${m.rap} mmHg` }),
    ] },
  ];
  if (eff) rows.push({ title: 'Pericardium', rows: [row('Effusion depth (PLAX)', f1(m.effusion), 'cm', { range: 'none', level: eff.level, word: eff.word })] });

  // ---- findings (BSE order: LV, RV, atria, valves, aorta, pericardium, IVC) ----
  const F = [];
  const lvSize = g.lvidd.level === 'normal' && g.lvedvi.level === 'normal' ? 'Normal LV size' : `LV ${g.lvedvi.level !== 'normal' ? g.lvedvi.word : g.lvidd.word} (LVEDVi ${f0(lvedvi)} mL/m², LVIDd ${f0(m.lvidd)} mm)`;
  const wall = g.ivsd.level === 'normal' && g.pwd.level === 'normal' ? 'normal wall thickness' : `increased wall thickness (IVSd ${f0(m.ivsd)} mm, PWd ${f0(m.pwd)} mm)`;
  F.push({ section: 'Left ventricle', text: `${lvSize} with ${wall}${geom && geom !== 'normal geometry' ? `; ${geom}` : ''}. LVEF ${f0(m.lvef)}% (${g.lvef.word})${m.gls != null ? `, GLS ${f1(m.gls)}%` : ''}. ${m.rwma ? `Regional wall motion abnormality: ${RWMA_WORDS[m.rwma] || m.rwma}.` : 'No regional wall motion abnormality.'} Mitral inflow E/A ${f1(m.ea)}.` });
  const rvOk = g.rvd1.level === 'normal' && g.tapse.level === 'normal' && g.rvfac.level === 'normal';
  F.push({ section: 'Right ventricle', text: rvOk ? `Normal RV size and systolic function (RVD1 ${f0(m.rvd1)} mm, TAPSE ${f0(m.tapse)} mm, FAC ${f0(m.rvfac)}%).`
    : `RV ${g.rvd1.level !== 'normal' ? 'dilated' : 'normal size'} (RVD1 ${f0(m.rvd1)} mm); TAPSE ${f0(m.tapse)} mm, FAC ${f0(m.rvfac)}%${g.tapse.level !== 'normal' || g.rvfac.level !== 'normal' ? ' — impaired RV systolic function' : ''}.` });
  F.push({ section: 'Atria', text: `LA ${g.lavi.word} (LAVi ${f0(lavi)} mL/m²). RA ${g.raArea.level === 'normal' ? 'normal size' : 'dilated'} (area ${f1(m.raArea)} cm²).${m.shunt && path.asd ? ` ${cap(m.shunt.label)} with left-to-right shunt across the interatrial septum.` : ''}` });
  F.push({ section: 'Aortic valve', text: as ? `${cap(as.word)}: Vmax ${f1(m.avVmax)} m/s, mean gradient ${f0(m.avMeanGrad)} mmHg, AVA ${f2(m.ava)} cm² (continuity).` : `Trileaflet, opens well. Vmax ${f1(m.avVmax)} m/s. No aortic stenosis or regurgitation.` });
  const mvParts = [];
  if (ms) mvParts.push(`${cap(ms.word)}: MVA ${f2(m.mva)} cm², mean gradient ${f1(m.mvMeanGrad)} mmHg`);
  if (mr) mvParts.push(`${cap(mr.word)} mitral regurgitation (regurgitant fraction ${f0(m.mrRf)}%, volume ${f0(m.mrRvol)} mL)`);
  F.push({ section: 'Mitral valve', text: mvParts.length ? `${mvParts.join('. ')}.` : 'Structurally normal; no significant regurgitation or stenosis.' });
  const trWord = path.tr ? `${path.rvpo ? 'Functional' : cap(m.grade)} tricuspid regurgitation` : 'Tricuspid regurgitation';
  F.push({ section: 'Tricuspid & pulmonary valves', text: m.trVmax ? `${trWord}, TR Vmax ${f1(m.trVmax)} m/s; estimated PASP ${f0(m.pasp)} mmHg (RAP ${m.rap} mmHg). ${cap(ph.word)} echocardiographic probability of pulmonary hypertension.` : 'Normal tricuspid and pulmonary valves; insufficient TR to estimate PASP.' });
  if (m.shunt && path.vsd) F.push({ section: 'Ventricular septum', text: `${cap(m.shunt.label)} with left-to-right flow.` });
  const aoOk = [g.sovh, g.stjh, g.asch].every((x) => x.level === 'normal');
  F.push({ section: 'Aorta', text: aoOk ? `Normal aortic root and proximal ascending aorta (sinus ${f0(m.sov)} mm, STJ ${f0(m.stj)} mm, ascending ${f0(m.asc)} mm).` : `Aortic dimensions: sinus ${f0(m.sov)} mm, STJ ${f0(m.stj)} mm, ascending ${f0(m.asc)} mm — check against height-indexed limits.` });
  F.push({ section: 'Pericardium', text: eff ? `${cap(eff.word)} pericardial effusion (${f1(m.effusion)} cm posteriorly in PLAX).` : 'No pericardial effusion.' });
  F.push({ section: 'IVC', text: `IVC ${f0(m.ivc)} mm${m.ivcCollapse != null ? ` with ${f0(m.ivcCollapse)}% inspiratory collapse` : ''}; estimated RA pressure ${m.rap} mmHg.` });

  // ---- conclusion: abnormal findings first, most important first ----
  const C = [];
  if (as) C.push(`${cap(as.word)} (Vmax ${f1(m.avVmax)} m/s, mean gradient ${f0(m.avMeanGrad)} mmHg, AVA ${f2(m.ava)} cm²).`);
  if (ms) C.push(`${cap(ms.word)} (MVA ${f2(m.mva)} cm², mean gradient ${f1(m.mvMeanGrad)} mmHg).`);
  if (mr && path.mr) C.push(`${cap(mr.word)} mitral regurgitation (RF ${f0(m.mrRf)}%, RVol ${f0(m.mrRvol)} mL).`);
  if (g.lvef.level !== 'normal') C.push(`${cap(g.lvef.word)} LV systolic function, LVEF ${f0(m.lvef)}%.`);
  if (g.lvedvi.level !== 'normal' && g.lvedvi.level !== 'low') C.push(`${cap(g.lvedvi.word)} LV (LVEDVi ${f0(lvedvi)} mL/m²).`);
  if (m.rwma) C.push(`Regional wall motion abnormality: ${RWMA_WORDS[m.rwma] || m.rwma}.`);
  if (geom && geom !== 'normal geometry') C.push(`LV ${geom}.`);
  if (mr && !path.mr) C.push(`${cap(mr.word)} functional mitral regurgitation (RF ${f0(m.mrRf)}%).`);
  if (g.lavi.level !== 'normal') C.push(`${cap(g.lavi.word)} left atrium (LAVi ${f0(lavi)} mL/m²).`);
  if (m.shunt) C.push(`${cap(m.shunt.label)}.`);
  if (ph && ph.level !== 'normal') C.push(`${cap(ph.word)} echocardiographic probability of pulmonary hypertension (TR Vmax ${f1(m.trVmax)} m/s, PASP ~${f0(m.pasp)} mmHg).`);
  if (path.tr && !path.rvpo) C.unshift(`${cap(m.grade)} tricuspid regurgitation (TR Vmax ${f1(m.trVmax)} m/s).`);
  if (!rvOk) {
    const rvBits = [];
    if (g.rvd1.level !== 'normal') rvBits.push(`RV dilatation (RVD1 ${f0(m.rvd1)} mm)`);
    if (g.tapse.level !== 'normal' || g.rvfac.level !== 'normal') rvBits.push(`impaired RV systolic function (TAPSE ${f0(m.tapse)} mm, FAC ${f0(m.rvfac)}%)`);
    C.push(`${cap(rvBits.join(' with '))}.`);
  }
  if (eff) C.push(`${cap(eff.word)} pericardial effusion (${f1(m.effusion)} cm).`);
  if (!C.length) C.push('Structurally normal heart. Normal biventricular size and systolic function. No significant valve disease.');

  // learner measurements, compared with the model's own values
  const L = [];
  if (learner.ef != null && Number.isFinite(learner.ef)) L.push({ label: 'Your Simpson EF', value: `${f0(learner.ef)}%`, truth: `${f0(m.lvef)}%`, diff: learner.ef - m.lvef, tol: 6 });
  if (learner.distCm != null && Number.isFinite(learner.distCm)) L.push({ label: 'Your caliper', value: `${f1(learner.distCm * 10)} mm`, truth: '', diff: null });

  return { patient: { ...patient, bsa: B }, rows, findings: F, conclusion: C, learner: L, source: REF_SOURCE, m };
}
const cap = (s) => (s ? s[0].toUpperCase() + s.slice(1) : s);

// ---- templates ---------------------------------------------------------------------
// Report scaffolds a learner completes in their own words. {KEY} placeholders
// are filled from the measurements; [square brackets] mark judgements the
// learner must make from the images (left for them on purpose).
export const TEMPLATES = {
  blank: { name: 'Blank structured report', text:
`Indication: [ ]
Rhythm: [ ]   Image quality: [ ]

Left ventricle: [size], [wall thickness]. LVEF {LVEF}% by [method]. [Regional wall motion]. [Diastolic function].
Right ventricle: [size]. TAPSE {TAPSE} mm, FAC {RVFAC}%.
Atria: LAVi {LAVI} mL/m². RA area {RAAREA} cm².
Aortic valve: [morphology]. Vmax {AVVMAX} m/s.
Mitral valve: [morphology], [regurgitation/stenosis].
Tricuspid / pulmonary: [ ]. TR Vmax [ ] m/s.
Aorta: sinus {SOV} mm, STJ {STJ} mm, ascending {ASC} mm.
Pericardium: [ ].   IVC: {IVC} mm, [collapse].

Conclusion:
1. [ ]` },
  normal: { name: 'Normal study', text:
`Left ventricle: normal size (LVIDd {LVIDD} mm, LVEDVi {LVEDVI} mL/m²) and wall thickness (IVSd {IVSD} mm, PWd {PWD} mm). LVEF {LVEF}% (Simpson's biplane). No regional wall motion abnormality. Normal diastolic function (E/A {EA}).
Right ventricle: normal size and systolic function (TAPSE {TAPSE} mm, FAC {RVFAC}%).
Atria: normal LA (LAVi {LAVI} mL/m²) and RA size.
Valves: structurally and functionally normal. No significant regurgitation or stenosis.
Aorta: normal root and proximal ascending aorta (sinus {SOV} mm).
Pericardium: no effusion. IVC {IVC} mm with > 50% collapse (RAP ~3 mmHg).

Conclusion:
1. Structurally normal heart with normal biventricular size and systolic function (LVEF {LVEF}%).` },
  as: { name: 'Aortic stenosis', text:
`Aortic valve: [calcified / restricted cusps, number of cusps]. Vmax {AVVMAX} m/s, mean gradient {AVMEAN} mmHg, AVA {AVA} cm² (continuity, LVOT [ ] mm).
Left ventricle: wall thickness IVSd {IVSD} mm, PWd {PWD} mm — [geometry]. LVEF {LVEF}%.
Check: [low-flow state? SVi = SV / BSA]  [aortic regurgitation]  [ascending aorta {ASC} mm].

Conclusion:
1. [Grade] aortic stenosis (Vmax {AVVMAX} m/s, mean gradient {AVMEAN} mmHg, AVA {AVA} cm²).
2. LV [geometry], LVEF {LVEF}%.` },
  mr: { name: 'Mitral regurgitation', text:
`Mitral valve: [leaflet morphology: prolapse / flail / restriction, which scallop]. Jet [central / eccentric, direction].
Quantification: regurgitant fraction {MRRF}%, regurgitant volume {MRRVOL} mL. [Vena contracta, PISA EROA, pulmonary-vein flow].
Left ventricle: LVEDVi {LVEDVI} mL/m², LVESVi {LVESVI} mL/m², LVEF {LVEF}%.
Left atrium: LAVi {LAVI} mL/m².

Conclusion:
1. [Grade] [primary / secondary] mitral regurgitation (RF {MRRF}%, RVol {MRRVOL} mL).
2. LV [size], LVEF {LVEF}%; LA [size].` },
  ms: { name: 'Mitral stenosis', text:
`Mitral valve: [rheumatic thickening, commissural fusion, doming anterior leaflet, subvalvular involvement].
MVA {MVA} cm² ([planimetry / pressure half-time]); mean gradient {MVMEAN} mmHg at HR [ ].
Left atrium: LAVi {LAVI} mL/m². [LAA thrombus?]
Pulmonary pressure: TR Vmax [ ] m/s, PASP [ ] mmHg.

Conclusion:
1. [Grade] rheumatic mitral stenosis (MVA {MVA} cm², mean gradient {MVMEAN} mmHg).
2. [LA dilatation]  [pulmonary hypertension probability].` },
  dcm: { name: 'Dilated cardiomyopathy', text:
`Left ventricle: [dilated], LVIDd {LVIDD} mm, LVEDVi {LVEDVI} mL/m². Global hypokinesis. LVEF {LVEF}% (Simpson's biplane), GLS {GLS}%.
Wall thickness IVSd {IVSD} mm, PWd {PWD} mm. [Thrombus at apex?]
Mitral valve: [functional MR from annular dilatation / tenting], RF {MRRF}%.
Right ventricle: TAPSE {TAPSE} mm, FAC {RVFAC}%.
Atria: LAVi {LAVI} mL/m².

Conclusion:
1. [Dilated] LV with [severely impaired] systolic function, LVEF {LVEF}%.
2. [Grade] functional mitral regurgitation.  3. [RV function]` },
  rwma: { name: 'Regional wall motion abnormality', text:
`Left ventricle: normal size. Regional wall motion: [segments, AHA 17-segment model, hypokinetic / akinetic / dyskinetic] — [coronary territory].
LVEF {LVEF}% (Simpson's biplane). GLS {GLS}% [pattern on bullseye].
[Complications: thrombus, aneurysm, VSD, MR, effusion]

Conclusion:
1. Regional wall motion abnormality in the [ ] territory, LVEF {LVEF}%.` },
  effusion: { name: 'Pericardial effusion', text:
`Pericardium: [global / localised] pericardial effusion, maximal diastolic depth {EFFUSION} cm [location].
Haemodynamics: [RA / RV diastolic collapse]  [respiratory variation in mitral / tricuspid inflow]  IVC {IVC} mm [plethoric?].
[Fibrin strands / echogenic content]

Conclusion:
1. [Size] pericardial effusion ({EFFUSION} cm), [with / without] echocardiographic features of tamponade.` },
  vsd: { name: 'Ventricular septal defect', text:
`Ventricular septum: [location: perimembranous / muscular] defect, [size] mm, left-to-right flow, peak velocity [ ] m/s (LV–RV gradient [ ] mmHg).
Left ventricle: LVEDVi {LVEDVI} mL/m² [volume load?]. LVEF {LVEF}%.
Right ventricle / pulmonary pressure: [ ].

Conclusion:
1. [Restrictive] [type] VSD with left-to-right shunt.` },
  asd: { name: 'Atrial septal defect', text:
`Interatrial septum: [secundum / primum / sinus venosus] defect, [size] mm, [rims], left-to-right flow on colour Doppler.
Right heart: RVD1 {RVD1} mm, RVED area {RVEDAREA} cm², RA area {RAAREA} cm² [volume load?]. TAPSE {TAPSE} mm.
Pulmonary pressure: TR Vmax [ ] m/s.

Conclusion:
1. [Type] atrial septal defect with left-to-right shunt.  2. [RV dilatation]` },
  tr: { name: 'Tricuspid regurgitation', text:
`Tricuspid valve: [leaflet morphology, annular diameter]. TR [grade] (vena contracta [ ] mm, hepatic-vein flow [ ]).
TR Vmax {TRVMAX} m/s; estimated PASP {PASP} mmHg (RAP from IVC {IVC} mm).
Right heart: RVD1 {RVD1} mm, RA area {RAAREA} cm², TAPSE {TAPSE} mm.

Conclusion:
1. [Grade] tricuspid regurgitation. TR Vmax {TRVMAX} m/s — [probability] of pulmonary hypertension.` },
  phtn: { name: 'Pulmonary hypertension', text:
`Right ventricle: [dilated], RVD1 {RVD1} mm, RVED area {RVEDAREA} cm². [Septal flattening / D-shaped LV in PSAX]. TAPSE {TAPSE} mm, FAC {RVFAC}%.
TR Vmax {TRVMAX} m/s; estimated PASP {PASP} mmHg (RAP {RAP} mmHg).
Additional PH signs: [RV/LV basal ratio > 1, PA diameter > 25 mm, RVOT acceleration time < 105 ms, IVC > 21 mm].
RA area {RAAREA} cm².

Conclusion:
1. [High] echocardiographic probability of pulmonary hypertension (TR Vmax {TRVMAX} m/s, PASP ~{PASP} mmHg).
2. [RV dilatation and systolic function].` },
};
// the template that matches each pathology key in the UI
export const TEMPLATE_FOR = { normal: 'normal', as: 'as', mr: 'mr', ms: 'ms', dcm: 'dcm', rwma: 'rwma', effusion: 'effusion', vsd: 'vsd', asd: 'asd', tr: 'tr', phtn: 'phtn' };

export function fillTemplate(key, rep) {
  const t = TEMPLATES[key] || TEMPLATES.blank, m = rep.m, B = rep.patient.bsa;
  const V = {
    LVEF: f0(m.lvef), LVIDD: f0(m.lvidd), IVSD: f0(m.ivsd), PWD: f0(m.pwd), LVEDVI: f0(m.lvedv / B), LVESVI: f0(m.lvesv / B),
    EA: f1(m.ea), GLS: f1(m.gls), TAPSE: f0(m.tapse), RVFAC: f0(m.rvfac), RVD1: f0(m.rvd1), RVEDAREA: f1(m.rvedArea),
    LAVI: f0(m.lavol / B), RAAREA: f1(m.raArea), AVVMAX: f1(m.avVmax), AVMEAN: f0(m.avMeanGrad), AVA: f2(m.ava),
    MVA: f2(m.mva), MVMEAN: f1(m.mvMeanGrad), MRRF: f0(m.mrRf), MRRVOL: f0(m.mrRvol), TRVMAX: f1(m.trVmax), PASP: f0(m.pasp),
    RAP: f0(m.rap), IVC: f0(m.ivc), SOV: f0(m.sov), STJ: f0(m.stj), ASC: f0(m.asc), EFFUSION: f1(m.effusion),
  };
  return t.text.replace(/\{([A-Z0-9]+)\}/g, (s, k) => (k in V ? V[k] : s));
}

// ---- plain-text rendering (for copy / download) --------------------------------------
export function reportText(rep, writeUp = '') {
  const p = rep.patient, out = [];
  out.push('TRANSTHORACIC ECHOCARDIOGRAM — EchoSim teaching report (simulated patient, not for clinical use)');
  out.push(`Patient: ${p.sex}, ${p.age} y, ${p.height} cm, ${p.weight} kg, BSA ${p.bsa.toFixed(2)} m²`);
  out.push(`Reference intervals: ${rep.source}`);
  out.push('');
  for (const grp of rep.rows) {
    out.push(grp.title.toUpperCase());
    for (const r of grp.rows) {
      const v = `${r.value}${r.unit ? ' ' + r.unit : ''}`;
      out.push(`  ${r.label.padEnd(26)}${v.padEnd(12)}${(r.idx || '').padEnd(16)}${r.ref ? `ref ${r.ref}`.padEnd(14) : ''.padEnd(14)}${r.word || ''}`);
    }
    out.push('');
  }
  out.push('FINDINGS');
  for (const f of rep.findings) out.push(`  ${f.section}: ${f.text}`);
  out.push('');
  out.push('CONCLUSION');
  rep.conclusion.forEach((c, i) => out.push(`  ${i + 1}. ${c}`));
  if (rep.learner.length) {
    out.push('');
    out.push('LEARNER MEASUREMENTS');
    for (const l of rep.learner) out.push(`  ${l.label}: ${l.value}${l.truth ? ` (model ${l.truth})` : ''}`);
  }
  if (writeUp && writeUp.trim()) { out.push(''); out.push('LEARNER REPORT'); out.push(writeUp.trim()); }
  return out.join('\n');
}
