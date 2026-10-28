// Reference intervals and severity grading for the structured echo report.
//
// Chamber quantification follows the British Society of Echocardiography 2020
// guideline and its poster "Echocardiography: Reference Intervals & Functional
// Assessment" (Harkness A, Ring L, Augustine DX, Oxborough D, Robinson S,
// Sharma V. Echo Res Pract 2020;7(1):G1-G18; right heart: Zaidi A et al. Echo
// Res Pract 2020;7(1):G19-G41). Those intervals are sex-specific and graded
// (normal / mild / moderate / severe), are derived from a largely Caucasian
// population, and are not valid for 3D measurements, athletes or pregnancy.
//
// Valve and haemodynamic grading uses the thresholds the BSE adopts from the
// ASE/EACVI valve guidelines (aortic stenosis: Baumgartner 2017; mitral
// stenosis: Baumgartner 2009; regurgitation: Zoghbi 2017; pulmonary
// hypertension probability: ESC/ERS 2022; RA pressure: Rudski 2010).
//
// Each graded entry lists band upper bounds in ascending order; a value at or
// below the first bound is normal. Units are those used in the report (mm for
// linear LV dimensions, mL/m2 for indexed volumes, g/m2 for indexed mass).

export const REF_SOURCE = 'BSE 2020 (Harkness et al., Echo Res Pract 2020;7:G1-G18)';

// ---- graded intervals: [normal max, mild max, moderate max]; above = severe ----
// Values are { male: {lo, bands}, female: {lo, bands} }; lo is the lower normal
// limit (a value below it is reported as "below normal").
export const GRADED = {
  lvidd:  { label: 'LVID diastole', unit: 'mm', male: { lo: 37, bands: [56, 61, 65] }, female: { lo: 35, bands: [51, 55, 59] } },
  lvids:  { label: 'LVID systole', unit: 'mm', male: { lo: 22, bands: [41, 45, 50] }, female: { lo: 20, bands: [37, 42, 46] } },
  lvedvi: { label: 'LVEDV indexed', unit: 'mL/m²', male: { lo: 30, bands: [79, 91, 103] }, female: { lo: 29, bands: [70, 81, 91] } },
  lvesvi: { label: 'LVESV indexed', unit: 'mL/m²', male: { lo: 9, bands: [31, 36, 42] }, female: { lo: 8, bands: [27, 32, 37] } },
  lvmi:   { label: 'LV mass indexed', unit: 'g/m²', words: 'increased', male: { lo: 40, bands: [110, 127, 145] }, female: { lo: 33, bands: [99, 115, 131] } },
};
// the severity word for each band index of a graded interval (dimensions and
// volumes are "dilated"; mass is "increased")
export const BAND_WORDS = ['normal', 'mildly dilated', 'moderately dilated', 'severely dilated'];
const bandWords = (g) => (g.words === 'increased' ? ['normal', 'mildly increased', 'moderately increased', 'severely increased'] : BAND_WORDS);

// ---- simple (two-sided) normal intervals, sex-specific -----------------------
export const NORMAL = {
  ivsd:   { label: 'IVS diastole', unit: 'mm', male: [6, 12], female: [5, 11] },
  pwd:    { label: 'PW diastole', unit: 'mm', male: [6, 12], female: [6, 12] },
  lvedv:  { label: 'LVEDV', unit: 'mL', male: [53, 156], female: [46, 121] },
  lvesv:  { label: 'LVESV', unit: 'mL', male: [15, 62], female: [13, 47] },
  lvmass: { label: 'LV mass', unit: 'g', male: [72, 219], female: [51, 173] },
  // aortic dimensions: inner edge to inner edge at end-diastole, indexed to HEIGHT
  sovh:   { label: 'Sinus of Valsalva / height', unit: 'mm/m', male: [13.8, 21.8], female: [13.1, 20.7] },
  stjh:   { label: 'Sino-tubular junction / height', unit: 'mm/m', male: [11.4, 18.6], female: [11.0, 17.8] },
  asch:   { label: 'Proximal ascending aorta / height', unit: 'mm/m', male: [11.5, 19.9], female: [11.4, 19.8] },
  // right heart
  rvotPlax: { label: 'RVOT proximal (PLAX)', unit: 'mm', male: [25, 43], female: [22, 40] },
  rvd1:   { label: 'RVD1 base', unit: 'mm', male: [26, 47], female: [22, 43] },
  rvd2:   { label: 'RVD2 mid', unit: 'mm', male: [19, 42], female: [17, 35] },
  rvd3:   { label: 'RVD3 length', unit: 'mm', male: [55, 87], female: [51, 80] },
};
// one-sided limits: { max } = upper limit of normal, { min } = lower limit
export const LIMIT = {
  rvedai: { label: 'RVED area indexed', unit: 'cm²/m²', male: { max: 13.6 }, female: { max: 12.6 } },
  rvfac:  { label: 'RV fractional area change', unit: '%', male: { min: 30 }, female: { min: 35 } },
  tapse:  { label: 'TAPSE', unit: 'mm', male: { min: 17 }, female: { min: 17 } },
  rvs:    { label: "RV s'", unit: 'cm/s', male: { min: 9 }, female: { min: 9 } },
  raArea: { label: 'RA area', unit: 'cm²', male: { max: 22 }, female: { max: 19 } },
  raAreai:{ label: 'RA area indexed', unit: 'cm²/m²', male: { max: 11 }, female: { max: 11 } },
};
// mean (septal + lateral) mitral annular s' lower limits by age band (cm/s)
export const MITRAL_S = [[40, 6.4], [60, 5.7], [Infinity, 4.9]];

// ---- LV function and LA size: BSE wording ----------------------------------
// The BSE asks for LVEF to be quoted as a number and the terms "mild" and
// "moderate" LV dysfunction to be avoided.
export const LVEF_BANDS = [
  { min: 55, word: 'normal', level: 'normal' },
  { min: 50, word: 'borderline low', level: 'mild' },
  { min: 36, word: 'impaired', level: 'moderate' },
  { min: -Infinity, word: 'severely impaired', level: 'severe' },
];
export const LAVI_BANDS = [
  { max: 34, word: 'normal', level: 'normal' },      // < 34
  { max: 38, word: 'borderline', level: 'mild' },     // 34-38
  { max: Infinity, word: 'dilated', level: 'severe' },
];
// LV geometry: relative wall thickness 2*PWd/LVIDd against indexed mass
export const RWT_CUTOFF = 0.42;

// ---- valves and haemodynamics ------------------------------------------------
// Aortic stenosis (peak velocity m/s, mean gradient mmHg, AVA cm2)
export const AS_GRADES = [
  { word: 'aortic sclerosis without stenosis', vmax: 2.5, level: 'normal' },
  { word: 'mild aortic stenosis', vmax: 3.0, mean: 20, ava: 1.5, level: 'mild' },
  { word: 'moderate aortic stenosis', vmax: 4.0, mean: 40, ava: 1.0, level: 'moderate' },
  { word: 'severe aortic stenosis', vmax: Infinity, level: 'severe' },
];
// Mitral stenosis by valve area (cm2): > 1.5 mild, 1.0-1.5 moderate, < 1.0 severe
export const MS_GRADES = [
  { min: 1.5, word: 'mild mitral stenosis', level: 'mild' },
  { min: 1.0, word: 'moderate mitral stenosis', level: 'moderate' },
  { min: -Infinity, word: 'severe mitral stenosis', level: 'severe' },
];
// Regurgitation by regurgitant fraction (%) and volume (mL)
export const REGURG_GRADES = [
  { rfMax: 30, rvolMax: 30, word: 'mild', level: 'mild' },
  { rfMax: 50, rvolMax: 60, word: 'moderate', level: 'moderate' },
  { rfMax: Infinity, rvolMax: Infinity, word: 'severe', level: 'severe' },
];
// Echocardiographic probability of pulmonary hypertension from peak TR velocity
export const PH_PROBABILITY = [
  { max: 2.8, word: 'low', level: 'normal' },
  { max: 3.4, word: 'intermediate', level: 'moderate' },
  { max: Infinity, word: 'high', level: 'severe' },
];
// Pericardial effusion size by maximal diastolic echo-free depth (cm)
export const EFFUSION_GRADES = [
  { max: 1.0, word: 'small', level: 'mild' },
  { max: 2.0, word: 'moderate', level: 'moderate' },
  { max: Infinity, word: 'large', level: 'severe' },
];
export const GLS_NORMAL = -18;    // % (more negative is normal)
export const EA_NORMAL = [0.8, 2.0];
export const IVC_NORMAL_MAX = 2.1; // cm

// ---- helpers -------------------------------------------------------------------
// Body surface area (Mosteller), m2, from height (cm) and weight (kg).
export const bsa = (hCm, wKg) => Math.sqrt((hCm * wKg) / 3600);

// Grade against a graded interval: { level, word, range } where level is one of
// normal / mild / moderate / severe / low (below the lower normal limit).
export function gradeGraded(key, v, sex = 'male') {
  const g = GRADED[key], s = g[sex];
  const range = `${s.lo}–${s.bands[0]}`;
  if (v == null || !Number.isFinite(v)) return { level: 'na', word: '—', range };
  if (v < s.lo) return { level: 'low', word: 'below normal', range };
  const i = s.bands.findIndex((b) => v <= b);
  const k = i < 0 ? 3 : i;
  return { level: ['normal', 'mild', 'moderate', 'severe'][k], word: bandWords(g)[k], range };
}
export function gradeNormal(key, v, sex = 'male') {
  const [lo, hi] = NORMAL[key][sex];
  const range = `${lo}–${hi}`;
  if (v == null || !Number.isFinite(v)) return { level: 'na', word: '—', range };
  if (v < lo) return { level: 'low', word: 'below normal', range };
  if (v > hi) return { level: 'high', word: 'above normal', range };
  return { level: 'normal', word: 'normal', range };
}
export function gradeLimit(key, v, sex = 'male') {
  const s = LIMIT[key][sex];
  const range = s.max != null ? `≤ ${s.max}` : `≥ ${s.min}`;
  if (v == null || !Number.isFinite(v)) return { level: 'na', word: '—', range };
  const ok = s.max != null ? v <= s.max : v >= s.min;
  return { level: ok ? 'normal' : 'abnormal', word: ok ? 'normal' : (s.max != null ? 'enlarged' : 'reduced'), range };
}
export function gradeLvef(ef) {
  if (ef == null || !Number.isFinite(ef)) return { level: 'na', word: '—', range: '≥ 55' };
  const b = LVEF_BANDS.find((x) => ef >= x.min);
  return { level: b.level, word: b.word, range: '≥ 55' };
}
export function gradeLavi(v) {
  if (v == null || !Number.isFinite(v)) return { level: 'na', word: '—', range: '< 34' };
  const b = LAVI_BANDS.find((x) => v < x.max || (x.max === 38 && v <= 38));
  return { level: b.level, word: b.word, range: '< 34' };
}
// LV geometry from relative wall thickness and indexed mass (normal mass
// limits from the graded table)
export function lvGeometry(rwt, lvmi, sex = 'male') {
  if (rwt == null || lvmi == null) return null;
  const massUp = lvmi > GRADED.lvmi[sex].bands[0], thick = rwt > RWT_CUTOFF;
  return massUp ? (thick ? 'concentric hypertrophy' : 'eccentric hypertrophy')
    : (thick ? 'concentric remodelling' : 'normal geometry');
}
export function gradeAs(vmax, mean, ava) {
  if (vmax == null) return null;
  return AS_GRADES.find((g) => vmax < g.vmax) || AS_GRADES[AS_GRADES.length - 1];
}
export function gradeMs(mva) {
  if (mva == null) return null;
  return MS_GRADES.find((g) => mva > g.min) || MS_GRADES[MS_GRADES.length - 1];
}
export function gradeRegurg(rfPct, rvol) {
  if (rfPct == null) return null;
  return REGURG_GRADES.find((g) => rfPct < g.rfMax && (rvol == null || rvol < g.rvolMax)) || REGURG_GRADES[REGURG_GRADES.length - 1];
}
export function phProbability(vtr) {
  if (vtr == null) return null;
  return PH_PROBABILITY.find((g) => vtr <= g.max);
}
export function effusionSize(depthCm) {
  if (depthCm == null) return null;
  return EFFUSION_GRADES.find((g) => depthCm < g.max) || EFFUSION_GRADES[EFFUSION_GRADES.length - 1];
}
