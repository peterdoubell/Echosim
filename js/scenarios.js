// Shared teaching-case definitions. The wall territory is used by anatomy.js;
// the injury vector is projected through the same VCG→12-lead transform as all
// other ECG morphology. These are representative infarction examples, not a
// claim that every RWMA produces acute ST elevation (see docs/SCENARIOS.md).
export const RWMA_SCENARIOS = {
  septal: {
    label: 'Septal', view: 'A4C', leads: ['V1', 'V2'],
    ecgText: 'Septal injury pattern · ST elevation most evident in V1–V2',
    cue: 'Compare reduced septal and apical thickening with the contracting lateral wall.',
    stVec: [-0.25, -0.12, -0.50],
  },
  anteroseptal: {
    label: 'Anteroseptal', view: 'PSAX', leads: ['V1', 'V2', 'V3', 'V4'],
    ecgText: 'Anteroseptal injury pattern · ST elevation in V1–V4',
    cue: 'Follow the anterior septum toward the apex; the opposite inferior/lateral wall keeps contracting.',
    stVec: [0.15, -0.12, -0.55],
  },
  anterior: {
    label: 'Anterior', view: 'A2C', leads: ['V3', 'V4'],
    ecgText: 'Anterior injury pattern · ST elevation most evident in V3–V4',
    cue: 'Compare the poorly thickening anterior wall with the inferior wall in the two-chamber view.',
    stVec: [0.35, -0.10, -0.35],
  },
  lateral: {
    label: 'Lateral', view: 'A4C', leads: ['I', 'aVL', 'V5', 'V6'],
    ecgText: 'Lateral injury pattern · ST elevation in I, aVL, V5–V6; inferior reciprocal depression',
    cue: 'Compare reduced lateral wall thickening with the contracting septum.',
    stVec: [0.50, -0.15, 0.30],
  },
  inferior: {
    label: 'Inferior', view: 'A2C', leads: ['II', 'III', 'aVF'],
    ecgText: 'Inferior injury pattern · ST elevation in II, III, aVF; reciprocal depression in I/aVL',
    cue: 'Compare reduced inferior wall thickening with the contracting anterior wall.',
    stVec: [-0.12, 0.45, 0.04],
  },
  inferolateral: {
    label: 'Inferolateral', view: 'PLAX', leads: ['II', 'III', 'aVF', 'I', 'V5', 'V6'],
    ecgText: 'Inferolateral injury pattern · inferior and lateral ST elevation; reciprocal depression in V1–V2',
    cue: 'Inspect the inferolateral wall in PLAX and PSAX; the opposing anteroseptal wall keeps contracting.',
    stVec: [0.50, 0.55, 0.40],
  },
};

// Use the same BSE 2020 terms as the report rather than importing a competing
// mild/moderate/severe EF scale. Targets describe the solved, illustrative cases.
export const DCM_SEVERITIES = {
  borderline: { label: 'Borderline low · EF ≈52%', targetEF: 52, minEF: 50, maxEF: 54, wallMul: 0.96 },
  impaired: { label: 'Impaired · EF ≈42%', targetEF: 42, minEF: 36, maxEF: 49, wallMul: 0.90 },
  severe: { label: 'Severely impaired · EF ≈27%', targetEF: 27, minEF: 0, maxEF: 35, wallMul: 0.82 },
};

export function dcmSeverityOf(path = {}) {
  return Object.hasOwn(DCM_SEVERITIES, path.dcmSeverity) ? path.dcmSeverity : 'severe';
}

export function rwmaScenario(path = {}) {
  if (!path.rwma) return null;
  const region = typeof path.rwma === 'object' ? path.rwma.region : path.rwma;
  return RWMA_SCENARIOS[region === 'posterior' ? 'inferolateral' : region] || RWMA_SCENARIOS.septal;
}
