// Learning pathway for BSE transthoracic echo (TTE) accreditation revision.
//
// Independent revision material: not produced by, endorsed by or affiliated
// with the British Society of Echocardiography. Learners should check the
// current syllabus, exam format and logbook rules at bsecho.org.
//
// Each module has objectives, key points (trusted HTML), practical tasks that
// load a case into the simulator, and a question bank. Questions are either
// plain multiple choice ({ q, options, answer, explain }) or case questions
// that also carry `case: { path, view, grade?, hideView? }`: the simulator
// plays that case with the diagnosis hidden while the learner answers, like
// the video questions in the written exam. `answer` indexes `options` as
// written; options are shuffled at display time.
//
// Numbers that also appear in the report come from js/references.js, so the
// pathway, the Measurements panel and the report always agree.

import { GRADED, NORMAL, LIMIT, LVEF_BANDS, RWT_CUTOFF, GLS_NORMAL } from './references.js';

const lvidM = GRADED.lvidd.male.bands[0], lvidF = GRADED.lvidd.female.bands[0];
const edviM = GRADED.lvedvi.male.bands[0];
const lvmiM = GRADED.lvmi.male;
const ivsF = NORMAL.ivsd.female, rvd1M = NORMAL.rvd1.male[1];
const tapse = LIMIT.tapse.male.min, rvs = LIMIT.rvs.male.min, facM = LIMIT.rvfac.male.min;
const efImp = LVEF_BANDS.find((b) => b.word === 'impaired').min;
const efBord = LVEF_BANDS.find((b) => b.word === 'borderline low').min;

export const PASS_MARK = 0.7;
export const EXAM = { questions: 40, minutes: 60 };
export const DISCLAIMER = 'Independent revision material, not produced by or affiliated with the British Society of Echocardiography. Check the current syllabus, exam format and logbook requirements at bsecho.org.';

export const MODULES = [
  {
    id: 'physics', title: 'Physics & instrumentation', icon: '〰️',
    objectives: [
      'Apply the Doppler equation and the effect of insonation angle',
      'Explain the Nyquist limit and how to reduce aliasing',
      'Relate frequency, resolution, penetration and frame rate',
      'Recognise common artefacts',
    ],
    learn: [
      'Doppler shift: <b>Δf = 2·f₀·v·cosθ / c</b>. A higher transmitted frequency gives a larger shift; at 20° velocity is underestimated by ~6%, at 60° by 50%.',
      'Pulsed-wave Doppler aliases when the shift exceeds the <b>Nyquist limit = PRF / 2</b>. Reduce aliasing by raising the scale, shifting the baseline, reducing depth, using a lower-frequency probe, or switching to <b>continuous-wave</b> (no aliasing, but no range resolution).',
      'Simplified Bernoulli: <b>ΔP = 4v²</b> (v in m/s, ΔP in mmHg). 4 m/s → 64 mmHg.',
      'Axial resolution improves with higher frequency / shorter pulses; penetration falls as frequency rises.',
      'Frame rate rises with a narrower sector, shallower depth and lower line density; colour Doppler lowers it.',
      'Tissue harmonic imaging receives at twice the transmitted frequency, reducing near-field clutter and side-lobe artefact.',
      'Artefacts: reverberation (equally spaced repeats behind a strong reflector), acoustic shadowing (behind calcium or prosthesis), mirror image, side lobes.',
      'Colour Doppler (BART): <b>B</b>lue <b>A</b>way, <b>R</b>ed <b>T</b>owards the transducer.',
    ],
    tasks: [
      { id: 'alias', text: 'Load mitral regurgitation in A4C. Lower the Nyquist (colour scale) in the Imaging panel and watch the jet alias into a mosaic; raise it again.', sim: { path: 'mr', view: 'A4C' } },
      { id: 'thi', text: 'Toggle tissue harmonic imaging (Imaging panel) in PLAX and compare near-field clutter and penetration.', sim: { path: 'normal', view: 'PLAX' } },
      { id: 'cw', text: 'Load severe aortic stenosis and read the peak velocity and gradient on the spectral trace; check ΔP ≈ 4v².', sim: { path: 'as', view: 'A5C', grade: 'severe' } },
    ],
    questions: [
      { q: 'For a given blood velocity, which change increases the Doppler frequency shift?', options: ['A higher transmitted frequency', 'An insonation angle closer to 90°', 'A lower transmitted frequency', 'A deeper sample volume'], answer: 0, explain: 'Δf = 2·f₀·v·cosθ / c: the shift is proportional to f₀ and to cosθ, which falls towards zero at 90°.' },
      { q: 'Pulsed-wave Doppler aliases when the Doppler shift exceeds:', options: ['Half the pulse repetition frequency', 'The pulse repetition frequency', 'Twice the pulse repetition frequency', 'The transmitted frequency'], answer: 0, explain: 'The Nyquist limit is PRF/2: above it the frequency cannot be sampled unambiguously and wraps around.' },
      { q: 'Which manoeuvre would make aliasing on pulsed-wave Doppler WORSE?', options: ['Switching to a higher-frequency transducer', 'Shifting the baseline', 'Reducing the sample depth', 'Switching to continuous-wave Doppler'], answer: 0, explain: 'A higher f₀ increases the Doppler shift for the same velocity. Reducing depth raises the maximum PRF; CW does not alias.' },
      { q: 'An angle of 20° between the beam and the flow underestimates velocity by approximately:', options: ['6%', '20%', '50%', 'Nothing — angles under 30° do not matter'], answer: 0, explain: 'cos 20° ≈ 0.94, so velocity is underestimated by ~6%. Keep the beam within 20° of the flow.' },
      { q: 'Peak aortic velocity is 4 m/s. The peak instantaneous gradient by the simplified Bernoulli equation is:', options: ['64 mmHg', '16 mmHg', '32 mmHg', '48 mmHg'], answer: 0, explain: 'ΔP = 4v² = 4 × 16 = 64 mmHg.' },
      { q: 'Axial resolution is improved by:', options: ['A higher transducer frequency', 'Increasing overall gain', 'Increasing imaging depth', 'Reducing the frame rate'], answer: 0, explain: 'Higher frequency shortens the spatial pulse length, improving axial resolution (at the cost of penetration).' },
      { q: 'Tissue harmonic imaging improves image quality mainly by:', options: ['Receiving at twice the transmitted frequency, reducing near-field clutter and side lobes', 'Transmitting at twice the usual frequency', 'Raising the pulse repetition frequency', 'Averaging several frames together'], answer: 0, explain: 'Harmonics are generated in tissue and are weaker from near-field reverberation and side lobes, so clutter falls.' },
      { q: 'Which change increases the 2D frame rate?', options: ['Narrowing the sector width', 'Increasing the depth', 'Adding colour Doppler', 'Increasing line density'], answer: 0, explain: 'Fewer scan lines (narrower sector or lower density) and a shallower depth take less time per frame.' },
      { q: 'Multiple equally spaced echoes deep to a strong reflector are:', options: ['Reverberation artefact', 'Acoustic shadowing', 'Mirror-image artefact', 'Side-lobe artefact'], answer: 0, explain: 'Sound bounces between the reflector and the transducer, producing repeats at multiples of the true depth.' },
      { q: 'Continuous-wave Doppler:', options: ['Measures high velocities without aliasing but has no range resolution', 'Has range resolution but aliases at high velocities', 'Is limited to velocities below 2 m/s', 'Is used to measure LVOT VTI'], answer: 0, explain: 'CW samples continuously along the whole beam, so velocities are unlimited but their depth is unknown. LVOT VTI uses PW.' },
      { q: 'On a standard colour Doppler map, red indicates flow:', options: ['Towards the transducer', 'Away from the transducer', 'Above the Nyquist limit', 'That is turbulent'], answer: 0, explain: 'BART: Blue Away, Red Towards. Turbulence shows as a mosaic (variance).' },
    ],
  },
  {
    id: 'views', title: 'Standard views & minimum dataset', icon: '🧭',
    objectives: [
      'Name each standard TTE view and the structures it shows',
      'Know which view is used for each key measurement',
      'Recognise views from the image alone',
    ],
    learn: [
      '<b>PLAX</b>: LV, LVOT, aortic root, mitral valve, LA. Linear LV measurements, LVOT diameter (zoomed, mid-systole), aortic root.',
      '<b>PSAX</b> levels: aortic valve (cusps, RVOT, pulmonary valve, TV), mitral valve (fish-mouth orifice, planimetry), papillary muscles (regional wall motion), apex.',
      '<b>A4C</b>: all four chambers, inferoseptal and anterolateral walls, mitral and tricuspid inflow, TAPSE (M-mode through the lateral tricuspid annulus), RV size (RVD1–3).',
      '<b>A2C</b>: anterior and inferior walls; with A4C for biplane Simpson\'s and LA volume. <b>A3C</b>: anteroseptal and inferolateral walls, LVOT, aortic valve.',
      '<b>A5C</b>: tilt anteriorly from A4C to open the LVOT and aortic valve — LVOT PW and AV CW Doppler.',
      '<b>Subcostal</b>: interatrial septum perpendicular to the beam (ASD), pericardial effusion, IVC (diameter ~1–2 cm from the RA junction, collapse with sniff).',
      '<b>Suprasternal</b>: aortic arch, coarctation, descending aorta flow reversal in aortic regurgitation.',
    ],
    tasks: [
      { id: 'tour', text: 'Step through PLAX, PSAX, PSAX AV, PSAX MV, A4C, A5C, A2C, A3C, Subcostal, SC IVC and SSN on the normal heart. Name the structures before reading the labels.', sim: { path: 'normal', view: 'PLAX' } },
      { id: 'labels', text: 'Turn off anatomy labels (L) and identify the chambers and valves in A4C yourself.', sim: { path: 'normal', view: 'A4C' } },
      { id: 'probe', text: 'Use the azimuth/tilt sliders to foreshorten A4C, then restore it. Note how the apex looks when the plane misses it.', sim: { path: 'normal', view: 'A4C' } },
    ],
    questions: [
      { q: 'Which view is used to measure the LVOT diameter for the continuity equation?', options: ['Zoomed parasternal long axis', 'Apical four-chamber', 'Parasternal short axis at papillary level', 'Subcostal'], answer: 0, explain: 'The LVOT is measured in zoomed PLAX in mid-systole, inner edge to inner edge, close to the aortic annulus.' },
      { q: 'Which view places the interatrial septum most perpendicular to the beam for assessing an ASD?', options: ['Subcostal', 'Apical four-chamber', 'Parasternal long axis', 'Suprasternal'], answer: 0, explain: 'In A4C the septum lies parallel to the beam and thins out (false dropout); subcostal images it side-on.' },
      { q: 'TAPSE is measured by M-mode through:', options: ['The lateral tricuspid annulus in A4C', 'The tricuspid leaflet tips in PSAX', 'The RV free wall in PLAX', 'The interventricular septum in A4C'], answer: 0, explain: 'TAPSE is the longitudinal excursion of the lateral tricuspid annulus in an RV-focused A4C.' },
      { q: 'Which apical view shows the anteroseptal and inferolateral walls?', options: ['Apical three-chamber (apical long axis)', 'Apical two-chamber', 'Apical four-chamber', 'Apical five-chamber'], answer: 0, explain: 'A3C: anteroseptal and inferolateral. A2C: anterior and inferior. A4C: inferoseptal and anterolateral.' },
      { q: 'The apical five-chamber view is obtained from the four-chamber view by:', options: ['Tilting the transducer anteriorly to open the LVOT and aortic valve', 'Rotating 90° anticlockwise', 'Sliding to the subcostal window', 'Tilting posteriorly towards the coronary sinus'], answer: 0, explain: 'Anterior tilt brings in the LVOT and aortic valve (posterior tilt shows the coronary sinus). Rotating anticlockwise gives A2C.' },
      { q: 'The suprasternal view is particularly used to assess:', options: ['The aortic arch, coarctation and descending aortic flow reversal', 'Mitral valve planimetry', 'Regional wall motion', 'TAPSE'], answer: 0, explain: 'SSN shows the arch and descending aorta: coarctation gradient and holodiastolic flow reversal in severe AR.' },
      { q: 'Name this standard view.', case: { path: 'normal', view: 'A2C', hideView: true }, options: ['Apical two-chamber', 'Apical four-chamber', 'Apical three-chamber', 'Parasternal long axis'], answer: 0, explain: 'Only the LV and LA are seen (no right heart), with the anterior and inferior walls: apical two-chamber.' },
      { q: 'Name this standard view.', case: { path: 'normal', view: 'PSAX_AV', hideView: true }, options: ['Parasternal short axis at aortic valve level', 'Parasternal short axis at mitral level', 'Apical five-chamber', 'Suprasternal'], answer: 0, explain: 'The three aortic cusps sit centrally ("Mercedes" sign), surrounded by the RA, tricuspid valve, RVOT and pulmonary valve.' },
      { q: 'Name this standard view.', case: { path: 'normal', view: 'SUBCOSTAL', hideView: true }, options: ['Subcostal four-chamber', 'Apical four-chamber', 'Parasternal long axis', 'Apical two-chamber'], answer: 0, explain: 'Liver in the near field and the heart lying across the screen with the septum perpendicular to the beam: subcostal.' },
      { q: 'Name this standard view.', case: { path: 'normal', view: 'PLAX', hideView: true }, options: ['Parasternal long axis', 'Apical three-chamber', 'Parasternal short axis', 'Subcostal'], answer: 0, explain: 'RV anteriorly, LV with the mitral valve and LVOT/aortic root, LA posteriorly, seen from the left parasternal window.' },
    ],
  },
  {
    id: 'lvsize', title: 'LV size, wall thickness & mass', icon: '📏',
    objectives: [
      'Measure LV dimensions correctly and know the BSE 2020 ranges',
      'Index volumes and mass to BSA, aortic dimensions to height',
      'Classify LV geometry from RWT and indexed mass',
    ],
    learn: [
      `Measure linear dimensions in PLAX, perpendicular to the LV long axis at the mitral leaflet tips: LVIDd at end-diastole, LVIDs at end-systole.`,
      `BSE 2020 LVIDd: men ≤ ${lvidM} mm, women ≤ ${lvidF} mm (graded mild / moderate / severe above). IVSd women ${ivsF[0]}–${ivsF[1]} mm.`,
      `Indexed LVEDV upper limit: men ${edviM} mL/m². Indexed LV mass normal: men ${lvmiM.lo}–${lvmiM.bands[0]} g/m², women ${GRADED.lvmi.female.lo}–${GRADED.lvmi.female.bands[0]} g/m².`,
      `<b>RWT = 2 × PWd / LVIDd</b>, cut-off ${RWT_CUTOFF}. Normal mass + RWT ≤ ${RWT_CUTOFF} = normal; normal mass + high RWT = concentric remodelling; high mass + high RWT = concentric hypertrophy; high mass + normal RWT = eccentric hypertrophy.`,
      'Index volumes and mass to body surface area; index aortic dimensions to height.',
      'BSE intervals do not apply to athletes, pregnancy or 3D measurements, and come from a largely Caucasian population.',
    ],
    tasks: [
      { id: 'calip', text: 'Freeze PLAX at end-diastole and measure LVIDd, IVSd and PWd with the calipers. Compare with the Echo report values.', sim: { path: 'normal', view: 'PLAX' } },
      { id: 'lvh', text: 'Load severe aortic stenosis. Measure IVSd and PWd, calculate RWT, and classify the geometry. Check against the report.', sim: { path: 'as', view: 'PLAX', grade: 'severe' } },
      { id: 'index', text: 'In the Echo report, change the patient from male to female and see how the same LV is regraded.', action: 'report', sim: { path: 'dcm', view: 'A4C' } },
    ],
    questions: [
      { q: 'BSE 2020 upper limit of normal for LVIDd in a man is:', options: [`${lvidM} mm`, '52 mm', '60 mm', '65 mm'], answer: 0, explain: `Men: ${GRADED.lvidd.male.lo}–${lvidM} mm; mild ${lvidM + 1}–${GRADED.lvidd.male.bands[1]}, moderate to ${GRADED.lvidd.male.bands[2]}, severe above.` },
      { q: 'BSE 2020 upper limit of normal for LVIDd in a woman is:', options: [`${lvidF} mm`, '56 mm', '45 mm', '59 mm'], answer: 0, explain: `Women: ${GRADED.lvidd.female.lo}–${lvidF} mm; mild ${lvidF + 1}–${GRADED.lvidd.female.bands[1]}, moderate to ${GRADED.lvidd.female.bands[2]}, severe above.` },
      { q: 'BSE 2020 normal range for IVSd in a woman is:', options: [`${ivsF[0]}–${ivsF[1]} mm`, '6–12 mm', '8–14 mm', '4–9 mm'], answer: 0, explain: `Women ${ivsF[0]}–${ivsF[1]} mm; men ${NORMAL.ivsd.male[0]}–${NORMAL.ivsd.male[1]} mm.` },
      { q: 'Relative wall thickness is calculated as:', options: ['2 × PWd / LVIDd', 'IVSd / PWd', '(IVSd + PWd) / LVIDs', 'PWd / LVIDs'], answer: 0, explain: `RWT = 2 × posterior wall thickness / LV internal diameter in diastole; > ${RWT_CUTOFF} is increased.` },
      { q: 'RWT 0.50 with a normal indexed LV mass indicates:', options: ['Concentric remodelling', 'Concentric hypertrophy', 'Eccentric hypertrophy', 'Normal geometry'], answer: 0, explain: 'Increased RWT without increased mass is concentric remodelling.' },
      { q: 'RWT 0.35 with an increased indexed LV mass indicates:', options: ['Eccentric hypertrophy', 'Concentric hypertrophy', 'Concentric remodelling', 'Normal geometry'], answer: 0, explain: 'Increased mass with normal RWT (dilated cavity) is eccentric hypertrophy, e.g. volume load from regurgitation.' },
      { q: 'Linear LV dimensions in PLAX should be measured:', options: ['Perpendicular to the LV long axis at the mitral leaflet tips', 'At the level of the papillary muscle heads in apical view', 'Parallel to the septum at the aortic annulus', 'At the LV apex'], answer: 0, explain: 'The ASE/BSE standard is perpendicular to the long axis at (or just below) the mitral leaflet tips.' },
      { q: 'BSE 2020 upper limit of normal for indexed LVEDV in a man is:', options: [`${edviM} mL/m²`, '58 mL/m²', '97 mL/m²', '120 mL/m²'], answer: 0, explain: `Men ${GRADED.lvedvi.male.lo}–${edviM} mL/m²; women ${GRADED.lvedvi.female.lo}–${GRADED.lvedvi.female.bands[0]} mL/m².` },
      { q: 'Indexed LV mass of 120 g/m² in a man is graded as:', options: ['Mildly increased', 'Normal', 'Moderately increased', 'Severely increased'], answer: 0, explain: `Men: normal ≤ ${lvmiM.bands[0]}, mild ${lvmiM.bands[0] + 1}–${lvmiM.bands[1]}, moderate ${lvmiM.bands[1] + 1}–${lvmiM.bands[2]}, severe > ${lvmiM.bands[2]} g/m².` },
      { q: 'BSE 2020 reference intervals should NOT be applied to:', options: ['Elite athletes', 'Women', 'Patients over 60', 'Patients with obesity'], answer: 0, explain: 'They are not valid for athletes, pregnancy or 3D measurements. They are sex-specific and indexed for body size.' },
    ],
  },
  {
    id: 'lvfunc', title: 'LV systolic function & wall motion', icon: '💓',
    objectives: [
      'Quote LVEF with the BSE bands',
      'Perform biplane Simpson\'s and know its pitfalls',
      'Map wall-motion abnormalities to coronary territories',
    ],
    learn: [
      `BSE LVEF bands: ≥ 55% normal, ${efBord}–54% borderline low, ${efImp}–49% impaired, ≤ 35% severely impaired. Quote the number; avoid "mild/moderate dysfunction".`,
      'Biplane Simpson\'s (method of discs) from A4C and A2C is the preferred 2D method. Foreshortened apical views underestimate volumes.',
      `GLS: normal is more negative than ${GLS_NORMAL}% (vendor-dependent).`,
      'Wall motion: normal, hypokinetic, akinetic (no thickening), dyskinetic (outward systolic motion), aneurysmal.',
      'Territories: <b>LAD</b> anterior wall, anterior septum, apex · <b>RCA</b> inferior wall, basal inferoseptum · <b>LCx</b> inferolateral and anterolateral walls.',
    ],
    tasks: [
      { id: 'simpson', text: 'On the dilated cardiomyopathy case, freeze and trace the endocardium at end-diastole and end-systole (Simpson tool) to calculate LVEF.', sim: { path: 'dcm', view: 'A4C' } },
      { id: 'rwma', text: 'Load the regional wall motion case. Step through systole and identify which segments fail to thicken; check the bullseye plot.', sim: { path: 'rwma', view: 'A4C' } },
    ],
    questions: [
      { q: 'In BSE terms, an LVEF of 42% is:', options: ['Impaired', 'Borderline low', 'Severely impaired', 'Normal'], answer: 0, explain: `BSE: ${efImp}–49% is impaired.` },
      { q: 'An LVEF of 52% should be reported as:', options: ['Borderline low', 'Normal', 'Mildly impaired', 'Impaired'], answer: 0, explain: `BSE: ${efBord}–54% is borderline low; "mildly impaired" wording is discouraged.` },
      { q: 'The preferred 2D method for calculating LVEF is:', options: ['Biplane Simpson\'s method of discs (A4C + A2C)', 'Teichholz from M-mode', 'Visual estimate only', 'Fractional shortening'], answer: 0, explain: 'Biplane disc summation does not assume a fixed LV shape, unlike Teichholz.' },
      { q: 'Foreshortening of the apical views typically causes:', options: ['Underestimation of LV volumes', 'Overestimation of LV volumes', 'No change in volumes', 'Overestimation of wall thickness only'], answer: 0, explain: 'Missing the true apex shortens the long axis, so disc volumes are too small.' },
      { q: 'In the 17-segment model, the LAD typically supplies:', options: ['The anterior wall, anterior septum and apex', 'The inferior wall and basal inferoseptum', 'The inferolateral wall only', 'The right ventricle only'], answer: 0, explain: 'LAD: anterior, anteroseptal and apical segments. Territories vary with coronary dominance.' },
      { q: 'Inferior wall akinesis is most often caused by occlusion of the:', options: ['Right coronary artery', 'Left anterior descending artery', 'First diagonal branch', 'Left main stem'], answer: 0, explain: 'The inferior wall is usually supplied by the RCA (or LCx if left-dominant).' },
      { q: 'Dyskinesis describes a segment that:', options: ['Moves outward in systole', 'Thickens normally', 'Does not move or thicken', 'Thickens less than normal'], answer: 0, explain: 'Hypokinesis = reduced, akinesis = absent, dyskinesis = paradoxical outward systolic motion.' },
      { q: 'This patient presents with breathlessness. What is the main finding?', case: { path: 'dcm', view: 'A4C' }, options: ['Dilated LV with globally reduced systolic function', 'Regional septal akinesis with normal LV size', 'Large pericardial effusion', 'Normal LV size and function'], answer: 0, explain: 'The LV is dilated and contracts poorly in every wall. The model measures LVEF ≈ 27% — severely impaired in BSE terms.' },
      { q: 'Describe the wall motion in this case.', case: { path: 'rwma', view: 'A4C' }, options: ['Septal and apical hypokinesis — LAD territory', 'Global hypokinesis of a dilated LV', 'Lateral wall akinesis — LCx territory', 'Normal wall motion'], answer: 0, explain: 'The septum and apex fail to thicken while the lateral wall contracts: a regional (LAD) abnormality. LVEF is borderline low (~51%).' },
    ],
  },
  {
    id: 'diastolic', title: 'Diastolic function & left atrium', icon: '🌊',
    objectives: [
      'Measure LA volume and know the BSE limits',
      'Interpret mitral inflow, e\' and E/e\'',
      'Apply the ASE/EACVI 2016 criteria',
    ],
    learn: [
      'LA volume: biplane disc summation or area–length from A4C and A2C at LV end-systole (just before mitral opening), avoiding the appendage and pulmonary veins.',
      'BSE 2020 LAVi: < 34 mL/m² normal, 34–38 borderline, > 38 dilated.',
      'Mitral inflow PW at the leaflet tips: E (early) and A (atrial) waves. Grade I (impaired relaxation) E/A ≤ 0.8; grade III (restrictive) E/A ≥ 2 with short deceleration time. No A wave in atrial fibrillation.',
      'Tissue Doppler e\': abnormal if septal < 7 cm/s or lateral < 10 cm/s. Average E/e\' > 14 suggests raised filling pressure.',
      'ASE/EACVI 2016 criteria with normal LVEF: average E/e\' > 14, septal e\' < 7 or lateral < 10 cm/s, TR velocity > 2.8 m/s, LAVi > 34 mL/m². More than half positive = diastolic dysfunction.',
    ],
    tasks: [
      { id: 'inflow', text: 'In A4C, read E and A on the mitral inflow trace in the normal heart, then in mitral stenosis.', sim: { path: 'normal', view: 'A4C' } },
      { id: 'la', text: 'Compare LA size and LAVi in the Echo report for the normal heart and severe mitral regurgitation.', action: 'report', sim: { path: 'mr', view: 'A4C', grade: 'severe' } },
    ],
    questions: [
      { q: 'BSE 2020: indexed LA volume is dilated above:', options: ['38 mL/m²', '28 mL/m²', '34 mL/m²', '48 mL/m²'], answer: 0, explain: '< 34 normal, 34–38 borderline, > 38 dilated.' },
      { q: 'LA volume should be measured:', options: ['At LV end-systole, just before mitral valve opening', 'At LV end-diastole', 'During atrial contraction', 'At peak E wave'], answer: 0, explain: 'The LA is largest at LV end-systole, before the mitral valve opens.' },
      { q: 'Which is NOT one of the four ASE/EACVI 2016 criteria for diastolic dysfunction with normal LVEF?', options: ['E/A ratio < 0.8', 'Average E/e\' > 14', 'TR velocity > 2.8 m/s', 'LAVi > 34 mL/m²'], answer: 0, explain: 'The four are average E/e\', septal or lateral e\', TR velocity and LAVi. E/A is used for grading, not detection.' },
      { q: 'Restrictive (grade III) filling is characterised by:', options: ['E/A ≥ 2 with a short deceleration time', 'E/A ≤ 0.8 with prolonged deceleration', 'Absent E wave', 'Normal E/A with a normal e\''], answer: 0, explain: 'High LA pressure drives a tall, rapidly decelerating E wave.' },
      { q: 'E/e\' is used to estimate:', options: ['LV filling pressure', 'LV ejection fraction', 'Pulmonary vascular resistance', 'Stroke volume'], answer: 0, explain: 'E reflects LA pressure and relaxation; dividing by e\' corrects for relaxation, giving a filling pressure estimate.' },
      { q: 'The mitral inflow PW sample volume should be placed:', options: ['At the mitral leaflet tips in A4C', 'At the mitral annulus in PLAX', 'In the LVOT in A5C', 'In the left atrium 2 cm from the valve'], answer: 0, explain: 'Leaflet tips give the peak inflow velocities; annular placement is used for TDI e\'.' },
      { q: 'The mitral A wave is absent in:', options: ['Atrial fibrillation', 'Sinus tachycardia', 'First-degree AV block', 'Aortic stenosis'], answer: 0, explain: 'No organised atrial contraction means no A wave.' },
      { q: 'In the ASE/EACVI 2016 criteria, septal e\' is abnormal below:', options: ['7 cm/s', '10 cm/s', '4 cm/s', '14 cm/s'], answer: 0, explain: 'Septal e\' < 7 cm/s or lateral e\' < 10 cm/s.' },
    ],
  },
  {
    id: 'aortic', title: 'Aortic valve', icon: '🔺',
    objectives: [
      'Grade aortic stenosis on concordant measures',
      'Use the continuity equation and know its pitfalls',
      'Recognise low-flow states and severe aortic regurgitation',
    ],
    learn: [
      'Severe AS: <b>Vmax ≥ 4 m/s, mean gradient ≥ 40 mmHg, AVA < 1.0 cm²</b>. Moderate: Vmax 3–3.9, mean 20–39, AVA 1.0–1.5. Mild: Vmax 2.5–2.9.',
      'Continuity: <b>AVA = (LVOT area × LVOT VTI) / AV VTI</b>. The LVOT diameter is squared, so it is the largest source of error.',
      'Dimensionless index (LVOT VTI / AV VTI) < 0.25 supports severe AS.',
      'Interrogate from apical, right parasternal and suprasternal windows; the highest velocity is the most parallel to the jet.',
      'Low flow: SVi < 35 mL/m². Low LVEF with AVA < 1.0 but low gradient → dobutamine stress echo to separate true-severe from pseudo-severe AS.',
      'Severe AR: vena contracta > 6 mm, pressure half-time < 200 ms, holodiastolic flow reversal in the descending aorta, RF ≥ 50%.',
    ],
    tasks: [
      { id: 'grade', text: 'Load aortic stenosis and switch between mild, moderate and severe. Watch Vmax, mean gradient and AVA change in the Measurements panel.', sim: { path: 'as', view: 'A5C', grade: 'moderate' } },
      { id: 'asreport', text: 'Open the Echo report on severe AS and write the conclusion yourself in the Write tab before reading the model\'s.', action: 'report', sim: { path: 'as', view: 'PLAX', grade: 'severe' } },
    ],
    questions: [
      { q: 'Severe aortic stenosis by peak velocity is:', options: ['≥ 4.0 m/s', '≥ 3.0 m/s', '≥ 2.5 m/s', '≥ 5.0 m/s'], answer: 0, explain: 'Vmax ≥ 4 m/s; 3.0–3.9 moderate; 2.5–2.9 mild.' },
      { q: 'Severe aortic stenosis by mean gradient is:', options: ['≥ 40 mmHg', '≥ 20 mmHg', '≥ 64 mmHg', '≥ 30 mmHg'], answer: 0, explain: 'Mean gradient ≥ 40 mmHg (64 mmHg is the peak gradient at 4 m/s).' },
      { q: 'The continuity equation calculates aortic valve area as:', options: ['LVOT area × LVOT VTI ÷ AV VTI', 'AV VTI ÷ LVOT VTI', 'LVOT diameter × AV Vmax', '220 ÷ pressure half-time'], answer: 0, explain: 'Flow through the LVOT equals flow through the valve: A₁V₁ = A₂V₂ using VTIs.' },
      { q: 'The largest source of error in the continuity equation is usually:', options: ['The LVOT diameter, because it is squared', 'The AV VTI', 'The LVOT VTI', 'Heart rate'], answer: 0, explain: 'Area = π(D/2)²: a 1 mm error in a 20 mm LVOT changes AVA by ~10%.' },
      { q: 'A dimensionless index (LVOT VTI / AV VTI) supporting severe AS is:', options: ['< 0.25', '< 0.5', '> 0.8', '< 0.9'], answer: 0, explain: 'DI < 0.25 supports severe AS and is independent of the LVOT diameter.' },
      { q: 'AVA 0.8 cm², mean gradient 25 mmHg, LVEF 30%. The most useful next test is:', options: ['Low-dose dobutamine stress echo', 'Exercise treadmill test', 'Repeat in two years', 'Transoesophageal echo for vegetations'], answer: 0, explain: 'Classical low-flow low-gradient AS: dobutamine shows whether AVA stays < 1.0 cm² as flow rises (true severe) and whether there is contractile reserve.' },
      { q: 'Low flow in aortic stenosis is defined as a stroke volume index below:', options: ['35 mL/m²', '50 mL/m²', '25 mL/m²', '65 mL/m²'], answer: 0, explain: 'SVi < 35 mL/m², even with normal LVEF (paradoxical low-flow).' },
      { q: 'Which finding indicates severe aortic regurgitation?', options: ['Holodiastolic flow reversal in the descending aorta', 'Pressure half-time of 600 ms', 'Vena contracta of 2 mm', 'Narrow central jet in the LVOT'], answer: 0, explain: 'Holodiastolic reversal, PHT < 200 ms and VC > 6 mm support severe AR.' },
      { q: 'Peak velocity 4.5 m/s, mean gradient 47 mmHg, AVA 0.74 cm² with normal LVEF. Grade the stenosis.', options: ['Severe', 'Moderate', 'Mild', 'Indeterminate'], answer: 0, explain: 'All three measures are concordant for severe AS.' },
      { q: 'Identify the main valve lesion in this case.', case: { path: 'as', view: 'PLAX', grade: 'severe' }, options: ['Aortic stenosis', 'Mitral stenosis', 'Mitral regurgitation', 'Normal valves'], answer: 0, explain: 'A thickened, restricted aortic valve with a high-velocity systolic jet and concentric LV hypertrophy.' },
    ],
  },
  {
    id: 'mitral', title: 'Mitral valve', icon: '🫀',
    objectives: [
      'Grade mitral regurgitation with an integrated approach',
      'Calculate EROA by PISA',
      'Assess mitral stenosis by planimetry and pressure half-time',
    ],
    learn: [
      'Severe MR (integrated): vena contracta ≥ 7 mm, EROA ≥ 0.40 cm², regurgitant volume ≥ 60 mL, regurgitant fraction ≥ 50%, systolic pulmonary-vein flow reversal.',
      'PISA: <b>EROA = 2πr² × Va / Vmax(MR)</b>; RVol = EROA × MR VTI.',
      'Primary MR: abnormal leaflets (prolapse, flail, rheumatic). Secondary (functional): LV or annular dilatation with tethering of normal leaflets.',
      'MS grading used here: MVA > 1.5 cm² mild, 1.0–1.5 moderate, < 1.0 severe. Report the mean gradient with the heart rate.',
      'Pressure half-time: <b>MVA = 220 / PHT</b>. Unreliable straight after balloon valvotomy, with significant AR or abnormal LV compliance.',
      'Rheumatic MS: commissural fusion, thickened tips, diastolic doming of the anterior leaflet ("hockey stick"), LA dilatation.',
    ],
    tasks: [
      { id: 'mrgrade', text: 'Load mitral regurgitation and change the grade. Compare the jet, LA size and regurgitant fraction in the report.', sim: { path: 'mr', view: 'A4C', grade: 'moderate' } },
      { id: 'msplan', text: 'Load severe mitral stenosis, go to PSAX MV and look at the fish-mouth orifice; then check the diastolic inflow on the trace.', sim: { path: 'ms', view: 'PSAX_MV', grade: 'severe' } },
    ],
    questions: [
      { q: 'A vena contracta indicating severe mitral regurgitation is:', options: ['≥ 7 mm', '≥ 3 mm', '≥ 5 mm', '≥ 12 mm'], answer: 0, explain: 'VC ≥ 7 mm severe, < 3 mm mild.' },
      { q: 'An EROA indicating severe primary MR is:', options: ['≥ 0.40 cm²', '≥ 0.10 cm²', '≥ 0.20 cm²', '≥ 1.0 cm²'], answer: 0, explain: 'EROA ≥ 0.40 cm² (40 mm²) for severe primary MR.' },
      { q: 'A regurgitant fraction indicating severe MR is:', options: ['≥ 50%', '≥ 30%', '≥ 20%', '≥ 70%'], answer: 0, explain: 'RF ≥ 50% severe, 30–49% moderate, < 30% mild.' },
      { q: 'PISA radius 1.0 cm, aliasing velocity 40 cm/s, peak MR velocity 5 m/s. EROA is approximately:', options: ['0.50 cm²', '0.25 cm²', '1.0 cm²', '0.08 cm²'], answer: 0, explain: 'EROA = 2π × 1² × 40 / 500 ≈ 0.50 cm².' },
      { q: 'Systolic flow reversal in the pulmonary veins suggests:', options: ['Severe mitral regurgitation', 'Mitral stenosis', 'Normal LA pressure', 'Aortic stenosis'], answer: 0, explain: 'A large regurgitant volume reverses pulmonary-vein flow in systole.' },
      { q: 'Secondary (functional) MR is caused by:', options: ['LV or annular dilatation with tethering of normal leaflets', 'Myxomatous leaflet prolapse', 'Rheumatic thickening', 'A flail posterior leaflet'], answer: 0, explain: 'In secondary MR the leaflets are normal; the ventricle and annulus are not.' },
      { q: 'Mitral valve area by pressure half-time is:', options: ['220 ÷ PHT', 'PHT ÷ 220', '4 × PHT²', '0.785 × PHT'], answer: 0, explain: 'MVA (cm²) = 220 / PHT (ms).' },
      { q: 'Pressure half-time is unreliable for MVA:', options: ['Immediately after balloon mitral valvotomy', 'In sinus rhythm', 'With a mean gradient over 5 mmHg', 'In the A4C view'], answer: 0, explain: 'Abrupt changes in LA/LV compliance after valvotomy, significant AR or stiff LVs make PHT inaccurate.' },
      { q: 'What does the colour jet in this case show?', case: { path: 'mr', view: 'A4C', grade: 'severe' }, options: ['A systolic jet into the LA — mitral regurgitation', 'Accelerated diastolic inflow — mitral stenosis', 'A systolic jet into the RA — tricuspid regurgitation', 'A shunt across the ventricular septum'], answer: 0, explain: 'A mosaic systolic jet fills the left atrium from the mitral valve. The model regurgitant fraction is ~54% — severe.' },
      { q: 'Identify the lesion in this case.', case: { path: 'ms', view: 'A4C', grade: 'severe' }, options: ['Mitral stenosis', 'Mitral regurgitation', 'Aortic stenosis', 'Normal mitral valve'], answer: 0, explain: 'Thickened, restricted mitral leaflets with accelerated, prolonged diastolic inflow and a dilated LA.' },
    ],
  },
  {
    id: 'rightheart', title: 'Right heart & pulmonary hypertension', icon: '🫁',
    objectives: [
      'Measure RV size and function and know the BSE limits',
      'Estimate RA pressure and PASP',
      'Assign the echocardiographic probability of PH',
    ],
    learn: [
      `RV function: TAPSE ≥ ${tapse} mm, RV s\' ≥ ${rvs} cm/s, FAC ≥ ${facM}% (men) / ${LIMIT.rvfac.female.min}% (women). RVD1 (base) men ≤ ${rvd1M} mm.`,
      '<b>PASP = 4 × TR Vmax² + RAP</b>.',
      'RAP from the IVC: ≤ 21 mm with > 50% sniff collapse ≈ 3 mmHg; > 21 mm with < 50% collapse ≈ 15 mmHg; otherwise ≈ 8 mmHg.',
      'PH probability (ESC/ERS 2022): TR Vmax ≤ 2.8 m/s low, 2.9–3.4 intermediate, > 3.4 high, raised by supporting signs (RV/LV > 1, septal flattening, RVOT AT < 105 ms, PA > 25 mm, dilated IVC, RA area > 18 cm²).',
      'Septal flattening in systole and diastole (D-shaped LV) = RV pressure overload; diastole only = volume overload.',
    ],
    tasks: [
      { id: 'tapse', text: 'In A4C, note TAPSE and RV size in the Echo report for the normal heart, then for pulmonary hypertension.', action: 'report', sim: { path: 'phtn', view: 'A4C' } },
      { id: 'dsign', text: 'Load pulmonary hypertension in PSAX and watch the septum flatten into a D-shape through the cycle.', sim: { path: 'phtn', view: 'PSAX' } },
      { id: 'ivc', text: 'Open the subcostal IVC view and estimate RA pressure.', sim: { path: 'tr', view: 'SC_IVC', grade: 'severe' } },
    ],
    questions: [
      { q: 'TAPSE is abnormal below:', options: [`${tapse} mm`, '25 mm', '10 mm', '21 mm'], answer: 0, explain: `TAPSE < ${tapse} mm indicates reduced RV longitudinal function.` },
      { q: 'Tissue Doppler RV s\' is abnormal below:', options: [`${rvs} cm/s`, '15 cm/s', '5 cm/s', '12 cm/s'], answer: 0, explain: `RV s' < ${rvs} cm/s indicates RV systolic dysfunction.` },
      { q: 'BSE upper limit of normal for RVD1 (basal diameter) in a man is:', options: [`${rvd1M} mm`, '35 mm', '55 mm', '41 mm'], answer: 0, explain: `RVD1 men ${NORMAL.rvd1.male[0]}–${rvd1M} mm, women ${NORMAL.rvd1.female[0]}–${NORMAL.rvd1.female[1]} mm.` },
      { q: 'Estimated pulmonary artery systolic pressure is:', options: ['4 × TR Vmax² + RA pressure', '4 × TR Vmax²', 'TR Vmax × 4 + 10', 'Mean TR gradient + 5'], answer: 0, explain: 'The RV–RA gradient (4v²) plus RA pressure gives RV, and so PA, systolic pressure (without RVOT obstruction).' },
      { q: 'TR Vmax 3.0 m/s with no other signs of PH gives an echocardiographic probability of PH that is:', options: ['Intermediate', 'Low', 'High', 'Diagnostic of PH'], answer: 0, explain: '2.9–3.4 m/s without other signs is intermediate. Echo gives a probability, not a diagnosis.' },
      { q: 'IVC 24 mm with less than 50% collapse on sniff. Estimated RA pressure:', options: ['15 mmHg', '3 mmHg', '8 mmHg', '25 mmHg'], answer: 0, explain: 'Dilated, non-collapsing IVC ≈ 15 mmHg.' },
      { q: 'IVC 18 mm with more than 50% collapse on sniff. Estimated RA pressure:', options: ['3 mmHg', '8 mmHg', '15 mmHg', '20 mmHg'], answer: 0, explain: 'Normal-size, collapsing IVC ≈ 3 mmHg.' },
      { q: 'A D-shaped LV in PSAX in both systole and diastole indicates:', options: ['RV pressure overload', 'RV volume overload only', 'LV hypertrophy', 'A normal finding'], answer: 0, explain: 'Pressure overload flattens the septum throughout the cycle; volume overload mainly in diastole.' },
      { q: 'What does this short-axis view show?', case: { path: 'phtn', view: 'PSAX' }, options: ['A dilated RV with a D-shaped LV — RV pressure overload', 'Normal biventricular size', 'Concentric LV hypertrophy', 'Pericardial effusion'], answer: 0, explain: 'The RV is dilated and the septum flattened. The model: TR Vmax 4.3 m/s, PASP ~89 mmHg — high probability of PH.' },
      { q: 'Identify the main finding in this case.', case: { path: 'tr', view: 'A4C', grade: 'severe' }, options: ['Tricuspid regurgitation', 'Mitral regurgitation', 'Atrial septal defect', 'Normal right heart'], answer: 0, explain: 'A systolic jet from the tricuspid valve into a dilated RA.' },
    ],
  },
  {
    id: 'pericongenital', title: 'Pericardium & shunts', icon: '🔄',
    objectives: [
      'Size a pericardial effusion and recognise tamponade physiology',
      'Distinguish pericardial from pleural fluid',
      'Know the types of ASD and VSD and how to assess a shunt',
    ],
    learn: [
      'Effusion size (maximal end-diastolic depth): < 1 cm small, 1–2 cm moderate, > 2 cm large.',
      'Tamponade signs: RA systolic collapse, RV diastolic collapse, plethoric IVC, respiratory variation in mitral E > 25% / tricuspid E > 40%. It is a clinical diagnosis.',
      'In PLAX, pericardial fluid lies <b>anterior</b> to the descending aorta; pleural fluid extends posterior to it.',
      'ASD: secundum (commonest, fossa ovalis), primum (part of AVSD, cleft mitral leaflet), sinus venosus (often with anomalous right pulmonary venous drainage), coronary sinus.',
      'VSD: perimembranous (commonest in adults), muscular, outlet, inlet. A restrictive VSD has a high-velocity jet (large LV–RV gradient, normal RV pressure).',
      'Qp:Qs = (RVOT area × RVOT VTI) / (LVOT area × LVOT VTI). > 1.5 suggests a significant shunt.',
    ],
    tasks: [
      { id: 'eff', text: 'Load pericardial effusion in PLAX; measure the echo-free depth at end-diastole and relate it to the descending aorta.', sim: { path: 'effusion', view: 'PLAX' } },
      { id: 'asd', text: 'Load the atrial septal defect in the subcostal view and find the shunt with colour; then look at RV size in A4C.', sim: { path: 'asd', view: 'SUBCOSTAL' } },
      { id: 'vsd', text: 'Load the VSD and catch the systolic jet crossing the septum.', sim: { path: 'vsd', view: 'A4C' } },
    ],
    questions: [
      { q: 'A pericardial effusion with a maximal end-diastolic depth of 2.5 cm is:', options: ['Large', 'Moderate', 'Small', 'Trivial'], answer: 0, explain: '> 2 cm is large; 1–2 cm moderate; < 1 cm small.' },
      { q: 'Which echo finding supports tamponade physiology?', options: ['RV diastolic collapse', 'A small IVC with > 50% collapse', 'LV hypertrophy', 'Mild aortic regurgitation'], answer: 0, explain: 'RV diastolic and RA systolic collapse support tamponade; a collapsing IVC argues against it.' },
      { q: 'Respiratory variation suggesting tamponade is an inspiratory fall in mitral E velocity of more than:', options: ['25%', '5%', '10%', '60%'], answer: 0, explain: 'Mitral E > 25% (tricuspid > 40%) respiratory variation.' },
      { q: 'In PLAX, pericardial (rather than pleural) fluid is recognised because it:', options: ['Lies anterior to the descending aorta', 'Lies posterior to the descending aorta', 'Is always echogenic', 'Is only seen in systole'], answer: 0, explain: 'Pericardial fluid tracks between the heart and the descending aorta; pleural fluid extends behind the aorta.' },
      { q: 'The most common type of atrial septal defect is:', options: ['Secundum', 'Primum', 'Sinus venosus', 'Coronary sinus'], answer: 0, explain: 'Secundum defects at the fossa ovalis account for about three-quarters.' },
      { q: 'A sinus venosus ASD is classically associated with:', options: ['Anomalous right upper pulmonary venous drainage', 'A cleft anterior mitral leaflet', 'Bicuspid aortic valve', 'Coarctation of the aorta'], answer: 0, explain: 'Superior sinus venosus defects usually come with partial anomalous pulmonary venous drainage.' },
      { q: 'A primum ASD is associated with:', options: ['A cleft anterior mitral leaflet', 'Anomalous pulmonary venous drainage', 'Ebstein anomaly', 'Pulmonary stenosis'], answer: 0, explain: 'Primum ASD is a partial atrioventricular septal defect, typically with a cleft mitral leaflet.' },
      { q: 'A VSD jet peak velocity of 5 m/s indicates:', options: ['A restrictive defect with normal RV systolic pressure', 'A large non-restrictive defect', 'Eisenmenger physiology', 'Severe pulmonary hypertension'], answer: 0, explain: '4 × 5² = 100 mmHg LV–RV gradient: RV pressure is low, so the defect is restrictive.' },
      { q: 'What is the main finding in this case?', case: { path: 'effusion', view: 'PLAX' }, options: ['Pericardial effusion', 'Pleural effusion only', 'Dilated cardiomyopathy', 'Normal study'], answer: 0, explain: 'An echo-free space surrounds the heart, deepest posteriorly — the model measures 2.3 cm, a large effusion.' },
      { q: 'What is the main finding in this case?', case: { path: 'asd', view: 'SUBCOSTAL' }, options: ['Atrial septal defect with left-to-right flow', 'Ventricular septal defect', 'Tricuspid regurgitation', 'Normal interatrial septum'], answer: 0, explain: 'Low-velocity flow crosses the atrial septum from LA to RA, best seen subcostally.' },
      { q: 'What is the main finding in this case?', case: { path: 'vsd', view: 'A4C' }, options: ['Ventricular septal defect', 'Atrial septal defect', 'Mitral regurgitation', 'Normal study'], answer: 0, explain: 'A high-velocity systolic jet crosses the ventricular septum from LV to RV.' },
    ],
  },
  {
    id: 'reporting', title: 'Reporting', icon: '📝',
    objectives: [
      'Structure a report and use BSE wording',
      'Write a conclusion that answers the clinical question',
      'Handle discordant measurements',
    ],
    learn: [
      'Report structure by structure: LV → RV → atria → valves → aorta → pericardium → IVC/RAP.',
      'Quote LVEF as a number with the BSE band; grade, then give the number ("moderately dilated LV, LVEDVi 88 mL/m²").',
      'Echo gives the <b>probability</b> of pulmonary hypertension, not a diagnosis.',
      'Conclusion: abnormal findings first, most important first, each with its number; compare with previous studies.',
      'When measures disagree (e.g. AVA < 1.0 cm² but Vmax < 4 m/s), check the LVOT measurement and SVi and state the discordance.',
      'The full guide and templates are in the Echo report dialog and docs/REPORTING.md.',
    ],
    tasks: [
      { id: 'write', text: 'Pick any pathology, open the Echo report, write your own report from a template, then compare with the model conclusion.', action: 'report', sim: { path: 'mr', view: 'A4C', grade: 'severe' } },
      { id: 'female', text: 'Report the dilated cardiomyopathy case for a 60-year-old woman, 160 cm, 60 kg, and note which gradings change.', action: 'report', sim: { path: 'dcm', view: 'A4C' } },
    ],
    questions: [
      { q: 'The BSE advises that LV systolic function is reported as:', options: ['LVEF as a number with the BSE band (e.g. "LVEF 42%, impaired")', '"Mild", "moderate" or "severe" LV dysfunction', 'Fractional shortening only', 'A visual description without numbers'], answer: 0, explain: 'Quote the LVEF and the BSE band; avoid mild/moderate dysfunction wording.' },
      { q: 'An echo study can establish:', options: ['The probability of pulmonary hypertension', 'A definitive diagnosis of pulmonary hypertension', 'Pulmonary vascular resistance precisely', 'The cause of PH in all cases'], answer: 0, explain: 'PH is diagnosed at right-heart catheterisation; echo assigns low, intermediate or high probability.' },
      { q: 'A good report conclusion should:', options: ['List abnormal findings in order of importance, each with its supporting number', 'Repeat every measurement from the report body', 'Only state whether the study is normal', 'Avoid numbers'], answer: 0, explain: 'The conclusion answers the clinical question with the key abnormal findings and their numbers.' },
      { q: 'Aortic root dimensions are indexed to:', options: ['Height', 'Body surface area', 'Weight', 'Age'], answer: 0, explain: 'BSE 2020 indexes aortic dimensions to height; volumes and mass to BSA.' },
      { q: 'AVA 0.9 cm², Vmax 3.5 m/s, normal LVEF. Before reporting, you should:', options: ['Check the LVOT diameter and SVi and report the discordance', 'Report severe AS based on AVA alone', 'Report moderate AS based on Vmax alone', 'Ignore the AVA'], answer: 0, explain: 'Discordant grading needs a check for LVOT under-measurement or low flow, and the discordance stated.' },
      { q: 'BSE 2020: an indexed LA volume of 36 mL/m² is reported as:', options: ['Borderline', 'Normal', 'Dilated', 'Severely dilated'], answer: 0, explain: '34–38 mL/m² is borderline.' },
      { q: 'An LVEF of 33% is reported as:', options: ['Severely impaired', 'Impaired', 'Moderately impaired', 'Borderline low'], answer: 0, explain: '≤ 35% is severely impaired in BSE terms.' },
      { q: 'When comparing with a previous study, the report should:', options: ['State whether there is a significant change and whether it could be technical', 'Copy the previous conclusion', 'Not mention the previous study', 'Average the two sets of measurements'], answer: 0, explain: 'Say what has changed, whether it is outside measurement variability, and the likely reason.' },
    ],
  },
];

export const MODULE_BY_ID = Object.fromEntries(MODULES.map((m) => [m.id, m]));

// every question gets a stable id "<module>-<index>"
MODULES.forEach((m) => m.questions.forEach((q, i) => { q.id = `${m.id}-${i}`; q.module = m.id; }));

// ---- scoring and quiz building (pure, testable) ---------------------------------
export function shuffle(arr, rand = Math.random) {
  const a = arr.slice();
  for (let i = a.length - 1; i > 0; i--) { const j = Math.floor(rand() * (i + 1)); [a[i], a[j]] = [a[j], a[i]]; }
  return a;
}
// a question prepared for display: options in random order, `correct` is the
// index of the right answer in that order
export function prepare(q, rand = Math.random) {
  const order = shuffle(q.options.map((_, i) => i), rand);
  return { ...q, shown: order.map((i) => q.options[i]), correct: order.indexOf(q.answer) };
}
export function moduleQuiz(id, rand = Math.random) {
  return shuffle(MODULE_BY_ID[id].questions, rand).map((q) => prepare(q, rand));
}
// mock exam: an even spread across modules, then shuffled
export function mockExam(n = EXAM.questions, rand = Math.random) {
  const pools = MODULES.map((m) => shuffle(m.questions, rand));
  const picked = [];
  for (let round = 0; picked.length < n; round++) {
    let added = false;
    for (const p of pools) if (round < p.length && picked.length < n) { picked.push(p[round]); added = true; }
    if (!added) break;
  }
  return shuffle(picked, rand).map((q) => prepare(q, rand));
}
// answers: array of chosen display indices (null = unanswered)
export function score(quiz, answers) {
  const byModule = {};
  let correct = 0;
  quiz.forEach((q, i) => {
    const ok = answers[i] === q.correct;
    if (ok) correct++;
    const b = byModule[q.module] || (byModule[q.module] = { correct: 0, total: 0 });
    b.total++; if (ok) b.correct++;
  });
  const pct = quiz.length ? correct / quiz.length : 0;
  return { correct, total: quiz.length, pct, passed: pct >= PASS_MARK, byModule };
}
// modules scoring below the pass mark in an attempt, weakest first
export function weakModules(result) {
  return Object.entries(result.byModule)
    .map(([id, b]) => ({ id, pct: b.correct / b.total }))
    .filter((x) => x.pct < PASS_MARK)
    .sort((a, b) => a.pct - b.pct)
    .map((x) => x.id);
}
