# Echo reporting guide and templates

EchoSim builds a structured transthoracic echo report for every simulated case.
To open it, click **📝 Echo report & templates** in the Measurements panel. The
dialog has three tabs:

- **Report** shows the full measurement table with each value graded against
  the reference interval, structure-by-structure findings, and a numbered
  conclusion. Everything is computed from the simulated heart.
- **Write your report** gives you a template pre-filled with this case's
  numbers. The judgements are left in `[square brackets]` for you to make from
  the images. Compare your conclusion with the Report tab when you have
  finished.
- **Reporting guide** is a short version of this document.

Enter the patient's sex, age, height and weight first. Reference intervals are
sex-specific, chamber volumes and LV mass are indexed to body surface area
(Mosteller), and aortic dimensions are indexed to height. The patient details
are saved in the browser. The report is not available in challenge mode,
because it would give the answer away.

Teaching use only. The patient is simulated, and the numbers come from a
lumped-parameter model, not from a real study.

---

## 1. How to write the report

1. **Patient and study details.** Record the indication, rhythm, heart rate,
   image quality, height, weight and BSA.
2. **Report in a fixed order, one structure at a time:**
   LV (size, wall thickness, geometry, systolic function, regional wall motion,
   diastolic function) → RV → atria → valves (aortic, mitral, tricuspid,
   pulmonary) → aorta → pericardium → IVC and RA pressure.
3. **Measure on the right frame.** End-diastole is the frame just before the
   mitral valve closes, or the largest LV cavity. End-systole is the frame just
   before the mitral valve opens, or the smallest cavity, after the aortic
   valve has closed. Measure linear LV dimensions in PLAX, perpendicular to the
   long axis, at the level of the mitral leaflet tips.
4. **Quote LVEF as a number** ("LVEF 42%"). The BSE asks reporters to avoid
   "mild" and "moderate" LV dysfunction and to use the LVEF bands in §2.2.
5. **Grade, then give the number:** "moderately dilated LV (LVEDVi 88 mL/m²)",
   not just "dilated".
6. **Grade valve disease on concordant measures.**
   - Aortic stenosis: Vmax, mean gradient and AVA.
   - Mitral stenosis: valve area and mean gradient, with the heart rate.
   - Regurgitation: an integrated assessment of regurgitant fraction and
     volume, vena contracta, jet and flow reversal.
   - When the measures disagree, say so and give the likely reason, such as
     low flow or a poor Doppler angle.
7. **Make the conclusion answer the clinical question.** List the abnormal
   findings first, most important first, each with the number that supports
   it. State the PH probability and the RA pressure estimate, and compare with
   any previous study.

### Words to use

| Say | Avoid |
|---|---|
| LVEF 42%, impaired | "mild–moderate LV dysfunction" |
| Moderately dilated LV (LVEDVi 88 mL/m²) | "big LV" |
| Severe aortic stenosis (Vmax 4.4 m/s, mean 46 mmHg, AVA 0.74 cm²) | "tight AS" |
| Intermediate echocardiographic probability of PH | "PH" (echo cannot diagnose PH, only estimate its probability) |
| No pericardial effusion | "no PE" (ambiguous) |

---

## 2. Reference intervals (BSE 2020)

Source: Harkness A, Ring L, Augustine DX, Oxborough D, Robinson S, Sharma V.
*Normal reference intervals for cardiac dimensions and function for use in
echocardiographic practice: a guideline from the British Society of
Echocardiography.* Echo Res Pract 2020;7(1):G1–G18. Right heart: Zaidi A et al.
Echo Res Pract 2020;7(1):G19–G41. The values below are transcribed from the BSE
poster *Echocardiography: Reference Intervals & Functional Assessment* (March
2020).

The intervals come from a largely Caucasian population. They do not apply to
3D measurements, athletes or pregnancy.

The app uses the same numbers from `js/references.js`. If you edit either
place, keep the two in step.

### 2.1 LV size and mass (graded)

A value at or below the upper limit of normal is normal. Above it, the value is
graded mild, moderate or severe. A value below the lower limit is "below
normal".

| Measure | Sex | Normal | Mild | Moderate | Severe |
|---|---|---|---|---|---|
| LVIDd (mm) | ♂ | 37–56 | 57–61 | 62–65 | > 65 |
| | ♀ | 35–51 | 52–55 | 56–59 | > 59 |
| LVIDs (mm) | ♂ | 22–41 | 42–45 | 46–50 | > 50 |
| | ♀ | 20–37 | 38–42 | 43–46 | > 46 |
| LVEDV/BSA (mL/m²) | ♂ | 30–79 | 80–91 | 92–103 | > 103 |
| | ♀ | 29–70 | 71–81 | 82–91 | > 91 |
| LVESV/BSA (mL/m²) | ♂ | 9–31 | 32–36 | 37–42 | > 42 |
| | ♀ | 8–27 | 28–32 | 33–37 | > 37 |
| LV mass/BSA (g/m²) | ♂ | 40–110 | 111–127 | 128–145 | > 145 |
| | ♀ | 33–99 | 100–115 | 116–131 | > 131 |

Normal, not graded:

| Measure | ♂ | ♀ |
|---|---|---|
| IVSd (mm) | 6–12 | 5–11 |
| PWd (mm) | 6–12 | 6–12 |
| LVEDV (mL) | 53–156 | 46–121 |
| LVESV (mL) | 15–62 | 13–47 |
| LV mass (g) | 72–219 | 51–173 |

**LV geometry.** Relative wall thickness is RWT = 2 × PWd / LVIDd, with a
cut-off of 0.42.

| | Indexed mass normal | Indexed mass increased |
|---|---|---|
| RWT ≤ 0.42 | Normal geometry | Eccentric hypertrophy |
| RWT > 0.42 | Concentric remodelling | Concentric hypertrophy |

### 2.2 LV systolic and diastolic function

| LVEF (%) | BSE wording |
|---|---|
| ≥ 55 | Normal |
| 50–54 | Borderline low |
| 36–49 | Impaired |
| ≤ 35 | Severely impaired |

- GLS: −18% or more negative is normal.
- MAPSE and mean mitral annular s′ lower limits (cm/s): < 40 y 6.4; 40–60 y
  5.7; > 60 y 4.9.
- Diastolic function: E/A 0.8–2.0. Interpret it with e′, E/e′, LAVi and TR
  Vmax, using the BSE/ASE algorithm.

### 2.3 Atria

| Measure | Normal | Borderline | Dilated |
|---|---|---|---|
| LA volume / BSA (mL/m²) | < 34 | 34–38 | > 38 |

| Measure | ♂ | ♀ |
|---|---|---|
| RA area (cm²) | ≤ 22 | ≤ 19 |
| RA area / BSA (cm²/m²) | ≤ 11 | ≤ 11 |

### 2.4 Right ventricle

| Measure | ♂ | ♀ |
|---|---|---|
| RVD1, base (mm) | 26–47 | 22–43 |
| RVD2, mid (mm) | 19–42 | 17–35 |
| RVD3, length (mm) | 55–87 | 51–80 |
| RVOT proximal, PLAX (mm) | 25–43 | 22–40 |
| RVED area / BSA (cm²/m²) | ≤ 13.6 | ≤ 12.6 |
| RV FAC (%) | ≥ 30 | ≥ 35 |
| TAPSE (mm) | ≥ 17 | ≥ 17 |
| RV s′ (cm/s) | ≥ 9 | ≥ 9 |

### 2.5 Aorta

Measure at end-diastole, inner edge to inner edge, and index to height.

| Measure (mm/m) | ♂ | ♀ |
|---|---|---|
| Sinus of Valsalva / height | 13.8–21.8 | 13.1–20.7 |
| Sino-tubular junction / height | 11.4–18.6 | 11.0–17.8 |
| Proximal ascending aorta / height | 11.5–19.9 | 11.4–19.8 |

---

## 3. Valve and haemodynamic grading

These are the thresholds the BSE adopts from the ASE/EACVI valve guidelines.

**Aortic stenosis** (Baumgartner 2017)

| | Sclerosis | Mild | Moderate | Severe |
|---|---|---|---|---|
| Vmax (m/s) | < 2.5 | 2.5–2.9 | 3.0–3.9 | ≥ 4.0 |
| Mean gradient (mmHg) | | < 20 | 20–39 | ≥ 40 |
| AVA (cm²) | | > 1.5 | 1.0–1.5 | < 1.0 |

If AVA is below 1.0 cm² but Vmax is below 4 m/s, check for low flow
(SVi < 35 mL/m²) before you call the stenosis moderate.

**Mitral stenosis** (MVA by planimetry or pressure half-time): mild above
1.5 cm², moderate 1.0–1.5 cm², severe below 1.0 cm². Report the mean gradient
with the heart rate.

**Regurgitation** (integrated assessment, Zoghbi 2017)

| | Mild | Moderate | Severe |
|---|---|---|---|
| Regurgitant fraction (%) | < 30 | 30–49 | ≥ 50 |
| Regurgitant volume (mL) | < 30 | 30–59 | ≥ 60 |

**Pulmonary hypertension probability** (ESC/ERS 2022). TR Vmax ≤ 2.8 m/s gives
low probability, 2.9–3.4 m/s intermediate, and > 3.4 m/s high. Supporting
signs move the probability up:

- RV/LV basal ratio > 1
- Septal flattening
- RVOT acceleration time < 105 ms
- PA diameter > 25 mm
- IVC > 21 mm with reduced collapse
- RA area > 18 cm²

Estimate PASP as 4 × TR Vmax² plus RA pressure.

**RA pressure from the IVC** (Rudski 2010)

| IVC | Collapse with sniff | RAP |
|---|---|---|
| ≤ 21 mm | > 50% | ~3 mmHg |
| ≤ 21 mm | < 50% | ~8 mmHg |
| > 21 mm | > 50% | ~8 mmHg |
| > 21 mm | < 50% | ~15 mmHg |

**Pericardial effusion** (maximal echo-free depth at end-diastole): small
< 1 cm, moderate 1–2 cm, large > 2 cm. Tamponade is a clinical and
haemodynamic diagnosis. Look for:

- RA systolic collapse or RV diastolic collapse
- A plethoric IVC
- Respiratory variation in mitral inflow > 25% or tricuspid inflow > 40%

---

## 4. How EchoSim measures

The Report tab takes its numbers from the simulated anatomy and circulation in
the same way a sonographer would. `js/report.js` → `measureModel()` does the
work:

| Measurement | Method in the simulator |
|---|---|
| LVIDd, IVSd, PWd, LVIDs | Calipers along a PLAX line at the mitral leaflet tips, at end-diastole and end-systole |
| LV mass | ASE cube formula: 0.8 × 1.04 × [(IVSd + LVIDd + PWd)³ − LVIDd³] + 0.6 g |
| LVEDV, LVESV, LVEF | Cavity volume from the circulation model, which matches Simpson's biplane |
| GLS, MAPSE, E/A | Myocardial strain and annular motion over the cycle; mitral inflow E and A peaks |
| LA volume | Voxel count of the LA cavity at LV end-systole |
| RA area, RVED area, RV FAC | Areas in the apical 4-chamber plane |
| RVD1 | Basal RV diameter in A4C |
| TAPSE | Lateral tricuspid annular excursion between RV end-diastole and end-systole |
| Aortic dimensions | Inner-edge calipers at the sinus, STJ and proximal ascending aorta |
| AV Vmax, mean gradient | Peak and ejection-averaged LV–aortic pressure difference (simplified Bernoulli) |
| AVA | Continuity equation: forward stroke volume / AV VTI |
| MVA, mitral mean gradient | Effective orifice and diastolic LA–LV gradient |
| MR RF and RVol | Regurgitant volume divided by total LV stroke volume |
| TR Vmax, PASP | RV–RA pressure difference, plus RAP estimated from the IVC |
| Effusion | Echo-free depth behind the LV posterior wall in PLAX at end-diastole |

Linear measurements are rounded to whole millimetres before grading, as they
would be on a report.

---

## 5. Templates

Each template is in the **Write your report** tab. Placeholders such as
`{LVEF}` are filled with the case's numbers. Text in `[brackets]` is for you to
decide from the images. The template text is in `js/report.js` (`TEMPLATES`).

### Blank structured report

```
Indication: [ ]
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
1. [ ]
```

### Normal study

```
Left ventricle: normal size (LVIDd {LVIDD} mm, LVEDVi {LVEDVI} mL/m²) and wall thickness (IVSd {IVSD} mm, PWd {PWD} mm). LVEF {LVEF}% (Simpson's biplane). No regional wall motion abnormality. Normal diastolic function (E/A {EA}).
Right ventricle: normal size and systolic function (TAPSE {TAPSE} mm, FAC {RVFAC}%).
Atria: normal LA (LAVi {LAVI} mL/m²) and RA size.
Valves: structurally and functionally normal. No significant regurgitation or stenosis.
Aorta: normal root and proximal ascending aorta (sinus {SOV} mm).
Pericardium: no effusion. IVC {IVC} mm with > 50% collapse (RAP ~3 mmHg).

Conclusion:
1. Structurally normal heart with normal biventricular size and systolic function (LVEF {LVEF}%).
```

### Aortic stenosis

```
Aortic valve: [calcified / restricted cusps, number of cusps]. Vmax {AVVMAX} m/s, mean gradient {AVMEAN} mmHg, AVA {AVA} cm² (continuity, LVOT [ ] mm).
Left ventricle: wall thickness IVSd {IVSD} mm, PWd {PWD} mm — [geometry]. LVEF {LVEF}%.
Check: [low-flow state? SVi = SV / BSA]  [aortic regurgitation]  [ascending aorta {ASC} mm].

Conclusion:
1. [Grade] aortic stenosis (Vmax {AVVMAX} m/s, mean gradient {AVMEAN} mmHg, AVA {AVA} cm²).
2. LV [geometry], LVEF {LVEF}%.
```

### Mitral regurgitation

```
Mitral valve: [leaflet morphology: prolapse / flail / restriction, which scallop]. Jet [central / eccentric, direction].
Quantification: regurgitant fraction {MRRF}%, regurgitant volume {MRRVOL} mL. [Vena contracta, PISA EROA, pulmonary-vein flow].
Left ventricle: LVEDVi {LVEDVI} mL/m², LVESVi {LVESVI} mL/m², LVEF {LVEF}%.
Left atrium: LAVi {LAVI} mL/m².

Conclusion:
1. [Grade] [primary / secondary] mitral regurgitation (RF {MRRF}%, RVol {MRRVOL} mL).
2. LV [size], LVEF {LVEF}%; LA [size].
```

### Mitral stenosis

```
Mitral valve: [rheumatic thickening, commissural fusion, doming anterior leaflet, subvalvular involvement].
MVA {MVA} cm² ([planimetry / pressure half-time]); mean gradient {MVMEAN} mmHg at HR [ ].
Left atrium: LAVi {LAVI} mL/m². [LAA thrombus?]
Pulmonary pressure: TR Vmax [ ] m/s, PASP [ ] mmHg.

Conclusion:
1. [Grade] rheumatic mitral stenosis (MVA {MVA} cm², mean gradient {MVMEAN} mmHg).
2. [LA dilatation]  [pulmonary hypertension probability].
```

### Dilated cardiomyopathy

```
Left ventricle: [dilated], LVIDd {LVIDD} mm, LVEDVi {LVEDVI} mL/m². Global hypokinesis. LVEF {LVEF}% (Simpson's biplane), GLS {GLS}%.
Wall thickness IVSd {IVSD} mm, PWd {PWD} mm. [Thrombus at apex?]
Mitral valve: [functional MR from annular dilatation / tenting], RF {MRRF}%.
Right ventricle: TAPSE {TAPSE} mm, FAC {RVFAC}%.
Atria: LAVi {LAVI} mL/m².

Conclusion:
1. [Dilated] LV with [severely impaired] systolic function, LVEF {LVEF}%.
2. [Grade] functional mitral regurgitation.  3. [RV function]
```

### Regional wall motion abnormality

```
Left ventricle: normal size. Regional wall motion: [segments, AHA 17-segment model, hypokinetic / akinetic / dyskinetic] — [coronary territory].
LVEF {LVEF}% (Simpson's biplane). GLS {GLS}% [pattern on bullseye].
[Complications: thrombus, aneurysm, VSD, MR, effusion]

Conclusion:
1. Regional wall motion abnormality in the [ ] territory, LVEF {LVEF}%.
```

### Pericardial effusion

```
Pericardium: [global / localised] pericardial effusion, maximal diastolic depth {EFFUSION} cm [location].
Haemodynamics: [RA / RV diastolic collapse]  [respiratory variation in mitral / tricuspid inflow]  IVC {IVC} mm [plethoric?].
[Fibrin strands / echogenic content]

Conclusion:
1. [Size] pericardial effusion ({EFFUSION} cm), [with / without] echocardiographic features of tamponade.
```

### Ventricular septal defect

```
Ventricular septum: [location: perimembranous / muscular] defect, [size] mm, left-to-right flow, peak velocity [ ] m/s (LV–RV gradient [ ] mmHg).
Left ventricle: LVEDVi {LVEDVI} mL/m² [volume load?]. LVEF {LVEF}%.
Right ventricle / pulmonary pressure: [ ].

Conclusion:
1. [Restrictive] [type] VSD with left-to-right shunt.
```

### Atrial septal defect

```
Interatrial septum: [secundum / primum / sinus venosus] defect, [size] mm, [rims], left-to-right flow on colour Doppler.
Right heart: RVD1 {RVD1} mm, RVED area {RVEDAREA} cm², RA area {RAAREA} cm² [volume load?]. TAPSE {TAPSE} mm.
Pulmonary pressure: TR Vmax [ ] m/s.

Conclusion:
1. [Type] atrial septal defect with left-to-right shunt.  2. [RV dilatation]
```

### Tricuspid regurgitation

```
Tricuspid valve: [leaflet morphology, annular diameter]. TR [grade] (vena contracta [ ] mm, hepatic-vein flow [ ]).
TR Vmax {TRVMAX} m/s; estimated PASP {PASP} mmHg (RAP from IVC {IVC} mm).
Right heart: RVD1 {RVD1} mm, RA area {RAAREA} cm², TAPSE {TAPSE} mm.

Conclusion:
1. [Grade] tricuspid regurgitation. TR Vmax {TRVMAX} m/s — [probability] of pulmonary hypertension.
```

### Pulmonary hypertension

```
Right ventricle: [dilated], RVD1 {RVD1} mm, RVED area {RVEDAREA} cm². [Septal flattening / D-shaped LV in PSAX]. TAPSE {TAPSE} mm, FAC {RVFAC}%.
TR Vmax {TRVMAX} m/s; estimated PASP {PASP} mmHg (RAP {RAP} mmHg).
Additional PH signs: [RV/LV basal ratio > 1, PA diameter > 25 mm, RVOT acceleration time < 105 ms, IVC > 21 mm].
RA area {RAAREA} cm².

Conclusion:
1. [High] echocardiographic probability of pulmonary hypertension (TR Vmax {TRVMAX} m/s, PASP ~{PASP} mmHg).
2. [RV dilatation and systolic function].
```

---

## 6. Worked examples: the simulator's conclusions

These are the conclusions the Report tab generates for each built-in case
(default patient: male, 45 y, 175 cm, 75 kg, BSA 1.91 m²). Use them to check
your own write-up.

| Case | Model conclusion |
|---|---|
| Normal | Structurally normal heart. Normal biventricular size and systolic function. No significant valve disease. |
| Aortic stenosis | 1. Severe aortic stenosis (Vmax 4.5 m/s, mean gradient 47 mmHg, AVA 0.74 cm²). 2. LV concentric hypertrophy. |
| Mitral regurgitation | 1. Severe mitral regurgitation (RF 54%, RVol 66 mL). 2. Moderately dilated LV (LVEDVi 102 mL/m²). 3. Dilated left atrium (LAVi 67 mL/m²). |
| Mitral stenosis | 1. Severe mitral stenosis (MVA 1.00 cm², mean gradient 13.6 mmHg). 2. Dilated left atrium (LAVi 89 mL/m²). 3. Intermediate echocardiographic probability of pulmonary hypertension (TR Vmax 3.4 m/s, PASP ~54 mmHg). |
| Dilated cardiomyopathy | 1. Severely impaired LV systolic function, LVEF 27%. 2. Severely dilated LV (LVEDVi 129 mL/m²). 3. Moderate functional mitral regurgitation (RF 31%). 4. Dilated left atrium (LAVi 70 mL/m²). 5. Impaired RV systolic function (TAPSE 12 mm, FAC 30%). |
| Regional wall motion abnormality | 1. Borderline low LV systolic function, LVEF 51%. 2. Regional wall motion abnormality: septal and apical hypokinesis (LAD territory). |
| Pericardial effusion | 1. Large pericardial effusion (2.3 cm). |
| Tricuspid regurgitation | 1. Severe tricuspid regurgitation (TR Vmax 2.5 m/s). 2. RV dilatation (RVD1 57 mm). |
| Pulmonary hypertension | 1. LV concentric remodelling. 2. High echocardiographic probability of pulmonary hypertension (TR Vmax 4.3 m/s, PASP ~89 mmHg). 3. RV dilatation (RVD1 55 mm) with impaired RV systolic function (TAPSE 14 mm, FAC 30%). |
| Atrial septal defect | 1. Secundum ASD. 2. RV dilatation (RVD1 54 mm). |
| Ventricular septal defect | 1. Restrictive muscular VSD. |

A few reporting points from these cases:

- **Mitral stenosis.** The mean gradient depends on heart rate, so always
  report it with the HR. Severity rests on the valve area.
- **Pulmonary hypertension.** The LV is small and under-filled, and the septum
  is flattened. "Concentric remodelling" there is a geometry label from RWT,
  not LV disease, so comment on septal shape in the PSAX view.
- **Tricuspid regurgitation.** In severe TR the RA v-wave blunts the RV–RA
  gradient: the CW envelope is dense, peaks early and is triangular, and its
  Vmax (2.5 m/s here) can underestimate PASP. Report the IVC and RAP estimate.
- **Mitral stenosis, right heart.** Look for the secondary pulmonary hypertension: a TR jet
  (Vmax 3.4 m/s), a mildly dilated RV and a dilated IVC. Measure the pressure
  half-time on the CW "ski slope" (MVA ≈ 220/PHT).
- **Dilated cardiomyopathy.** With LAP ~18 mmHg the inflow is restrictive
  (E/A > 2, DT < 150 ms), not impaired relaxation.

`node tools/verify-report.mjs` checks that every case's conclusion names its
lesion, that the BSE grading and indexing are correct, and that every template
fills.
