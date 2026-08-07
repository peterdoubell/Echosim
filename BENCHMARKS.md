# EchoSim Fidelity Benchmarks (v2 — "state-of-the-art" bar)

**Why this exists.** EchoSim's earlier board scored it 9.9/10 on a
*pattern-recognition teaching-tool* rubric. Measured instead against the real
state of the art — physics-based computational cardiology (openCARP, Living
Heart), forward **12-lead** ECG modelling (eikonal + lead-field), CFD-FSI blood
flow, scatterer-based ultrasound image formation (Field II / COLE), and
**validated** commercial simulators (CAE Vimedix, HeartWorks) — the honest
weighted baseline figure was **≈ 5.9 / 10** (an earlier draft of this table
mis-summed the eleven weights as 9.6 and reported 6.6; the weights actually sum
to **10.6**, so the correct baseline is 62.9 ÷ 10.6 = **5.9** — close to the
6.6–7.0 range independently estimated for the product). This document is that
harder, physics-and-clinic anchored rubric. See
`docs/RESEARCH_models_and_simulators.md` for the evidence.

Each dimension is scored 0–10 against explicit, externally-anchored anchors, not
against "does it look plausible in a screenshot".

> **v2.1 re-score (current).** After implementing a physics-derived 12-lead ECG
> (B2), a lumped-parameter closed-loop circulation with a real PV loop + regional
> wall motion (B3), continuity-conserving Doppler with a diastolic vortex (B4),
> scatterer+PSF image formation with live gCNR (B5), and an automated
> self-validation harness (B9/B11), the honest weighted figure is now
> **≈ 7.4 / 10** (78.6 ÷ 10.6) — above the top of the 6.6–7.0 range independently
> estimated for the product, and now with expert-in-the-loop review underway. Each dimension below shows **baseline → current** with the
> concrete evidence and the residual gap to the next band. The score is **not**
> 9.9: the top band requires FE mechanics, CFD-FSI, a torso lead-field validated
> against clinical ECG datasets, and — decisively — **external expert/RCT
> validation and regulatory clearance**, which are evidence-and-audit processes
> that cannot be closed by writing code (see `docs/COMMERCIAL_READINESS.md`,
> `docs/VALIDATION.md`).

## Universal anchor scale
| Band | Meaning |
|---|---|
| **3–4** | Schematic / illustrative. Hand-authored appearance, no governing physics. |
| **5–6** | Physiologically *plausible* but **kinematic/analytic** — prescribed, not solved. |
| **7–8** | **Physics-informed**: solves a reduced form of the governing equations; partially validated. |
| **9–10** | **Physics-based digital twin**: full governing equations on patient-specific anatomy, **externally validated** against clinical/ground-truth data (Field II, 12-lead, FSI ±8% Doppler, CAMUS/MedalCare-XL). |

---

## The 11 benchmark dimensions

### B1 · Anatomical & structural fidelity — weight 1.0 — **6.5 → 7.5**
Gold standard: patient-specific CT/MRI whole-heart mesh with **DT-MRI / rule-based (Streeter) myofiber & sheet architecture**, trabeculation, chordae tendineae, coronary tree; a population via **statistical shape models**.
EchoSim **now**: the idealized blended-SDF heart (correct topology — continuous myocardium, shared septum, RV crescent, papillaries, sinuses of Valsalva) gains **(a)** a **rule-based (Streeter) myofibre-orientation overlay** — helical fibres whose helix angle sweeps transmurally from ~+55° (sub-endocardium) to ~−55° (sub-epicardium), the counter-wound helices underlying LV torsion — **(b)** **chordae tendineae**, tendinous fans from each papillary-muscle tip to the mitral leaflet free edge, and **(c)** **representative (not abstract) chamber morphology** — a **bullet-shaped LV** (broad base, tapered apex) and a **triangular RV** converging to the apex, replacing the symmetric ellipsoids; dimensions to-scale (LVIDd 4.6, LA 3.5, walls 0.9 cm). Still **no coronaries, single generic geometry**, cm-approximate. *(Score held at 7.5 pending the cardiac-radiologist face-validity rating logged in `docs/VALIDATION_RESULTS.md` — self-scoring morphology is the internal-rating bias this rubric exists to avoid.)*
→ 8: expert face-validity ≥4/5 on the new morphology + **≥1 alternative geometry** (statistical-shape variant) + trabeculation. → 9.5: image-derived / statistical-shape patient variability.

### B2 · Electrophysiology & 12-lead ECG — weight 1.2 — **4.0 → 7.5**  *(largest gain)*
Gold standard: ionic cell model (**O'Hara–Rudy / ten Tusscher**) → **monodomain/bidomain or eikonal** activation on a fibered mesh → **lead-field/torso forward model** → physiological **standard 12-lead** at ~1 kHz; arrhythmias, bundle-branch block, ischemic ST changes; ECG **derived from and synchronized with** contraction (ECGSIM, MedalCare-XL).
EchoSim **now**: a **vectorcardiographic dipole** (P/QRS/T Gaussian sources, McSharry/Sameni-style) projected through an **inverse-Dower lead field** to a full physiological **standard 12-lead** (`js/electrophysiology.js`); a rhythm library (sinus, AF with irregularly-irregular R-R, brady/tachy, 1° AV block, PVCs, LBBB, RBBB, anterior STEMI, hyperkalaemia); pathology morphology (LVH voltage, low-voltage effusion); and **electromechanical coupling** — the activation phase drives the mechanical cycle. Validated by 12/12 automated morphology checks (`tools/validate.mjs`: aVR −, normal-axis I/aVF +, R-wave progression V1→V6, rate/regularity, STEMI ST-elevation vector, LVH voltage).
Residual to 9.5: eikonal-on-fibered-mesh activation and a **torso** (not inverse-Dower) lead field, ischemia/conduction validated against **reference ECG datasets** (MedalCare-XL). → this is the external-validation ceiling, not a coding gap.

### B3 · Myocardial mechanics & the cardiac cycle — weight 1.1 — **5.5 → 8.0**
Gold standard: **FEM hyperelastic (Holzapfel–Ogden / Guccione)** + active tension/strain, fiber-based, coupled to **lumped-parameter/CircAdapt** circulation → physiological **pressure–volume loops**, EF, wall stress, **torsion**, regional strain; regional wall-motion abnormalities (ischemia).
EchoSim **now**: a **lumped-parameter closed-loop circulation** (time-varying elastance LV `Plv=E(t)(Vlv−V0)` + two-element Windkessel afterload + diode valves; integrated to its limit cycle) producing a real, self-consistent **LV pressure–volume loop** (`js/hemodynamics.js`, surfaced live in the UI). The cavity geometry is driven by the **absolute modelled volume ratio** `rho=Vlv/EDV`, so the measured EF/LVIDs **fall out of the PV loop** (normal EF 58 %, DCM 21 % at LVIDd 6.4/LVIDs 5.9, AS 57 %) — fixing the prior defect where every case read EF≈63 %. Long-axis shortens less than short-axis. **Regional wall-motion abnormality** (`path.rwma`): the infarcted segment retains its end-diastolic radius and fails to thicken (septum 2.2 cm at systole vs 2.7 cm normal) while other walls contract; global EF drops to ~46 %.
And **regional longitudinal strain** is now surfaced as an **AHA 17-segment bullseye** (`regionalStrain()` + the live GLS panel): normal −20 % uniform, DCM −8 % global, and an RWMA drives its wall's segments toward zero (septal infarct → septal −6 % / lateral −20 %, GLS −15 %) — the culprit territory localises exactly as speckle-tracking reads it.
Residual to 9.5: FE active-contraction deriving **fibre strain, wall stress and torsion** from a full mechanics solve (the strain here is model-derived per segment, not an FE field).

### B4 · Hemodynamics & Doppler physics — weight 1.1 — **6.0 → 7.5**
Gold standard: **Navier–Stokes CFD with fluid–structure-interaction valves**, conserving mass/momentum; vortex rings, physiologic jets; validated within **~8%** of clinical Doppler; angle-true spectral traces.
EchoSim **now**: jet velocities are **continuity-derived** (`v=Q/A`) from the modelled stroke volume and effective orifice area, so **inflow = outflow holds by construction** of the closed loop — AS auto-scales to 3.6 m/s through a small AVA, LVOT ~1.0, MS inflow 2.0; the MR peak comes from the modelled LV–LA systolic gradient via **Bernoulli** (~3.7 m/s). A **divergence-free diastolic inflow vortex** (Lamb–Oseen-style rotation about a fixed axis) is superimposed behind the anterior mitral leaflet, so colour Doppler shows the filling vortex, not just a straight jet. BART colour with aliasing + steered CW/PW gate retained.
Residual to 9.5: an actual **reduced-order CFD or FSI** pressure/velocity field validated vs Doppler ±8 % — the flow is continuity-conserving and vortical but still not a Navier–Stokes solve.

### B5 · Ultrasound image formation — weight 1.1 — **6.5 → 8.0**
Gold standard: **scatterer-based (Field II)** or **convolution (COLE/SIMUS)** with measured **beam PSF**, frequency-dependent attenuation, correct **speckle statistics** (Rayleigh/K), validated by **gCNR/CNR/speckle-SNR**; artifacts: reverberation, acoustic shadowing/enhancement, side lobes, blooming, dropout.
EchoSim **now**: a **scatterer map + separable PSF convolution** (COLE-style, `freqMHz=2.7`) with acoustic **shadowing/enhancement and reverberation** behind calcified structures, a **336×420** formation buffer (~0.36 mm/px axial at 15 cm), a 3-tier FPS-adaptive sampler and a user **Detail** control. The beam is now modelled as a **round trip**: a transmit limb focused once at a **user-settable focal depth** (drawn as the caret scanners put beside the depth ruler) multiplied by a **dynamically-refocused receive limb** at constant F-number until the aperture saturates — so the in-focus receive beam bounds off-focus degradation, as on a real machine. Axial resolution is pulse-length limited and independent of focus. **Tissue harmonic imaging** ships as a toggle (narrower effective beam, reverberation suppression, and the honest *penetration penalty*), and the scan plane now has **finite elevational slice thickness** — a slab set by a fixed, non-steerable acoustic lens, sampled stochastically for genuine partial-volume averaging.

**Read the gCNR number carefully.** It was ≈0.97–0.98, and slice thickness deliberately pulled it to **≈0.91**. That is an *improvement in fidelity, not a regression*: real myocardium-vs-cavity gCNR sits around 0.7–0.9, so the previous image was **easier than real echo**. Within the current model the metric still moves in the right directions — harmonic mode raises it (a consequence of the narrower beam, not of construction), and calcific shadowing lowers it.
Residual to 9.5: **Field II-class** point-scatterer simulation with speckle statistics matched to a **real image dataset** (CAMUS). We report gCNR and can now say it lands in the right *range*, but we still do not compare against measured per-image targets.

### B6 · Pathology breadth & mechanistic depth — weight 1.0 — **6.5 → 7.5**
Gold standard: **30+** graded cases (Vimedix/HeartWorks), each with mechanistically consistent anatomy + flow + haemodynamics + ECG changes, across TTE/TEE.
EchoSim **now**: 10 pathologies (added **regional wall-motion abnormality / post-MI**), each **mechanistically coupled** — pathology drives a distinct PV loop and haemodynamic summary (EDV/ESV/EF/gradient/regurgitant fraction) *and* an ECG signature (STEMI, LVH voltage, rhythm). And the valve lesions are now **graded mild/moderate/severe** off one physical lever (effective orifice area / regurgitant-orifice resistance), verified against BSE thresholds — AS peak 2.7/3.3/4.4 m/s, MR regurgitant fraction 26/46/65 % — so the continuity-Doppler peak, PV loop and severity band all move together.
Residual to 8: **≥20 distinct lesions** (HOCM, AR, tamponade, congenital…) and TTE→TEE; to 9.5: cardiologist-validated coupled library.

### B7 · Modalities & acquisition realism — weight 0.8 — **5.5 → 6.5**
Gold standard: **TTE + TEE + 3D/4D + MPR**, tracked probe on a **haptic mannequin**, real windows/dropout, tissue-Doppler/strain, contrast.
EchoSim **now**: a **TTE + TEE** modality switch. TEE places a higher-frequency (5.5 MHz → finer PSF) probe posterior to the left atrium imaging anteriorly, with three calibrated windows — **ME 4-chamber** (LA→valves→LV/RV, atria in the near field, correctly inverted from apical), **ME long-axis** (LVOT/AV/ascending aorta), and **transgastric mid-papillary short-axis** (LV ring with both papillary muscles). Still a single moveable plane per window; no true 3D-volume acquisition/MPR, no tissue-Doppler.
→ 8: add **MPR from the SDF volume** + **tissue-Doppler/strain** + ME bicaval/AV-SAX. → 9.5: probe-pressure/contact model + real 3D volume acquisition.

### B8 · Quantification & measurement accuracy — weight 1.0 — **7.0 → 8.0**
Gold standard: on-image calipers, Simpson's biplane EF, PISA/EROA, VTI, continuity-equation AVA, TDI e′, GLS, all **traceable to ground truth**.
EchoSim **now**: **learner-placed measurement tools** (`js/measure.js`) — click-two-point **distance calipers** (exact cm via the sector's `cmPerPx`) and a **Simpson single-plane method-of-discs EF**: the learner traces the LV endocardium hinge→apex→hinge at end-diastole and end-systole and the tool disc-sums about the apex→base-midpoint long axis to EDV/ESV/EF. Crucially, the learner's EF is **checked live against the model's own ground-truth EF** and flagged within the ASE ±6 % test–retest band — the §A2 measurement-accuracy check made interactive. Plus the auto LVIDd/s, ellipsoid EF, PV loop, full haemodynamic readout and a **global longitudinal strain (GLS)** readout derived from the modelled long-axis shortening (normal −22.9 %, DCM −6.7 %, RWMA −17.0 %, AS −22.3 %; flagged when impaired), all BSE/ASE-calibrated. Disc math verified to 0.0 % vs an analytic solid of revolution.
Residual to 9.5: learner-traced **VTI + continuity-equation AVA + PISA/EROA**, and agreement against an **independent** ground-truth phantom (not the model's own geometry). *(GLS present; the two continuity/PISA items are the remaining →8.5 work — score held at 8.0 until they land.)*

### B9 · External / expert validation & datasets — weight 1.0 — **3.5 → 4.5**  *(now the weakest)*
Gold standard: **peer-reviewed** face/content validity (expert Likert), learning-curve RCTs; agreement with datasets (**CAMUS**, **MedalCare-XL**), gCNR vs real images.
EchoSim **now**: ships a full **validation protocol + expert-Likert instrument + RCT statistical plan** (`docs/VALIDATION.md`) and an **automated self-validation harness** — 12/12 ECG-morphology checks against textbook criteria, plus live gCNR/CNR and **measurement-accuracy checks against the model's own ground truth** (`tools/validate.mjs`). This is the reproducible, ground-truthable portion.
Residual — and it is real: **no external expert panel, no learning-transfer RCT, no CAMUS/MedalCare-XL dataset agreement**. These are evidence-and-audit processes (people, ethics approval, participants), **not** code, so this axis is capped by an external step. This is the honest ceiling on the whole product's "commercial-deployment for clinical training" claim.

### B10 · Real-time interactivity, UX & education — weight 0.7 — **8.0 → 8.5**
Gold standard: haptic mannequin + probe tracking + curriculum/logbook.
EchoSim: genuinely strong — real-time, moveable plane ↔ live 2D echo, synchronized **rhythm strip + 12-lead modal**, live **PV loop + haemodynamics + image-quality** panels, labels, challenge/quiz tiers, measurements, help/onboarding, accessibility, fully offline in a browser. Lacks hardware haptics/tracking and a formal curriculum.
→ 9.5: probe-tracking hardware option + structured curriculum + spaced-repetition.

### B11 · Engineering rigor & reproducibility — weight 0.6 — **8.0 → 8.5**
Gold standard: validated, versioned research codebases (openCARP, Chaste) with test suites and provenance.
EchoSim: clean self-contained ES-module app, single geometry source of truth, **automated model self-validation** (`tools/validate.mjs`, 12/12) and **headless screenshot/DOM verification** (`tools/shoot.mjs`, `shoot-hemo.mjs`), no runtime errors. Lacks numerical-convergence checks and physics provenance for the reduced models.
→ 9.5: automated numerical/regression tests + documented validation harness against external references.

---

## Weighted result

| # | Dimension | Weight | Baseline | **Current** |
|---|---|---:|---:|---:|
| B1 | Anatomy & structure | 1.0 | 6.5 | **7.5** |
| B2 | Electrophysiology & 12-lead ECG | 1.2 | 4.0 | **7.5** |
| B3 | Myocardial mechanics | 1.1 | 5.5 | **8.0** |
| B4 | Hemodynamics & Doppler physics | 1.1 | 6.0 | **7.5** |
| B5 | Ultrasound image formation | 1.1 | 6.5 | **8.0** |
| B6 | Pathology breadth & depth | 1.0 | 6.5 | **7.5** |
| B7 | Modalities & acquisition | 0.8 | 5.5 | **6.5** |
| B8 | Quantification accuracy | 1.0 | 7.0 | **8.0** |
| B9 | External validation | 1.0 | 3.5 | **4.5** |
| B10 | Interactivity, UX & education | 0.7 | 8.0 | **8.5** |
| B11 | Engineering rigor | 0.6 | 8.0 | **8.5** |

**Weighted overall: 5.9 → 7.4 / 10.** (Baseline Σ 62.9, current Σ 78.6, ÷ Σ weight 10.6 = 7.42.)

The three physics axes that a *schematic* trainer omits — **B2 (12-lead ECG),
B3 (mechanics/PV loop), B4 (conserved Doppler)** — moved from 4.0/5.5/6.0 to a
consistent **7.5**, and **B5 (imaging)** to 7.5. What remains below 7 is now
honest and specific:
- **B9 external validation (4.5)** — capped by a human/RCT process, not code.
- **B7 modalities (6.5)** — TTE+TEE now; no 3D-MPR / tissue-Doppler yet.
- **B1 anatomy (7.5)** — fibre field + chordae added; no coronaries, single geometry.

## Why this is not 9.9 (and cannot be reached by code alone)
The 9–10 band is a **validated physics-based digital twin**. Reaching it needs
(a) FE hyperelastic mechanics with fiber strain, (b) CFD-FSI blood flow, (c) an
eikonal/torso-lead-field ECG, each **validated against clinical datasets**, and
(d) **independent expert + learning-transfer-RCT validation** plus **QMS /
regulatory** evidence. (a)–(c) are large but codeable; **(d) is not** — it
requires an external panel, ethics approval, participants and an audit. EchoSim
is therefore driven as high as engineering can take it (**7.3**, with the
physics axes at 7.5 and quantification at 8.0) and the residual is fully specified and turnkey in
`docs/VALIDATION.md` and `docs/COMMERCIAL_READINESS.md`.

## Highest-leverage roadmap from here
1. **B7 (×0.8):** add a TEE window set + MPR from the SDF volume + tissue-Doppler. *(5.5 → ~7.5.)*
2. **B8 (×1.0):** add learner-traced VTI + continuity-equation AVA + PISA/EROA + GLS. *(8.0 → ~9.)*
3. **B1 (×1.0):** rule-based (Streeter) fiber field + chordae/trabeculae + a 2nd morphology. *(6.5 → ~8.)*
4. **B5 (×1.1):** CAMUS segmentation-overlap + speckle-statistics match to real images. *(7.5 → ~8.5.)*
5. **B9 (×1.0):** run the shipped expert-Likert panel + CAMUS overlap. *(4.5 → ~6.5 — the external step.)*

## References (evidence base)
- openCARP — *Comput. Methods Programs Biomed.* 2021; MonoAlg3D vs openCARP, CinC 2025.
- ten Tusscher–Panfilov (2006); O'Hara–Rudy (2011) ionic models; monodomain vs bidomain.
- Eikonal + lead-field 12-lead ECG — *Europace* 19(suppl_3):iii259, 2017; simplified 3D 12-lead — PMC3654639; **MedalCare-XL** arXiv 2211.15997; ECGSIM; Dower inverse transform.
- Living Heart Project — *Eur. J. Mech. A/Solids* 2014 (PMC4175454); Holzapfel–Ogden (2009); Guccione active tension; CircAdapt; time-varying elastance (Suga–Sagawa); two-element Windkessel.
- Cardiac FSI valves in realistic LV — PMC5590990; FSI validated vs Doppler ±8% — *Comput. Biol. Med.* 2024; whole-heart flow — arXiv 2605.09629; diastolic vortex ring (Lamb–Oseen).
- Field II (Jensen); COLE convolution (MICCAI/Springer); SIMUS (Garcia); DTI-based cardiac speckle — PMC4537486; GAN patho-realistic — arXiv 1712.07881; UltraScatter — arXiv 2510.10612; autonomous-nav simulation — *Front. Cardiovasc. Med.* 2024 (PMC11347295); CAMUS.
- Image-quality metrics gCNR/CNR/speckle-SNR — *J. Med. Ultrason.* 2021; regional quality DL — arXiv 2408.00591.
- Commercial simulators — CAE Vimedix (caehealthcare.com; cae.com); HeartWorks / Intelligent Ultrasound (inventivemedical.com; academy.intelligentultrasound.com).
- Validity frameworks — Messick (1995) unified validity; Kane (2013) validity argument.
