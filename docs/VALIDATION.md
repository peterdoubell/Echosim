# EchoSim — Validation Protocol & Self-Validation Harness (B9)

Validity for a simulator is an **argument** (Messick/Kane), built from several
evidence sources. This document specifies (A) the **automated self-validation**
we run in-repo/in-app, (B) the **expert-review instrument**, and (C) the
**learning-transfer study** — the last two are the external steps and are fully
protocolised here so the study is turnkey.

---

## A. Automated self-validation (shipped, quantitative)

### A1. Image quality — gCNR / CNR
Real ultrasound realism is quantified with the **generalized contrast-to-noise
ratio (gCNR)**, which is robust to dynamic-range changes. We compute, live, the
gCNR and CNR between **myocardium (target)** and **LV-cavity blood (background)**
from the rendered B-mode and surface them in the app's Validation panel.
- **Acceptance target:** gCNR in the range reported for real cardiac B-mode
  tissue-vs-lumen contrast (typically ~0.7–0.95 for well-separated targets).
- **Method:** histogram-overlap form, gCNR = 1 − ∫ min(p_target, p_background).

### A2. Measurement accuracy vs ground truth
Because EchoSim owns its geometry, we can check that what the **image/caliper
pathway reports** matches the **model's own ground truth** (the true chamber
dimensions/volumes/velocities in `cardiac-model.js`). The harness reports the
error for: LVIDd, LVIDs, Simpson's-disc EF, LVOT/AV peak velocity and gradient,
MR/AS severity band.
- **Acceptance target:** |measured − ground truth| within clinical repeatability
  (e.g. LVIDd ≤ 2 mm, EF ≤ 5 %, peak V ≤ 0.2 m/s), i.e. the tool does not add
  error beyond its own discretisation.

### A3. ECG morphology checks (automated, 12-lead)
The 12-lead engine is checked against textbook criteria for each rhythm:
| Check | Expectation |
|---|---|
| Axis (sinus) | Normal: I + / aVF + ; aVR always − |
| R-wave progression | V1→V6 R amplitude increases, transition ~V3–V4 |
| Intervals | PR 120–200 ms, QRS ≤ 110 ms (normal), QT rate-appropriate (Bazett) |
| AF | No organised P, irregularly-irregular R-R, fibrillatory baseline |
| LBBB/RBBB | QRS > 120 ms with the correct terminal-force axis (rSR′ V1 for RBBB) |
| STEMI | ST-segment elevation vector projecting to the expected lead group |
| LVH | Voltage criteria increased; strain T-wave changes |
`tools/validate.mjs` runs A1–A3 headless and prints a pass/fail table. This is
the reproducible, ground-truthable portion of validity.

---

## B. Expert face/content-validity instrument (external panel)

**Panel:** ≥ 6 raters — cardiologists, accredited sonographers, echo educators.
**Design:** each rates a fixed set of scenarios (5 standard views × normal + 8
pathologies + 4 rhythms), blinded to development. **Instrument** (5-point Likert,
1 = unrealistic … 5 = indistinguishable from clinical):
1. Anatomical accuracy of each standard view.
2. Realism of B-mode texture (speckle, resolution, artifacts).
3. Correctness & realism of colour/spectral Doppler.
4. Physiological correctness of the cardiac cycle / wall motion / PV behaviour.
5. **12-lead ECG realism and correct pathology/rhythm signature.**
6. Correctness of measurements & severity grading vs guidelines.
7. Educational value / curricular fit.
**Analysis:** item means + 95 % CI, inter-rater reliability (ICC / Krippendorff's
α), free-text thematic analysis → prioritised defect list. **Pass bar:** median
≥ 4 on every item, no item with > 20 % ratings ≤ 2.

---

## C. Learning-transfer study (external, gold standard)

**Design:** prospective **randomised controlled trial**. Novice trainees →
EchoSim training arm vs conventional-teaching control.
- **Primary outcome:** performance on a **real scanner / high-fidelity mannequin**
  OSCE (view acquisition, image optimisation, interpretation) scored by blinded
  experts.
- **Secondary:** time-to-competency, knowledge test, retention at 3 months,
  eye-tracking gaze efficiency, learner confidence.
- **Statistics:** power for a moderate effect (d≈0.6) → ~30/arm; ANCOVA on
  post-scores with baseline covariate; pre-registered.
- **Framework:** map evidence to **Kane's inferences** (scoring → generalisation
  → extrapolation → implications). This is the evidence CAE Vimedix / HeartWorks
  published and the bar for commercial credibility.

---

## Current status
- **A (automated):** implemented in-app (Validation panel) and in
  `tools/validate.mjs` — **runs today**.
- **B, C (human/RCT):** protocol + instrument + statistical plan **shipped here**;
  execution requires an external panel/ethics approval/participants — scoped and
  turnkey, but not completable by software alone. This is the honest residual on
  benchmark **B9**.
