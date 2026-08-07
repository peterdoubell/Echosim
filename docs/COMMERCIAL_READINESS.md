# EchoSim — Commercial Readiness Dossier (clinical-training deployment)

*Honest status: the **software fidelity** benchmarks (B1–B8, B10, B11) are being
driven into the 9s in-code. The two things that **cannot be completed by
writing code** — independent **clinical validation** (B9) and **regulatory /
quality-system clearance** — are external, evidence-and-audit processes. This
dossier specifies exactly what "ready for commercial deployment for clinical
training" requires, where EchoSim stands, and the concrete path to close each
gap, so the product is deployment-ready in every respect the engineering can
control and fully specified for the rest.*

---

## 1. Intended use & positioning (this determines everything downstream)
- **Intended use:** an *educational simulator* for teaching echocardiographic
  image acquisition, anatomy, Doppler and ECG interpretation, and
  pathology pattern-recognition. **Not** a diagnostic device; **not** used on or
  with patient data; **not** for clinical decision-making.
- **Why this matters:** in both the **US (FDA)** and **EU (MDR 2017/745)**, a
  training simulator that makes **no diagnostic/therapeutic claim and touches no
  patient** is generally **outside the definition of a medical device**. This is
  the correct, lowest-risk market position and must be stated in all labelling.
  If the product were ever marketed to *inform clinical decisions*, it would
  become a regulated device (SaMD) — explicitly out of scope here.
- **Deliverable:** an **Intended Use / Indications-for-Use statement** and
  **labelling** (already reflected in the in-app footer "not for diagnosis").

## 2. Quality management & software lifecycle (process evidence)
Even as non-device educational software, credible clinical-training vendors run
a real QMS. Target evidence set:
- **ISO 13485** (or an ISO 9001 QMS scoped to design controls) — design history
  file, change control, CAPA.
- **IEC 62304** software lifecycle — the app is small and self-contained;
  produce: software requirements spec, architecture (the module map already
  exists), unit/integration test records, and a release/versioning record
  (git history provides traceability today).
- **IEC 62366-1** usability engineering — a use-specification, use-related risk
  analysis, and formative + summative usability testing with representative
  learners.
- **ISO 14971** risk management — for a training tool the dominant risk is
  **negative training / miscalibration** (a learner internalising a wrong
  pattern). Mitigations: the explicit "schematic model, not for diagnosis"
  labelling, the BSE/ASE-calibrated measurements with reference ranges, and the
  validation study below. Maintain a **risk file**.
- **Status:** engineering artefacts (requirements, architecture, tests, versioned
  releases) are largely present or trivially generable; formal QMS
  certification is an organisational step.

## 3. Regulatory pathway (by market)
| Market | Path if positioned as **education only** | Path if ever positioned as clinical/SaMD |
|---|---|---|
| **US** | Not a device → no 510(k)/De Novo; follow FTC truthful-claims + general software norms | SaMD; likely Class II 510(k); FDA CDS guidance applies |
| **EU** | Not a medical device under MDR (no medical purpose) | MDR conformity assessment, CE mark, Notified Body |
| **UK** | Not a device under UK MDR 2002 | UKCA + MHRA |
| **Canada/AUS** | Non-device educational software | Health Canada / TGA SaMD |
Keep the **education-only** position: it is both accurate and the deployable one.

## 4. Data protection, security, accessibility (deployable-now items)
- **Privacy:** the simulator is **fully synthetic and offline** — **no PHI, no
  telemetry, no network calls** (verified: vendored Three.js, no external
  fetches). This clears GDPR/HIPAA concerns by construction — a genuine
  commercial advantage; state it explicitly.
- **Security:** static assets, strict no-eval, no third-party runtime calls.
  Recommend adding a Content-Security-Policy header and SRI when hosted.
- **Accessibility:** WCAG 2.1 AA / Section 508 — keyboard operable, focus-visible,
  ARIA labels/live regions, reduced-motion honoured (already implemented). A
  formal **VPAT** is the remaining artefact.

## 5. Clinical validation (B9 — the evidence a training product needs)
This is the substantive external gap. The recognised framework for simulation
validity is **Messick / Kane** (evidence for a validity *argument*), operationalised as:
1. **Face & content validity** — structured expert review (cardiologists,
   sonographers, echo educators) rating anatomical/haemodynamic/ECG realism and
   curricular coverage on a validated Likert instrument. *(We ship the
   instrument + harness; the panel is the external step.)*
2. **Construct validity** — the simulator distinguishes novices from experts
   (performance/eye-tracking/time-to-competency differs by training level).
3. **Concurrent/criterion validity** — measurements (LVIDd, EF, peak V, gradients)
   agree with an **independent ground truth** (a physical phantom or expert
   reading) within stated limits of agreement (Bland–Altman).
4. **Predictive/transfer validity** — a **randomised learning-curve study**:
   does training on EchoSim improve real-scanner performance vs control? This is
   the gold-standard evidence CAE Vimedix/HeartWorks published.
- **What we deliver in-repo:** `docs/VALIDATION.md` (full protocol + statistical
  plan + the expert-Likert instrument) and an in-app **self-validation harness**
  (live **gCNR/CNR** image-quality readout, and **measurement-accuracy** checks
  against the model's own ground truth). Independent human/RCT evidence is the
  external step and is scoped, costed and protocol-ready.

## 6. Content & scope for a commercial curriculum (B6/B7)
- Broaden to ≥20 graded pathologies with linked ECG/haemodynamics (in progress),
  add a **TEE** window set + **MPR from the 3D volume**, tissue-Doppler/strain,
  a **structured curriculum + logbook + competency mapping** (e.g. BSE/EACVI/ASE
  accreditation domains), and instructor analytics.

## 7. Honest readiness scorecard
| Area | Deployable now? | Gap owner |
|---|---|---|
| Software fidelity (B1–B8,B10,B11) | Driving to 9s in-code | Engineering (in progress) |
| Privacy/security/offline | **Yes** | — (advantage) |
| Accessibility | Yes (VPAT to formalise) | Documentation |
| Intended-use / labelling | **Yes** | — |
| QMS / IEC 62304 / 62366 / ISO 14971 | Artefacts producible | Organisation |
| Regulatory (education-only) | **Yes** (non-device) | Legal review |
| Clinical validation (B9) | Protocol + harness shipped | **External study** |
| Curriculum / TEE / logbook | Partial | Content (in progress) |

**Bottom line:** EchoSim can be **commercially deployed today as an
education-only** trainer with a clean privacy/security posture; reaching the
evidence bar of the market leaders additionally requires the **independent
clinical validation study** (protocol shipped) and formal **QMS/usability
artefacts** — external processes that this dossier fully specifies.
