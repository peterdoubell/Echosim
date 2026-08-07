# EchoSim — Expert Validation Results (B9, in progress)

This file records the **executed** expert review — the live evidence for benchmark
**B9 (external / expert validation)**, whose protocol and instrument are in
`VALIDATION.md §B`. It is updated as expert ratings and defect reports come in.

> **Status:** expert review **in progress** with a **cardiac radiologist / anatomist**
> (self-identified). This is genuine content-validity evidence and moves B9 off its
> internal-only floor. Full B9 (→ higher band) still requires the **≥6-rater panel**,
> a **learning-transfer RCT**, and **dataset agreement** (CAMUS/MedalCare-XL) per the
> protocol — a single expert's face/content review is the first, not the last, step.

---

## Panel
| Rater | Role | Domain strengths | Session |
|---|---|---|---|
| R1 | Cardiac radiologist / anatomist | Anatomical/structural accuracy vs CT/MRI ground truth | current |

## Review sequence (as directed by R1)
R1 set the order: **anatomy first** — "representative and to-scale rather than
abstract morphology, with state-of-the-art modality-like pixel and temporal
resolution/fidelity" — then the remaining domains in sequence.

## Domain 1 · Anatomy & imaging resolution

### Changes made in response (committed)
| # | Change | Evidence | Commit |
|---|---|---|---|
| 1 | B-mode formation buffer 240×300 → **336×420** (~0.36 mm/px axial), 3-tier adaptive sampler, user **Detail** control | gCNR 0.98 / speckle-SNR 5.8 preserved; PSF rescaled to hold physical beam width | `cd72242` |
| 2 | **Bullet-shaped LV** (apical taper) replacing the symmetric ellipsoid | EF/LVIDd unchanged (58 %, 4.6 cm); 10 pathologies no-throw | `7a55a90` |
| 3 | **Triangular RV** apex (matches the bullet LV) in A4C | PSAX crescent + papillaries preserved; 12/12 ECG; 0 console errors | `1b0572e` |
| 4 | **AHA 17-segment regional-strain bullseye** + GLS panel | normal −20 % uniform, DCM −8 % global, RWMA-septal septal −6 %/lateral −20 %; 0 errors | `a4ad354` |
| 5 | **Diastolic-filling fix** (E-dominant) + live **E/A** | normal E/A 1.15, DCM 0.52; every EF/PV-loop preserved | `7256179` |
| 6 | **Apex-anchored contraction** + **MAPSE** (14/6/11 mm normal/DCM/RWMA) | EF/LVIDd intact; valve tracks annulus; 12/12 ECG | `e549b56` |
| 7 | **PASP** from TR jet (PH 84 / TR 41 mmHg) + **D-sign** PH pathology | 0 errors | `af6a0e7`,`f118d90` |
| 8 | **UX**: collapsible console panels (−33% height, persisted) | 0 errors | `40f587d` |

### R1 ratings (to be entered from the review)
| Item (VALIDATION.md §B, 1–5) | Score | Notes |
|---|---:|---|
| 1. Anatomical accuracy of standard views (PLAX/A4C/PSAX) | _pending_ | |
| 2. B-mode texture realism (speckle/resolution/artifacts) | _pending_ | |
| 6. Correctness of measurements & severity grading | _pending_ | |
| 7. Educational value / curricular fit | _pending_ | |

### R1 defect list (to be entered)
- _pending — specific dimensional corrections (cm targets) and morphology priorities:
  triangular-RV moderator band, D-shaped septum in PSAX, mitral leaflet/chordae/papillary
  detail, atrial appendage shaping._

## Domains 2–7 (Doppler, ECG, physiology, measurements, teaching)
_Not yet reviewed — scheduled after anatomy is signed off by R1._

---

### How scores feed the benchmark
B1/B5 improvements are **implemented and objectively verified** (resolution,
morphology). The **B1 score is deliberately held at 7.5** until R1's face-validity
rating is entered here — self-scoring one's own morphology is the internal-rating
bias the v2 rubric exists to correct. B5 is raised to 8.0 on the **objective**
resolution gain (measured mm/px + preserved gCNR), independent of subjective rating.
