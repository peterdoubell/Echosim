# EchoSim — Roadmap from 7.3 to 9.9 (and to clinical-training deployment)

This document is the honest, actionable plan to take EchoSim from its current
**weighted 7.3 / 10** (v2 "state-of-the-art" rubric, see `../BENCHMARKS.md`) to
the **9.9 / validated-clinical-training** bar. It separates what **code** can do
from what only an **external evidence-and-audit process** can do, because the two
are fundamentally different kinds of work and conflating them is how simulators
overclaim.

## Where 7.3 comes from
Every dimension was driven up with verified, physics-grounded implementations
(12-lead ECG, PV-loop-coupled mechanics, continuity Doppler, scatterer imaging,
graded pathology, learner Simpson-EF, TEE, myofibre + chordae anatomy, GLS).
Six axes sit at 7.5, quantification at 8.0, UX/engineering at 8.5. The two that
hold the weighted score down are **B7 (6.5)** and **B9 (4.5)** — and B9 is the
decisive one.

## Two categories of remaining work

### A. Codeable — large physics builds (raise the ceiling to ~8.0–8.3)
These are real, achievable, and each is a **multi-day-to-multi-week** effort for
roughly **+0.05–0.15 weighted** apiece. None, alone or together, reaches 9.9 —
the weighted ceiling with B9/regulatory unchanged is ~8.3.

| Build | Axis (wt) | Payoff | Notes / empirically-found constraints |
|---|---|---:|---|
| Reduced-order CFD / FSI (Navier–Stokes or potential+jets with a real pressure field) | B4 (1.1) | 7.5→~8.5 | Largest single physics lift; validate jets within ±8% of the analytic continuity values already in the model. |
| FE hyperelastic mechanics (Holzapfel–Ogden) with fibre strain + regional/segmental strain bullseye | B3 (1.1) | 7.5→~8.5 | The GLS field (`longitudinalStrain`) is the seed; a 17-segment strain map showcases RWMA. |
| Eikonal activation on a fibred mesh → torso lead-field 12-lead | B2 (1.2) | 7.5→~8.5 | Replaces the VCG→inverse-Dower forward model; validate vs MedalCare-XL. |
| Field II-class scatterer sim + CAMUS speckle/segmentation match | B5 (1.1) | 7.5→~8.5 | Needs the external **CAMUS** dataset for the comparison (semi-external). |
| Statistical-shape morphology variants + trabeculation | B1 (1.0) | 7.5→~8 | Fibre field + chordae already done; add ≥1 alternative geometry. |
| MPR from the SDF volume + tissue-Doppler | B7 (0.8) | 6.5→~7.5 | **Empirical caveat:** a naïve TDI shows e′<a′ for a normal heart (impaired-relaxation artifact) because the diastolic filling waveform needs re-tuning — do the re-tune *before* shipping TDI, or it is negative training. |
| Learner VTI + continuity-AVA / DVI | B8 (1.0) | 8.0→~8.5 | **Empirical caveat:** a clinically-consistent AVA needs AS forward SV ≈ 90 mL; the calibrated model produces ~50 mL (a low-flow state), so continuity-AVA reads discordantly "severe" for mild AS. Requires an AS-hemodynamics re-tune (raise compensated SV without breaking the PV loop) as a prerequisite. |

**Empirical boundary log (attempted at the code level this program):**
1. *Tissue-Doppler* — model gives correct s′ (0.11 m/s) but e′<a′ for normal;
   shipping it would teach a false diastolic pattern. Not shipped.
2. *Continuity-AVA re-tune* — raising AS Emax lifted SV 50→60 and EF to a more
   physiologic 68–70%, but velocities fell out of band and a clinical AVA still
   needs SV≈90; reverted to preserve the verified calibration.
3. *Weight arithmetic* — the original rubric divided Σ(w·score) by 9.6 when the
   weights sum to 10.6; corrected (baseline 6.6→5.9).
4. **Diastolic-filling calibration — ✅ RESOLVED (commit `7256179`).** Was A-wave
   dominant (normal E/A 0.86). Fixed by trimming the a-wave (`aWave` 5→4.4) and
   raising resting LA pressure (`Pla0` 7→8.0): normal filling is now **E-dominant
   (E/A 1.15)**, DCM 0.52 (impaired relaxation), MR 1.45 — and **every case's EF/PV-loop
   is preserved** (normal 59, AS 57, DCM 20, MR 90, MS 58; 12/12 ECG still pass). This
   unblocked and shipped the live **mitral E/A** readout. Residual: **tissue-Doppler
   e′/a′** still needs the annular-velocity timing checked, and **MAPSE/TAPSE** still
   needs #5.
5. **Apex-anchored contraction — ✅ RESOLVED (commit `e549b56`).** The LV now
   contracts toward a fixed apex (centre shifts apically as the long axis shortens),
   so the mitral annulus descends the FULL shortening. Implemented volume-preservingly
   (only the centre moves; radii unchanged → EF/LVIDd/LVIDs intact: normal EF 59, DCM
   20, AS 57). The AV valve planes + papillary muscles follow the same material map
   (verified the valve tracks the descending annulus in A4C diastole vs systole), and
   the coarse proxy shifts with the SDF. Shipped the **MAPSE** readout (flagged <10 mm):
   normal 14 mm, DCM 6 mm, RWMA 11 mm. This also moves absolute TDI e′ toward physiologic
   and is the geometry basis for TAPSE. Residual: expose TAPSE (RV annulus, same idea)
   and absolute e′/E-over-e′ TDI numbers.
6. **Valve + subvalvular apparatus — ✅ RESOLVED.** The valves were a flat annular
   disc that irised open concentrically, and there was no subvalvular apparatus at
   all. Now: leaflets hinge at the annulus and swing as curved, tapered-thin
   curtains (S-profile, thickest at the base, fine free edge); the mitral is
   bileaflet-asymmetric with a scalloped (P1/P2/P3) posterior leaflet; chordae
   tendineae fan from the papillary tips to the leaflet free edges (taut in
   systole, slack in diastole, fused in rheumatic MS); the tricuspid has its short
   septal leaflet and an RV anterior papillary continuous with the moderator band;
   and the semilunar valves are individuated into **three cusps** with a trilobed
   orifice and commissural coaptation seams — the short-axis "Mercedes" Y. Shipped
   with a **PSAX-AV** view (cusps named R/L/N by anatomy, not by index) and a
   matching 3D subvalvular overlay.
7. **Atrial anatomy — ✅ RESOLVED, and it surfaced a real defect.** Probing the
   model showed the **LA and RA blood pools abutted directly with no interatrial
   septum**: every simulated heart in effect had a wide-open ASD. The septum is now
   explicit — a disc normal to the LA→RA axis, thin over the **fossa ovalis**
   (1.1 mm) and thicker at the limbus (2.7 mm) — and the ASD pathology perforates
   the fossa, so the defect sits exactly where secundum ASDs occur. Also added:
   a narrow-ostium multi-lobed **chicken-wing LAA** (vs the broad-based triangular
   RAA) with a dedicated **ME LAA** TEE view, the **crista terminalis**, the
   **Eustachian valve**, SVC/IVC stubs, the four **pulmonary veins**, and the
   LSPV/LAA ("warfarin") ridge. Total per-pixel cost held to **+7 %** by gating the
   atrial sub-anatomy behind bounding tests.
8. **Imaging fidelity — ✅ PARTLY RESOLVED (step two of the expert's sequence).**
   The PSF was monotonic in depth — narrowest at the transducer face, which no
   focused beam is. It is now narrowest **at the focus** and diverges either side,
   with a user-settable focal depth and the caret scanners draw beside the depth
   ruler; this alone lifted fundamental A4C gCNR **0.963 → 0.986**. **Tissue
   harmonic imaging** ships as a toggle (narrower effective beam, reverberation
   suppression, and an honest *penetration penalty*), improving gCNR a further
   0.971 → 0.985 (PLAX) **as a consequence of the narrower beam, not by
   construction**. Persistence and an **honest measured frame rate** are exposed.
9. **Slice thickness, two-way beam, PV Doppler — ✅ RESOLVED.** The renderer
   sampled an infinitely thin plane through a transmit-only beam; both flattered
   the image. **Elevational slice thickness** is now finite — a slab set by a
   FIXED (non-steerable) acoustic lens, ~1.8 mm at the elevation focus widening to
   ~5–6 mm either side — sampled stochastically for genuine partial-volume
   averaging. It pulls gCNR from an unrealistically clean ~0.98 down to **~0.91**,
   into the 0.7–0.9 band real myocardium-vs-cavity contrast occupies: the
   simulator had been *easier than real echo*. The beam is now the **round trip** —
   a once-focused transmit multiplied by a dynamically-refocused receive limb at
   constant F-number — so the in-focus receive beam bounds off-focus degradation
   (lateral radius 3→2 at the transducer face, 7→5 at 20 cm vs transmit-only),
   which is why real scanners stay usable away from focus. Axial resolution is
   decoupled from focus entirely (pulse-length limited, degrading with depth only
   through attenuation-driven downshift). Also shipped **pulmonary-vein Doppler**:
   S/D waves with S ≥ D, atrial reversal, and **systolic flow reversal in severe
   MR** — an ASE-specific criterion, correctly grade-dependent.
   Residual: decoupling model rate from render rate; 3-D / biplane acquisition.
These are why the codeable increments below 8 are "re-calibrate-first," not quick wins.

### B. Non-codeable — the actual gate to 9.9 / clinical deployment
The 9–10 band is a **validated digital twin**, and "ready for commercial
deployment for clinical training" is a **regulatory/validation claim**. Neither
can be produced by writing code:

1. **B9 external validation** (weight 1.0, currently 4.5). Requires, per
   Messick/Kane (see `VALIDATION.md`):
   - a **face/content-validity expert panel** (≥6 cardiologists/sonographers,
     the Likert instrument is shipped) — needs people + time;
   - **construct validity** (novice-vs-expert discrimination) — needs subjects;
   - **concurrent validity** vs an independent phantom/ground truth — needs a
     physical phantom;
   - **predictive/transfer validity** — a **randomised learning-curve trial**
     (~30/arm, pre-registered, IRB-approved) — needs ethics approval,
     participants, months of follow-up.
2. **Regulatory / QMS** (see `COMMERCIAL_READINESS.md`): ISO 13485 QMS,
   IEC 62304/62366 lifecycle + usability files, ISO 14971 risk file, and an
   organisational audit. The education-only positioning keeps this out of
   device regulation, but a *clinical-training* deployment claim still needs the
   QMS + usability evidence.

**These are audit-and-evidence processes with people in the loop. No commit
closes them.** Reporting 9.9 without them would violate the product's own
ISO 14971 risk control against miscalibration/overclaiming.

## The honest bottom line
- **Deployable today:** as an **education-only, non-device** trainer (offline, no
  PHI) at a verified **7.3/10** — a genuinely physics-grounded tool.
- **To ~8.3 (code):** execute the Category-A builds above; each is scoped and
  achievable, none reaches 9.9.
- **To 9.9 / validated clinical-training deployment:** the Category-B external
  validation study + QMS/regulatory evidence — the only path, and it is funded
  human work, not engineering.

If you are picking this up: choose **one Category-A build** (recommend B4 CFD or
B3 FE-strain for the biggest lift) and I will execute it as a dedicated effort,
labelled for exactly what it delivers — or stand up the **Category-B** package,
for which everything code can contribute (protocol, instrument, analysis harness,
CAMUS comparison tooling) is already in the repo.
