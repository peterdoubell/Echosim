# EchoSim — Anatomical Model and Verification

This document describes how the heart is built, how it is checked, and what an
internal expert-style review found. The numbers below come from
`node tools/verify-anatomy.mjs` at the commit that added this file.

## Construction

The heart is a signed-distance field assembled from **landmarks** (`LM` in
`js/anatomy.js`) placed at adult ASE/EACVI 2015 dimensions. Every standard view
(`js/views.js`) is derived from the same landmarks and displayed in the ASE
orientation, so geometry and views cannot drift apart.

Motion comes from a lumped-parameter circulation (`js/hemodynamics.js`):
time-varying elastance for both ventricles, and flow-driven valve opening.
Contraction is apex-anchored, so the base descends (MAPSE / TAPSE). The atria
fill as reservoirs, empty in conduit and booster phases, and hang between their
moving annuli and roofs.

Key modelling choices, each made to match what a sonographer sees:

| Structure | Model |
|---|---|
| LV | Half-prolate body with basal shoulder, mitral and LVOT funnels, apical thinning and trabeculation; remodelling from the modelled EDV |
| RV | Crescentic body with inflow, RVOT and apex taper; moderator band in the four-chamber plane, anterior papillary muscle |
| Atria | Ellipsoids on their annuli, stretched to one shared **interatrial septum** (limbus ~3 mm, fossa ovalis ~1 mm; secundum ASD = real defect); LAA on the anterior LA wall (seen in A2C), broad-based RAA, pulmonary veins, SVC/IVC, crista, Eustachian valve, coronary sinus |
| Mitral / tricuspid | Hinged leaflets on a D-shaped, saddle-shaped annulus: coaptation ~0.45 cm below the hinge line with ~0.5 cm apposition; fish-mouth orifice in PSAX-MV; MS commissural fusion and doming; MR P2 prolapse graded with severity; chordae to the live free edges |
| Tricuspid annulus | Hinges: lateral rim moves by TAPSE, septal rim ~45 % of it |
| Aortic valve / root | Three cusps in their sinuses (RCC anterior, LCC left, NCC toward the septum), hammock-shaped when closed (Y closure in PSAX-AV), sinuses, STJ, ascending aorta, arch, descending aorta |
| Papillary muscles | Wall-based cones at ~4 and ~8 o'clock, free tips under the commissures |
| Coronaries | LM, LAD, LCx, RCA and PDA as smooth tubes in their grooves |
| Pericardium & surroundings | Parietal pericardium (not between the atria), anterior epicardial fat pad, transverse sinus; a continuous diaphragm whose height field rises to meet the pericardial sac (one fused bright line on subcostal views); an ellipsoidal liver with IVC and hepatic veins; aerated lung with parasternal and apical windows; parasternal views look through ~1.6 cm of chest wall |

## Image formation (`js/echo.js`)

- Speckle comes from a complex scatterer field fixed to the tissue (material
  coordinates that undo the ventricular contraction and annular descent). It is
  convolved with the same PSF and envelope-detected, so the Rayleigh grain moves
  with the walls and is elongated laterally along the depth arcs. Blood
  scatterers are redrawn every frame.
- Ventricular myocardial backscatter depends on the fibre angle (~6 dB).
- Backscatter is specified per tissue in dB: blood −60, myocardium −28,
  valve −15, calcium −4, pericardium −10 plus a specular term. It is shown with
  log compression over 55 dB.
- Elevation slice thickness uses a stochastic slab that is resolved before the
  PSF, and TEE uses a thinner slab.
- The PSF is about 2.2 mm laterally at focus and wider elsewhere; axially it is
  about 0.7 mm.
- Attenuation accumulates by tissue along each beam and is offset by TGC, which
  produces posterior enhancement. Calcified valves cast acoustic shadows and
  normal valves never do.
- Specular interface echoes depend on incidence angle, so walls parallel to the
  beam drop out.
- Lung appears as a pleural line with A-lines.

## Quantitative audit — 176 of 187 checks pass, 11 tagged known open

Selected results (normal case, ED = phase 0, ES = minimum LV volume):

| Measurement | Model | Adult reference |
|---|---:|---|
| LVIDd / LVIDs (PLAX), FS | 4.88 / 3.53 cm, 27.7 % | 4.2–5.8 / 2.5–4.0, 25–45 % |
| IVSd / PWd | 0.89 / 0.87 cm | 0.6–1.0 |
| EDV / ESV / EF (voxel) | 114 / 49 mL / 56 % (circulation 121 / 47 / 61 %) | 62–150 / 21–61 / 52–72 |
| LV length, sphericity | 8.35 cm, 1.71 | 7.2–9.4, 1.5–2.2 |
| Aortic annulus / sinus / STJ | 2.41 / 3.03 / 2.94 cm | 2.0–2.9 / 2.9–3.7 / 2.4–3.2 |
| LA AP (PLAX) / LA volume | 3.32 cm / 67 mL | 3.0–4.0 / 30–70 |
| RA major / minor | 4.97 / 3.39 cm | 3.4–5.3 / 2.6–4.4 |
| RVD1 / RVD2 / RVD3 | 3.38 / 3.38 / 6.39 cm | 2.5–4.1 / 1.9–3.5 / 5.9–8.3 |
| MAPSE / TAPSE | 1.70 / 2.09 cm | 1.0–2.0 / 1.7–2.8 |
| Mitral tenting height (ES) | 0.46 cm | < 0.6 |
| Fossa ovalis thickness (ED / ES) | 0.10 / 0.11 cm | 0.1–0.2 |
| TV septal offset | 0.85 cm | 0.3–1.2 |
| Aortoseptal angle | 131° | 120–150 |
| First cardiac structure, parasternal | 1.6 cm | ≥ 1.5 |
| Echo-free gap liver → heart (subcostal) | 0.26 cm | ≤ 0.3 |

**Views (21):**
- TTE: PLAX, PSAX (mid), PSAX-MV, PSAX-AV, A4C, A5C, A2C, A3C, RV inflow,
  subcostal 4C, subcostal IVC and suprasternal arch (with the head-and-neck
  branches).
- TEE: ME 4C, ME 2C, ME LAX, ME AV SAX, ME bicaval, ME RV inflow–outflow, ME
  LAA, TG mid SAX and descending-aorta SAX.

All are derived from landmarks.

View checks cover orientation and content for PLAX, PSAX (papillary muscles at
4:00 and 8:00), PSAX-AV, PSAX-MV (no atria at any phase, leaflets in diastole),
A4C (not foreshortened, no LVOT), A2C, A5C, subcostal and TEE ME4C / ME LAX.

**Close to the limits:** LA minor 4.79 cm (upper limit 4.8), RVD2 3.38 cm
(upper limit 3.5) and LA volume 67 mL (upper limit 70). All pass, but a small
geometric change could push them out of range. The voxel LV volume is ~6 %
below the circulation's EDV, so the geometric EF (56 %) reads lower than the
circulation's (61 %).

## Internal expert-panel review

The model at commit `6a57fb8` went to a multi-agent expert panel (AI reviewers
role-playing a senior sonographer and an ultrasound physicist, each with
adversarial verification). This is **not** external validation. Two of the
planned eight lenses completed; the rest were cut off by a usage limit.

| Lens | Grade at 6a57fb8 | Main findings | Status now |
|---|---:|---|---|
| TTE views (sonographer) | 5 / 10 | Mitral/tricuspid tented 1.5–2.7 cm; aortic cusps rotated ~60°; no interatrial septum; PSAX-MV at AV-groove level; PLAX plane off the mitral centre; RA sweeping into mid-PSAX; moderator band out of plane; papillary muscles as detached rods; no chest wall; subcostal false-effusion gap; LAA not in any TTE view; aortic "brick"; SVC pseudo-mass | All addressed and covered by audit checks or visual checks |
| Image physics | 5 / 10 | No thin IAS; normal cusps shadowing; TEE dither hatch; plastic texture and too-fine lateral resolution; no attenuation / linear-ish grey map; outline artefacts; no lung; no chest wall; enhancement streaks; bright AV-groove fat | All addressed (see `js/echo.js`), including fibre anisotropy and tissue-locked speckle |

Neither lens has re-graded the current model. The grades above describe the
earlier snapshot, and nobody has independently confirmed that the fixes raise
them.

### Third review (commit `b9fa357`)

All six lenses completed, and each finding went to adversarial verifiers.
15 findings were confirmed and 5 were refuted. Grades out of 10:

| Lens | Grade |
|---|---:|
| TTE views | 7 |
| Motion / timing | 7 |
| 3D anatomy | 7 |
| Image physics | 6.5 |
| Pathology | 6.5 |
| TEE / subcostal | 6 |

The panel counts 8/10 as commercial clinical-training grade, so the model is
not there yet.

Fixed since that review (each is covered by an audit check or by the
haemodynamic summary):
- **A3C:** now in the apex–mitral–aortic plane, so the RV drops out.
- **PSAX-MV:** now at the leaflet tips, showing the fish-mouth orifice in the
  normal heart and in MS.
- **Grey scale:** a narrower, weaker specular lobe; mediastinal soft tissue
  darker than myocardium; posterior enhancement capped at 6 dB.
- **SSN:** a longer PA trunk and a raised arch, so the RPA crosses beneath the
  arch above the LA.
- **Aortic stenosis:** the mild, moderate and severe grades give 2.7, 3.5 and
  4.5 m/s, and the Doppler peak comes from the modelled gradient.
- **DCM:** now has functional MR (regurgitant fraction ~29 %).
- **Pulmonary hypertension:** the RV free wall is hypertrophied.
- **Atria:** the drawn atria empty physiologically (voxel LAEF 49 %,
  passive 33 %, booster 24 %).

Still open:
- **ME bicaval (critical):** the aortic root intrudes into the plane.
- ME LAA morphology (no LUPV / coumadin ridge).
- Subcostal IVC obliquity.
- LA outline lobulation.
- Right-heart valve timing identical to the left heart.
- The crux thickening in systole.

### Fourth review (commit `1ae49a7`) and the heart-position rework

Grades out of 10:

| Lens | Grade |
|---|---:|
| TTE views | 7 |
| Motion | 6.5 |
| 3D anatomy | 6 |
| Image physics | 6.5 |
| Pathology | 6.5 |
| TEE / subcostal | 6 |

A usage limit cut the adversarial verification short, so these findings are
unverified.

The 3D-anatomy lens marked the heart's rotation about its long axis as
critical. The heart is now rolled 45° into its true position:
- **Chambers:** the RA is level with the LA and anterior-right of it, and the
  RV is anterior.
- **Valves:** the pulmonary valve is anterior, superior and to the left of the
  aortic valve.
- **Great vessels:** the ascending aorta runs up just right of the sternum.
  The SVC is right-posterior of the ascending aorta, and the caval orifices sit
  on the posterior RA.

Body-referenced parts were made robust to the orientation:
- a near-horizontal diaphragm whose height comes from the heart, with a
  cleft-bridging lift;
- fat in the interatrial groove;
- pulmonary veins and the aortic-root frame that stay rigid with the heart.

Other fixes:
- **ME bicaval:** now passes through both caval orifices, with no aortic root.
- **A4C:** aimed just posterior to the crux, so the false crux "mass" is gone.
- **Lung:** gated per beam, giving a smooth pleural line with no dither
  ladder.
- **TEE near field:** the fixed reverberation bands are removed.
- **Mitral stenosis:** Doppler follows the grade.
- **Tricuspid septal offset:** now persists through systole.

### Fifth and sixth reviews (commits `53c64a6`, `0b4bddf`)

Grades out of 10:

| Lens | Review 5 | Review 6 |
|---|---:|---:|
| TTE views | 6.5 | 6.5 |
| Image physics | 6 | 6 |
| TEE / subcostal | 5.5 | 5.5 |
| Motion | 6.5 | 7 |
| 3D anatomy | 6 | 6.5 |
| Pathology | 6 | 6.5 |

A usage limit cut the adversarial verification short in both reviews, so
their findings are partly unverified.

Fixed since:
- **Atria:**
  - The atria empty by free-wall motion toward a fixed septum, so the IAS is
    thin at ED and ES; both phases are audited.
  - Groove fat appears only in the true interatrial groove.
  - Drawn LAEF is 49 %, and the RA keeps a chamber at ED.
- **Venous tubes:** the IVC and pulmonary veins keep their lumen, with no solid
  plugs. The RPA runs beneath the arch without being cut (audited).
- **Subcostal views:**
  - Subcostal 4C passes through the sub-xiphoid beam and the apex, with the
    septum across the beam.
  - Subcostal IVC runs across the image into the RA.
- **PSAX:** the standard orientation, with the papillary muscles at 4 and
  8 o'clock.
- **TEE:**
  - The probe sits against the LA wall.
  - The ME bicaval view shows the LA near field, the IAS and the RA with the
    SVC.
- **Pathology:**
  - Atrial remodelling: giant LA in MS (~120 mL), MR and DCM; dilated RA in
    TR, PH and ASD.
  - Nodular calcific AS cusps.
  - Domed MS leaflets.
  - Eccentric P2-prolapse MR jet.
  - The VSD is seen in both PSAX and A4C.
  - MS Doppler follows the grade.
- **Motion:** TAPSE 2.0 cm, above MAPSE.
- **B-mode:**
  - Lung reads as air.
  - Pericardium and vessel walls are dimmer (no angle-independent white
    outlines).

Still open:
- A too-fast aortic upstroke. Adding inertance destabilised the lumped model,
  so it was reverted.
- Liver visible in the parasternal far field.
- A boxy suprasternal arch.
- The TEE ME4C RA is small.
- No IVC plethora and no TR annular dilatation.
- The aortic-root frame is heart-rigid, so the root stays more horizontal than
  real.

### Seventh and eighth reviews (commits `e6751e1`, `70963c9`)

Grades out of 10:

| Lens | Review 7 | Review 8 |
|---|---:|---:|
| TTE views | 6.5 | 6.5 |
| Image physics | 6.5 | 7 |
| TEE / subcostal | 6 | 6 |
| Motion | 6.5 | 6.5 |
| 3D anatomy | 6.5 | 6 |
| Pathology | 7 | 7 |

Review 8 verified all 40 of its findings.

Fixed since:
- **Atria:**
  - Emptying is split between free-wall collapse toward the septum and roof
    descent toward the annulus. The ED atria stay rounded chambers, with no
    slit and no pericardial slab flashing inside the LA in A2C.
  - The RA body lies right-posterior of the aortic root; before, the root
    carved it into a crescent.
  - A dilating atrium balloons within the four-chamber plane.
- **Caval system:**
  - The sinus venarum joins the caval orifices, so SVC → RA → IVC is one
    channel.
  - The caval junctions are tethered to the mediastinum.
  - ME bicaval, subcostal IVC and RV inflow show the cavae entering the RA.
- **Planes:**
  - ME4C passes through both AV valves and the apex.
  - ME bicaval passes through both cavae and the fossa.
  - ME LAA passes through the LAA neck and the LSPV.
  - ME RV inflow–outflow passes through the TV and PV.
  - Every ME probe keeps the oesophageal wall in front of the LA.
  - PSAX-MV, PSAX-AV and ME AV SAX follow the annulus or root through the
    cycle.
- **Remodelling:** a dilated LV grows toward the apex, with the fibrous
  skeleton fixed. The apical windows follow the live apex.
- **Great vessels:**
  - Candy-cane arch with the RPA snug beneath it.
  - The RVOT moves as a unit with the pulmonary valve, so it stays in front
    of the root in systole.
- **Circulation:**
  - Aortic outflow has inertance and characteristic impedance: peak flow
    about 545 mL/s at about 85 ms, pressure 131/70 mmHg, LVET 288 ms.
  - The A wave ends at mitral closure: IVCT 64 ms, IVRT 77 ms, E/A 1.2.
- **B-mode:**
  - The pleural echo depends on incidence angle and is blurred by the pulse.
  - Near-field reverberation scales with the chest wall, so there is none in
    TEE.

### Fixes after review 8

These fixes have **not** been independently re-reviewed. The grades above
stand until a further panel review.

- **Audit:** `tools/verify-anatomy.mjs` now samples any view at any phase and
  pathology, prints a pass count per section, and tags eleven checks
  KNOWN-FAIL: open defects that are printed but do not break the exit code.
  `tools/verify-timing.mjs` covers heart-rate warping and valve timing.
- **Pericardium and liver:**
  - The pericardial layer is measured in true distance, so no thick slab
    appears where two parts of the heart meet.
  - The left liver lobe is a curved wedge with a gastric impression.
  - The diaphragm domes rise laterally.
  - The subcostal window is closer.
- **Caval and RV inflow:** hepatic-vein confluence in the subcostal IVC view,
  and the CS and IVC orifices in the RV inflow view.
- **Motion:**
  - The mitral annulus displaces by angle, so the aorto-mitral curtain no
    longer lengthens.
  - The epicardial apex stays quasi-stationary, and the papillary muscles
    thicken.
  - Systolic intervals scale physiologically with heart rate, with diastole
    absorbing most of the change.
  - Right-heart valve timing is offset from the left.
- **Pathology:**
  - The effusion no longer surrounds the coronary sinus.
  - The VSD is labelled restrictive and the ASD secundum.
  - RWMA territories follow the coronary distribution.
  - One mitral valve area per grade drives both geometry and Doppler.
  - CW Doppler reads the modelled peak velocity.
  - A thin membranous septum marks the perimembranous VSD site.
- **TEE:**
  - The transgastric view sits on a gastric wall.
  - The descending-aorta view has lung behind the aorta and no LA in the
    sector.
  - The ME2C plane is rolled to include the LAA ostium.

Known open (tagged in the audit):
- The SVC reaches the RA only through the sinus venarum, and leaves the ME
  bicaval sector.
- The pericardium can touch LA blood at some phases in A4C, A2C and A3C.
- In DCM at ES, the tricuspid annulus still drops the normal amount while the
  mitral annulus barely moves.
- The ME4C RV is 0.65 cm short of the LV apex at ES (target 0.8).
- A marginal straight liver-edge stretch in PSAX-MV at ES.
- The LSPV is not a separate tube beside the LAA in the ME LAA view.
- No left brachiocephalic vein, and no Qp:Qs shunt model.
- Speckle and PSF stay horizontal rather than following the sector arcs.

## Known limitations

- It is one idealised adult heart. There is no patient-to-patient variation,
  no sex, body-size or age scaling, and no congenital variants beyond the modelled
  pathologies.
- Speckle follows a simplified deformation (axial shortening and radial
  scaling), not true myocardial strain; out-of-plane motion decorrelation is
  not modelled.
- LAA morphology is simplified and not seen in PSAX-AV.
- The lung acoustic windows are fixed cones and do not respond to probe
  position or respiration.
- Pathology geometry (DCM, AS/LVH, MR, MS, PH, effusion, ASD, VSD, RWMA, TR)
  comes from parameters. It is not image-derived.
- External expert validation (VALIDATION.md §B) remains outstanding.
