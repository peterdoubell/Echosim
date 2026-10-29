# Regional wall motion and DCM subcases

The regional selector changes the wall-motion territory, its regional-strain
plot, and the ECG injury vector together. Each ECG is projected from the same
three-dimensional dipole as the rhythm strip and 12-lead display, preserving
the limb-lead identities and producing reciprocal changes.

| Wall-motion example | Principal ECG lead group | Suggested echo view |
| --- | --- | --- |
| Septal | V1–V2 | A4C |
| Anteroseptal | V1–V4 | PSAX |
| Anterior | V3–V4 | A2C |
| Lateral | I, aVL, V5–V6 | A4C |
| Inferior | II, III, aVF | A2C |
| Inferolateral | II, III, aVF with I, V5–V6; reciprocal V1–V2 depression | PLAX/PSAX |

These are illustrative **acute infarction-associated** patterns. A wall-motion
abnormality does not by itself imply an acute infarction or ST elevation; it can
also reflect old scar, ischaemia without ST elevation, conduction abnormalities,
or other disease. Coronary territories overlap and vary between individuals.
The simplified dipole is not a diagnostic ECG model, and neighbouring leads can
show changes as well as the principal group. The existing LAD wall deformation
includes apical involvement; the septal label is not an isolated basal-septal
lesion. See the [AHA/ACCF/HRS ECG standardization statement, part VI](https://www.jacc.org/doi/10.1016/j.jacc.2008.12.016)
and the [Fourth Universal Definition of Myocardial Infarction](https://www.jacc.org/doi/10.1016/j.jacc.2018.08.1038).

DCM uses the project's existing **BSE 2020** EF terminology:

| Selection | BSE LVEF band | Solved example LVEF |
| --- | --- | --- |
| Borderline low | 50–54% | 52.3% |
| Impaired | 36–49% | 41.9% |
| Severely impaired | ≤35% | 26.8% |

Contractility, ventricular remodelling, filling pressure and secondary MR vary
between these examples. The circulation still solves EDV and ESV, and EF is
calculated from them; the selector does not overwrite a displayed EF number.
Both the rendered ventricular shortening and the strain plot use that solved
model. These are illustrative degrees of LV dysfunction in a dilated ventricle,
not a statement that EF alone establishes DCM. The default retains the original
severe case used in the learning pathway. The reference intervals and wording
come from the [British Society of Echocardiography guideline](https://pmc.ncbi.nlm.nih.gov/articles/PMC7040881/).

Run `node tools/verify-scenarios.mjs` to verify lead groups, reciprocal changes,
regional wall deformation, immediate case switching, EF bands, geometry/EF
agreement and model-cache separation between DCM severities.
