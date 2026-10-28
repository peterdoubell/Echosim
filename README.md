# EchoSim

An interactive 3D **in-silico echocardiography trainer** that runs entirely in the
browser. A signed-distance-field heart, a lumped-parameter closed-loop
circulation, a physics-derived 12-lead ECG and a scatterer/PSF B-mode renderer
all share one model, so the 3D anatomy view and the 2D echo are always
consistent with each other.

> **Education only.** This is a teaching model, not a medical device, and not for
> diagnosis. Measurements are illustrative.

## Running it

It is a static site with no build step and no runtime dependencies — three.js is
vendored. Serve the directory with anything:

```bash
python3 -m http.server 8000
# then open http://localhost:8000
```

ES modules require a real HTTP origin, so opening `index.html` from the
filesystem will not work.

## What it models

**Anatomy** — a continuous myocardium with a crescentic RV wrapping a shared
septum; an explicit interatrial septum thin over the fossa ovalis; a
narrow-ostium multi-lobed left atrial appendage against the broad-based
triangular right; crista terminalis, Eustachian valve, SVC/IVC, four pulmonary
veins and the LSPV/LAA ridge; papillary muscles, moderator band and chordae.

**Valves** — leaflets hinge at the annulus and swing as curved, tapered-thin
curtains. The mitral is bileaflet-asymmetric with a scalloped posterior leaflet;
the semilunar valves are individuated into three cusps with a trilobed orifice
and commissural coaptation seams (the short-axis "Mercedes" Y).

**Physiology** — a closed-loop circulation (time-varying elastance + Windkessel)
drives a real LV pressure–volume loop, and the cavity geometry is volume-exact
against it, so measured EF and LVIDd fall out of the modelled volume rather than
a prescribed curve. The atria run anti-phase to the ventricles on an explicit
reservoir/conduit/booster curve — peak volume at AV-valve opening, emptying on
the E wave and again on the A wave. Continuity-equation and Bernoulli Doppler, AHA 17-segment
regional strain, MAPSE, mitral E/A, and pulmonary-vein S/D waves with systolic
flow reversal in severe MR.

**Imaging** — a round-trip beam (a once-focused transmit limb multiplied by a
dynamically-refocused receive limb), tissue harmonic imaging with its honest
penetration penalty, finite elevational slice thickness giving genuine
partial-volume averaging, persistence, and a frame rate reported as measured.

**Reporting** — a structured echo report for every case (📝 in the
Measurements panel): measurements taken from the model, each graded against the
sex-specific BSE 2020 reference intervals and indexed to the patient's BSA or
height; findings structure by structure; a numbered conclusion; and fill-in
templates for writing your own report and comparing. The full guide, reference
tables and templates are in [`docs/REPORTING.md`](docs/REPORTING.md).

## Verifying it

```bash
npm install playwright          # browsers are expected to be preinstalled
node tools/validate.mjs         # 12 ECG morphology checks
node tools/shoot.mjs out/       # screenshots across views/pathologies, asserts 0 console errors
node tools/verify-measure.mjs   # end-to-end caliper + Simpson-EF check
node tools/verify-atrial-phase.mjs  # atrio-ventricular phase relationship
node tools/verify-anatomy.mjs   # anatomy audit: ASE dimensions and relationships, swept
                                # over views, phases and pathologies (open defects tagged KNOWN-FAIL)
node tools/verify-timing.mjs    # valve timing, heart-rate warping, right-heart offsets
node tools/verify-report.mjs    # structured echo report: BSE 2020 grading, conclusions, templates, dialog
node tools/verify-feedback.mjs  # feedback dialog and api/feedback.js (validation, rate limit, fallbacks)
```

`docs/ANATOMY.md` describes the anatomical model, the audit and the expert-panel
review history.

## Feedback (bug reports and feature requests)

The **💬 Feedback** button in the header opens a dialog with *Report a bug* and
*Request a feature* tabs. Messages are emailed to the maintainer by the Vercel
function `api/feedback.js`. Bug reports include diagnostics by default: the
case, view, phase, heart rate, browser, screen, WebGL renderer and recent
JavaScript errors. The user can see these and switch them off.

Set these in Vercel → Project → Settings → Environment Variables, then redeploy:

| Variable | Purpose |
|---|---|
| `FEEDBACK_TO` | Recipient address. It is kept out of the public code. |
| `RESEND_API_KEY` | API key from [Resend](https://resend.com) (free tier). Without it, the dialog opens the user's own mail app with the message pre-filled, addressed to `FEEDBACK_TO`. |
| `FEEDBACK_FROM` | Optional sender, e.g. `EchoSim <feedback@your-domain>`. By default Resend's shared sender is used, which only delivers to the address the Resend account was created with. |

If neither variable is set, or the app runs on a static local server, the
message is copied to the user's clipboard so nothing is lost. The function
limits each IP to 5 messages per 10 minutes and silently drops bots that fill
a hidden honeypot field.

## Honest limits

`BENCHMARKS.md` scores the simulator against a weighted rubric, and
`docs/ROADMAP_TO_9.9.md` is explicit about what code can and cannot close. The
remaining gap to any *validated clinical-training* claim is **external validation
and regulatory/QMS evidence** — an expert panel, construct and transfer validity,
a risk file — which is funded human work, not something a commit can deliver.
Anatomical face validity is still pending expert sign-off; see
`docs/VALIDATION_RESULTS.md`.

Development before the move to this repository happened in
`peterdoubell/desktop-tutorial`; its history is in pull requests
[#2](https://github.com/peterdoubell/desktop-tutorial/pull/2) and
[#69](https://github.com/peterdoubell/desktop-tutorial/pull/69) there. This
repository is now the only home of EchoSim.
