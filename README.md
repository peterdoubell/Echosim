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

## Verifying it

```bash
npm install playwright          # browsers are expected to be preinstalled
node tools/validate.mjs         # 12 ECG morphology checks
node tools/shoot.mjs out/       # screenshots across views/pathologies, asserts 0 console errors
node tools/verify-measure.mjs   # end-to-end caliper + Simpson-EF check
node tools/verify-atrial-phase.mjs  # atrio-ventricular phase relationship
```

## Honest limits

`BENCHMARKS.md` scores the simulator against a weighted rubric, and
`docs/ROADMAP_TO_9.9.md` is explicit about what code can and cannot close. The
remaining gap to any *validated clinical-training* claim is **external validation
and regulatory/QMS evidence** — an expert panel, construct and transfer validity,
a risk file — which is funded human work, not something a commit can deliver.
Anatomical face validity is still pending expert sign-off; see
`docs/VALIDATION_RESULTS.md`.

Development history for the work up to this point lives in the pull request on
the repository this was extracted from.
