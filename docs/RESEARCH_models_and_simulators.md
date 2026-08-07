# Deep Research — Mathematical Heart Models & Echo Simulators (with ECG)

*Prepared to re-baseline EchoSim against the genuine state of the art. The
headline finding: measured against physics-based computational cardiology and
validated commercial/research echo simulators, EchoSim's previous "9.9/10"
scores were graded on a **pattern-recognition teaching-tool** bar. On a
**fidelity-vs-state-of-the-art** bar the honest figure is ~6.6–7.0. This
document is the evidence base; `BENCHMARKS.md` is the resulting rubric.*

---

## 1. Mathematical / computational models of the heart

### 1.1 Electrophysiology (the electrical substrate of the ECG)
The gold standard is a **reaction–diffusion PDE** over an anatomical mesh:
- **Cell (ionic) models** — human ventricular myocyte models such as
  **ten Tusscher–Panfilov (2006)** and **O'Hara–Rudy (ORd, 2011)** integrate
  dozens of ion-channel ODEs per cell; atrial (Courtemanche, Maleckar) and
  nodal models exist too. Minimal models (Bueno–Orovio, Aliev–Panfilov,
  FitzHugh–Nagumo) trade detail for speed.
- **Tissue propagation** — the **bidomain** equations are the most complete
  description (separate intra/extra-cellular potentials); the **monodomain**
  simplification is standard when conductivities are proportional, for large
  compute savings. The **eikonal** model approximates just the activation
  wavefront for real-time use.
- **Solvers** — **openCARP**, MonoAlg3D, Chaste, CARPentry run these on
  image-derived meshes with **fiber/sheet anisotropy**.

Sources: openCARP simulation environment (*Comput. Methods Programs Biomed.*,
2021); monodomain-vs-bidomain comparison; ten Tusscher–Panfilov & O'Hara–Rudy
ionic models; MonoAlg3D vs openCARP benchmark (CinC 2025).

### 1.2 Forward ECG modeling (the 12-lead the user asked about)
A real simulated ECG is **derived from** the electrical activation, not drawn:
1. Compute the 3D activation sequence (monodomain/bidomain/eikonal) with
   fiber-anisotropic conduction.
2. Map cardiac current sources to body-surface potentials via the
   **lead-field / bidomain torso** method (heart + lungs + torso as volume
   conductors), or the cheaper **pseudo-lead-field / pseudo-ECG**.
3. Read off the **standard 12 leads at ~1 kHz**.

State of the art simulates 3D ventricular activation at 1 mm and the full
12-lead ECG in **<1 s on GPU**, reproducing pathological morphologies (bundle
branch block, ischemia, ectopy). Large synthetic datasets exist
(**MedalCare-XL**: 16,900 healthy + pathological 12-lead ECGs from EP
simulation; **ECGSIM** for interactive source→ECG teaching).

Sources: eikonal + lead-field 12-lead (*Europace*, 2017); simplified 3D
whole-heart 12-lead model (PMC3654639); MedalCare-XL (arXiv 2211.15997);
ECGSIM.

> **EchoSim gap:** the ECG is a hand-authored analytic P-QRS-T bump with a
> phase cursor — a single rhythm strip, not derived from any electrical
> propagation, no 12-lead, no arrhythmias, no rate/conduction-dependent
> morphology. This is the single weakest axis.

### 1.3 Myocardial mechanics & the cardiac cycle
- **Passive tissue** — nonlinear, anisotropic hyperelastic constitutive laws:
  **Holzapfel–Ogden (2009)** (fiber/sheet/normal invariants) and **Guccione**
  (transversely isotropic). Solved by **FEM**.
- **Active contraction** — active-tension or active-strain models driven by
  intracellular Ca²⁺ (myofilament models), transmurally varying.
- **Circulation coupling** — chamber pressure/volume from **lumped-parameter
  (Windkessel) / CircAdapt** closed-loop models → physiological **pressure–
  volume loops**, ejection fraction, wall stress, **torsion** and regional
  strain.
- **Integrated simulators** — Dassault's **Living Heart Project** couples
  monodomain EP + Holzapfel–Ogden mechanics + circulation in one FE model.

Sources: Living Heart Project (*Eur. J. Mech. A/Solids*, 2014; PMC4175454);
Holzapfel–Ogden constitutive law; Guccione active tension; PV-loop coupling via
lumped circulation.

> **EchoSim gap:** contraction is **kinematic** — chamber ellipsoids/SDF
> scaled by a phase function. No stress/strain, no PV loop, no torsion, no
> regional wall-motion abnormality beyond a global EF scalar.

### 1.4 Intracardiac blood flow & valves (the Doppler substrate)
Gold standard is **computational fluid dynamics with fluid–structure
interaction (FSI)**: Navier–Stokes blood coupled to deforming myocardium and
**valve leaflets** (immersed-boundary or ALE-FEM), producing spatially and
temporally resolved velocity fields, **vortex rings**, and jets that conserve
mass/momentum. Validation against clinical Doppler shows transvalvular peak
velocities within **~8%**.

Sources: mitral/aortic FSI in a realistic LV (PMC5590990); FSI validated vs
Doppler echo, 7.9–8.4% under-prediction (*Comput. Biol. Med.*, 2024);
whole-heart image-based flow (arXiv 2605.09629).

> **EchoSim gap:** flow is a **prescribed piecewise-Gaussian bulk-velocity
> field** (kinematic jets/PISA), not a conserved solution of the flow
> equations — no vortices, no true continuity, no pressure field.

### 1.5 Anatomy & fiber architecture
Research/commercial models use **patient-specific or statistical-shape**
whole-heart meshes segmented from CT/MRI, with **DT-MRI–measured** or
**rule-based (Streeter)** myofiber orientation, trabeculation, papillary
muscles, chordae, and coronary trees (Zygote Solid Heart, Living Heart,
openCARP atlases).

> **EchoSim gap:** one idealized blended-SDF geometry, no fibers, no chordae,
> stylized trabeculation/papillaries, single generic (non-patient) heart.

---

## 2. Echocardiography simulators (commercial + research)

### 2.1 Commercial
- **CAE / Elevate Healthcare VIMEDIX** — mannequin + tracked transducer;
  **real-time volumetric (3D/4D) scanning** with multiplanar reconstruction for
  **TTE, TEE and abdominal**; normal + advanced pathologies; synchronized
  animated 3D anatomy; peer-reviewed realism/learning studies.
- **HeartWorks (Inventive Medical / Intelligent Ultrasound)** — built by UCLH
  cardiac anaesthetists; **anatomically accurate interactive 3D heart**;
  realistic 2D **and** 3D echo with simultaneous 3D model; **PW, CW and colour
  Doppler**; **~30 pathology cases**; TTE/TEE/lung packages.
- Others: 3D Systems / Simbionix **U/S Mentor**, MedaPhor/Intelligent Ultrasound
  **ScanTrainer**, **Echocom**.

These run on a **haptic mannequin with probe tracking**, include **TEE and 3D/
MPR**, and are validated in the education literature. They pair the image with a
scripted/animated **ECG and haemodynamic** context.

Sources: CAE Vimedix product & studies (caehealthcare.com; cae.com press
releases); HeartWorks (inventivemedical.com; academy.intelligentultrasound.com;
morethansimulators.com).

### 2.2 Research-grade image simulation
- **Field II** (Jensen) — spatial-impulse-response scatterer simulation, the de
  facto reference; **FOCUS**, **k-Wave** (full-wave) similar class.
- **COLE / SIMUS** — fast **convolution** of a (Field II-derived) point-spread
  function with a scatterer map; COLE adds realistic anisotropic **beam
  profiles**; used to build synthetic cardiac datasets.
- **DT-MRI-driven** cardiac speckle simulation; **GAN / diffusion** patho-
  realistic synthesis; **ray-based** scattering (UltraScatter, 2025).
- Benchmarks/datasets: **CAMUS** (segmented 2-/4-chamber), autonomous-navigation
  slice simulators from CT meshes.
- Quantitative image-quality metrics: **gCNR**, CNR, speckle-SNR, resolution —
  used to validate that simulated speckle statistics match real B-mode.

Sources: Field II vs FOCUS; COLE convolution methodology (Springer, MICCAI);
SIMUS (Garcia); DTI-based cardiac simulation (PMC4537486); GAN patho-realistic
(arXiv 1712.07881); UltraScatter (arXiv 2510.10612); autonomous-nav simulation
(Frontiers, 2024); gCNR & image-quality assessment (*J. Med. Ultrason.*, 2021;
arXiv 2408.00591).

> **EchoSim gap:** the B-mode is an **analytic tissue-echogenicity + procedural
> correlated-speckle** render with heuristic TGC and BART colour — no wave
> physics, no scatterer/PSF model, speckle statistics not validated against
> Rayleigh/K-distributions, and few artifacts (reverberation, shadowing, side
> lobes, blooming) beyond a bright pericardial line.

---

## 3. Where EchoSim genuinely stands (honest, evidence-anchored)

| Capability | State of the art | EchoSim | Honest band |
|---|---|---|---|
| Anatomy & fibers | Patient-specific CT/MRI mesh + DT-MRI fibers | Idealized SDF, no fibers | 6.5 |
| Electrophysiology + ECG | Ionic→monodomain→lead-field **12-lead** | Hand-drawn 1-lead PQRST | 4.0 |
| Myocardial mechanics | FEM Holzapfel–Ogden + active tension + PV loop | Kinematic scaling | 5.5 |
| Blood flow / Doppler | CFD-FSI, valves, validated ±8% vs Doppler | Prescribed Gaussian jets | 6.0 |
| Ultrasound image formation | Field II / COLE scatterer+PSF, gCNR-validated | Analytic echo + speckle | 6.5 |
| Pathology breadth & depth | 30+ cases, TEE/3D, severity physics | 9 stylised cases, TTE-only | 6.5 |
| External validation | Peer-reviewed Likert/RCT, CAMUS/MedalCare | LLM-board only | 3.5 |
| Real-time interactivity/UX | Haptic mannequin + probe tracking + 3D/MPR | Browser, moveable plane, strong UX | 8.0 |
| Education scaffolding | Curriculum, logbook | Labels, quiz, measurements, BSE ranges | 8.0 |
| Engineering/code | Validated research codebases | Clean self-contained web app | 8.0 |

A physics-and-clinical-weighted mean lands at **~6.6–7.0** — consistent with the
user's recalibration. EchoSim is an excellent *real-time schematic trainer*; it
is **not** a physics-based digital twin, and the ECG in particular is
illustrative rather than simulated.

*(All claims above are drawn from the cited literature; URLs collected in the
commit message and BENCHMARKS.md references section.)*
