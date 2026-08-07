// echo.js
// Renders the 2D ultrasound "machine" view by sampling the analytic cardiac
// model along the imaging plane. Produces a fan-shaped B-mode sector with
// speckle, depth attenuation + time-gain compensation and log compression, a
// translucent colour-Doppler overlay obeying the "Blue Away, Red Towards"
// (BART) convention with a wall filter and aliasing, on-image anatomical
// labels, live measurements, and a scrolling spectral-Doppler trace.

import { classify, velocityAt, geometryAt, TISSUE, VALVE_DEFS, FLOW, AORTIC_CUSPS } from './cardiac-model.js';
import { clamp } from './mathutils.js';

// Human-readable labels for the flow compartment that produced the frame's peak
// velocity (surfaced to the UI as metrics.peakLabel for teaching context).
const FLOW_LABELS = {
  [FLOW.MITRAL_IN]: 'LV inflow (E)',
  [FLOW.TRICUSPID_IN]: 'RV inflow',
  [FLOW.LVOT]: 'LVOT',
  [FLOW.RVOT]: 'RVOT',
  [FLOW.MR_JET]: 'MR jet',
  [FLOW.AS_JET]: 'AS jet',
  [FLOW.VSD_JET]: 'VSD jet',
  [FLOW.ASD_JET]: 'ASD jet',
  [FLOW.TR_JET]: 'TR jet',
  [FLOW.PV_FLOW]: 'pulmonary vein',
};

// Structures expected to lie in each standard imaging plane. When the UI sets
// this.viewName, _labels() only annotates structures on this list (plus the
// tightened geometric plane test); when unset it falls back to geometry alone.
const VIEW_WHITELIST = {
  PLAX: ['LV', 'LA', 'Ao', 'AV', 'MV', 'RV'],
  PSAX: ['LV', 'RV'],
  PSAX_AV: ['RCC', 'LCC', 'NCC', 'RV', 'RA', 'LA'],
  MELAA: ['LAA', 'LA', 'LV'],
  A4C: ['LV', 'RV', 'LA', 'RA', 'MV', 'TV'],
  A2C: ['LV', 'LA', 'MV'],
  SUBCOSTAL: ['LV', 'RV', 'LA', 'RA'],
};

// sampling resolution of the offscreen buffer (upscaled to the visible canvas).
// Raised from 240x300 toward modality-grade spatial resolution: a real adult TTE
// cine is ~0.3-0.5 mm axial / 1-2 mm lateral. At depth 15 cm this buffer gives
// ~0.36 mm/px axially, so endocardial borders and speckle grain resolve finely on
// capable hardware; the FPS-adaptive stride (up to 3) coarsens the SAMPLING on slow
// machines while the larger buffer still upscales less blockily than before.
const SW = 336;
const SH = 420;
const WALL_FILTER = 0.06; // m/s colour dead-zone (suppresses near-zero flow)

// Elevational (out-of-plane) beam geometry. A cardiac phased array focuses the
// elevation plane with a FIXED acoustic lens, so the slab is thinnest at a set
// depth (~8 cm) and thickens both nearer and deeper: roughly 2.5 mm full
// thickness at the elevation focus out to ~9 mm in the near field, matching
// published slice-thickness measurements for adult cardiac probes.
// Receive-beam constants. RX_W is the (depth-independent) receive beam width the
// scanner maintains by growing the aperture at constant F-number; RX_APERTURE_CM
// is the depth at which the physical aperture is fully open, beyond which the
// F-number can no longer be held and the receive beam widens in proportion to depth.
const RX_W = 2.6;
const RX_APERTURE_CM = 8.0;

const ELEV_FOCUS = 8.0;   // cm — fixed by the lens, NOT steerable
const ELEV_MIN = 0.09;    // cm half-thickness at the elevation focus (~1.8 mm slab)
const ELEV_DIV = 0.025;   // cm half-thickness added per cm away from that focus
// Stratified 2x2 dither over the slab. Four offsets spanning the thickness, laid
// out so that ANY 2x2 neighbourhood sums to zero: the smallest averaging window
// the PSF can supply already cancels the local bias, which is what keeps the
// Monte-Carlo estimate from reading as mottle rather than as partial-volume
// softening. A larger (e.g. 4x4 Bayer) pattern samples the slab more finely but
// needs a wider kernel to cancel, and visibly blotches where the beam is tight.
const ELEV_DITHER = Float32Array.from([-1, 0.333, -0.333, 1]); // sums to zero

export class EchoView {
  constructor(canvas, spectralCanvas) {
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d');
    this.buf = document.createElement('canvas');
    this.buf.width = SW; this.buf.height = SH;
    this.bctx = this.buf.getContext('2d');
    this.img = this.bctx.createImageData(SW, SH);

    this.spectral = spectralCanvas;
    this.sctx = spectralCanvas ? spectralCanvas.getContext('2d') : null;
    this.spectralHistory = [];

    this.gain = 1.0;
    this.depthCm = 17;         // display depth
    this.sectorHalf = 0.66;    // ~38 degrees each side
    this.nyquist = 0.62;       // m/s colour scale limit
    this.colorOn = true;
    this.showColorBox = true;
    this.showLabels = true;    // on-image anatomical labels (toggled by UI)
    this.viewName = undefined; // optional standard-view name set by the UI (PLAX/PSAX/A4C/A2C/SUBCOSTAL)

    // live measurements consumed by the UI each frame. gcnr/cnr/speckleSNR are a
    // live image-quality self-validation readout computed from the rendered
    // B-mode each frame (see _computeImageQuality).
    this.metrics = { lvidd: null, lvids: null, ef: null, peakVel: null, peakGrad: null, effusion: null, peakLabel: null, severity: null, severityLabel: null,
      gcnr: null, cnr: null, speckleSNR: null, pasp: null, mapse: null, fps: null };
    this._lastFrameT = 0;      // wall-clock of the previous rendered frame
    this._emaFrameMs = 0;      // smoothed frame interval -> displayed frame rate

    // transducer centre frequency (MHz). Real adult TTE phased arrays run
    // ~2-3.5 MHz; the PSF widths scale with 1/freq (higher freq => tighter beam,
    // finer speckle). The probe object carries no frequency, so it lives here.
    this.freqMHz = 2.7;
    // Tissue harmonic imaging: transmit at f, receive at 2f. Off by default so the
    // fundamental image (and the gCNR/speckle figures validated against it) is
    // unchanged; turning it on is a genuine, measurable contrast improvement.
    this.harmonic = false;
    // Transmit focal depth (cm). null = auto, tracking 55% of the display depth
    // the way a scanner's default focus sits in the middle of the field. Setting
    // it explicitly is a real skill: put the focus on what you are measuring.
    this.focusCm = null;
    // Persistence (frame averaging), 0..0.85. Scanners blend each new frame with
    // the previous ones to suppress speckle noise; the cost is temporal blurring,
    // so high persistence smears fast-moving structures — valves especially. It
    // is a genuine trade-off control, not a quality slider.
    this.persistence = 0;
    // Elevational slice thickness (out-of-plane partial-volume averaging). On by
    // default because every real probe has it; exposed so it can be switched off
    // to show learners exactly what artefact it is responsible for.
    this.elevation = true;
    this._prevGrey = null;     // lazily allocated when persistence is first used

    // running peak jet velocity since the pathology last changed — a stable
    // source for velocity-based severity grading (see _computeMetrics).
    this._peakVelMax = 0;

    // steered spectral gate: world point + line-of-sight velocity of this frame's
    // peak flow (null when no significant flow — _spectral falls back to centre).
    this._peakGate = null;

    // memoised fixed-phase geometry for _computeMetrics (recomputed only when the
    // path/pathology object identity changes).
    this._lastMetricsPath = null;
    this._ged = null;
    this._ges = null;

    // reusable scratch for beam direction (avoid per-pixel allocation)
    this._bd = [0, 0, 0];

    // FPS-adaptive sampling: sample the sector on a coarser grid (stride) and
    // block-fill when a frame's sampling cost exceeds the budget, restoring full
    // resolution when there's headroom. Keeps the app smooth on slow hardware.
    this._stride = 1;
    this._emaMs = 6;

    // deterministic speckle texture
    this.noise = new Float32Array(SW * SH);
    let s = 12345;
    for (let i = 0; i < this.noise.length; i++) {
      s = (s * 1103515245 + 12345) & 0x7fffffff;
      this.noise[i] = (s / 0x7fffffff);
    }
    // Pre-correlate the speckle so it has a realistic grain instead of 1-pixel
    // uncorrelated white noise: a small separable 3-tap blur (correlation length
    // ~2px) followed by a contrast stretch that restores the variance the blur
    // removed. Deterministic (operates on the seeded buffer above).
    this._correlateNoise();

    // --- scatterer + PSF image-formation scratch (COLE-style) -----------------
    // The B-mode is formed as tissue reflectivity (with sub-resolution scatterer
    // speckle) convolved with a depth/frequency-dependent separable PSF, then
    // display-processed (TGC, artifacts, log compression). These reusable buffers
    // hold the intermediate reflectivity field and the per-pixel context needed
    // for the post-passes, so nothing is allocated per frame.
    this._refl = new Float32Array(SW * SH);   // linear reflectivity (pre-PSF, then blurred in place)
    this._tmp = new Float32Array(SW * SH);    // scratch for the separable blur
    this._depth = new Float32Array(SW * SH);  // cm depth per pixel (0 outside fan)
    this._abin = new Uint16Array(SW * SH);    // artifact beam bin per pixel
    this._kind = new Uint8Array(SW * SH);     // tissue kind per pixel (TISSUE.*)
    this._mask = new Uint8Array(SW * SH);     // 1 = inside sector, 0 = background
    this._dcol = new Float32Array(SW * SH * 3); // Doppler colour (kept crisp, composited after PSF)
    this._dalpha = new Float32Array(SW * SH);   // Doppler overlay alpha

    // acoustic-artifact fields marched in polar (angle x depth) space each frame.
    this._NA = 128;                           // artifact beam count across the sector
    this._NS = 96;                            // depth samples per beam
    this._shadowGain = new Float32Array(this._NA * this._NS); // distal attenuation (shadowing)
    this._enhGain = new Float32Array(this._NA * this._NS);    // posterior enhancement

    // Precomputed normalised 1-D Gaussian PSF kernels indexed by integer radius
    // 0..RMAX. The blur picks a kernel per pixel from its depth/frequency, so no
    // exp() runs in the hot loop. r=0 is the identity (near-field crispness).
    this._RMAX = 8;
    this._kern = [];
    for (let r = 0; r <= this._RMAX; r++) {
      const k = new Float32Array(2 * r + 1);
      const sig = Math.max(0.5, r * 0.62);
      let sum = 0;
      for (let t = -r; t <= r; t++) { const wv = Math.exp(-(t * t) / (2 * sig * sig)); k[t + r] = wv; sum += wv; }
      for (let t = 0; t < k.length; t++) k[t] /= sum;
      this._kern.push(k);
    }
  }

  // Lateral (across-beam) PSF radius in px for a given depth. Lateral resolution
  // is worse than axial (wider beam) and both degrade with depth (beam
  // divergence) and improve with frequency (freqScale = ref/freq, <1 = tighter).
  // Coefficients are scaled to the SW=336/SH=420 buffer so the PHYSICAL beam
  // width (mm) matches the prior validated image (gCNR/speckle unchanged); the
  // sharpness gain comes from the denser sampling, not a tighter point spread.
  // A focused beam is NARROWEST AT ITS FOCUS and diverges on both sides of it —
  // it does not simply widen with depth from the transducer face. Modelling the
  // focus properly is what makes "set your focus at the structure of interest"
  // a real, learnable skill here rather than a decorative control. The residual
  // `depth` term carries the aperture/attenuation penalty that still degrades the
  // far field even at focus.
  _focusCm() {
    return this.focusCm == null ? this.depthCm * 0.55 : this.focusCm;
  }

  // TWO-WAY beam. The round-trip point spread is the PRODUCT of the transmit and
  // receive beam patterns, and the two behave completely differently:
  //
  //  - Transmit is focused ONCE, at a depth the operator picks, so it is narrow
  //    at that focus and diverges either side of it.
  //  - Receive is DYNAMICALLY refocused as each echo returns, so it is in focus
  //    at every depth. Its width is set by the F-number, which the scanner holds
  //    constant by expanding the active aperture with depth — giving essentially
  //    depth-independent receive resolution until the physical aperture runs out,
  //    after which the F-number (and the beam) grows linearly with depth.
  //
  // Combining them as Gaussians (1/w² = 1/w_tx² + 1/w_rx²) means the round-trip
  // beam is never wider than the narrower of the two, so the in-focus receive
  // beam limits how badly the image degrades away from the transmit focus. That
  // is why real scanners stay usable off-focus, and it is the physics a
  // transmit-only model gets wrong: it over-punishes the near and far field.
  _latRadius(depth, freqScale, focus) {
    const wTx = 1.4 + Math.abs(depth - focus) * 0.42;            // focused once
    const wRx = depth <= RX_APERTURE_CM                          // dynamically focused
      ? RX_W : RX_W * (depth / RX_APERTURE_CM);                  // aperture-limited beyond
    const w = (wTx * wRx) / Math.sqrt(wTx * wTx + wRx * wRx);    // round trip
    return clamp(Math.round(w * freqScale), 0, this._RMAX);
  }

  // Axial (along-beam) PSF radius in px — always tighter than lateral, and NOT
  // focus-dependent at all: axial resolution is set by the pulse length, not by
  // beam geometry. It still degrades with depth because attenuation preferentially
  // removes the high-frequency content, downshifting the centre frequency and so
  // stretching the effective pulse.
  _axRadius(depth, freqScale) {
    return clamp(Math.round((0.25 + depth * 0.17) * freqScale), 0, this._RMAX);
  }

  // Tissue-harmonic beam narrowing. Harmonic energy is generated in tissue in
  // proportion to the SQUARE of the local fundamental pressure, so it only forms
  // in the intense core of the beam: the effective harmonic beam is substantially
  // narrower than the fundamental one, and the effective pulse is shorter. The
  // lateral gain is the larger of the two, which is why harmonic imaging cleans up
  // side-lobe smear across a cavity so much more than it sharpens along the beam.
  // Returns [lateralFactor, axialFactor] applied on top of freqScale.
  _harmonicScale() {
    return this.harmonic ? [0.62, 0.80] : [1, 1];
  }

  // Separable [0.25,0.5,0.25] blur of the noise buffer (treated as SW x SH),
  // then a stretch about the 0.5 mean to restore speckle contrast. Run once.
  _correlateNoise() {
    const n = this.noise, w = SW, h = SH;
    const tmp = new Float32Array(n.length);
    for (let y = 0; y < h; y++) {           // horizontal pass
      const row = y * w;
      for (let x = 0; x < w; x++) {
        const a = n[row + (x > 0 ? x - 1 : x)];
        const b = n[row + x];
        const c = n[row + (x < w - 1 ? x + 1 : x)];
        tmp[row + x] = 0.25 * a + 0.5 * b + 0.25 * c;
      }
    }
    for (let x = 0; x < w; x++) {            // vertical pass + contrast restore
      for (let y = 0; y < h; y++) {
        const a = tmp[(y > 0 ? y - 1 : y) * w + x];
        const b = tmp[y * w + x];
        const c = tmp[(y < h - 1 ? y + 1 : y) * w + x];
        const blur = 0.25 * a + 0.5 * b + 0.25 * c;
        n[y * w + x] = clamp(0.5 + (blur - 0.5) * 2.66, 0, 1); // ~restore std
      }
    }
  }

  // beam direction in heart space for a given fan angle (theta from centre),
  // written into the reusable scratch array this._bd.
  beamDir(probe, theta) {
    const c = Math.cos(theta), s = Math.sin(theta);
    this._bd[0] = c * probe.dir[0] + s * probe.lat[0];
    this._bd[1] = c * probe.dir[1] + s * probe.lat[1];
    this._bd[2] = c * probe.dir[2] + s * probe.lat[2];
    return this._bd;
  }

  render(G, path, probe, shimmer) {
    const data = this.img.data;
    const cx = SW * 0.5;
    const apexY = SH * 0.04;
    const pxPerCm = (SH * 0.92) / this.depthCm;
    const half = this.sectorHalf;
    const P = probe.pos;
    // subtle live speckle shimmer: rotate the noise lookup slowly over time
    const shim = ((shimmer || 0) * 7) | 0;
    let peakSpeed = 0, peakFlow = FLOW.NONE;
    // world point + line-of-sight velocity of the peak-flow sample this frame,
    // used to steer the spectral (PW) gate onto the actual jet.
    let pkx = 0, pky = 0, pkz = 0, pkLos = 0, pkTurb = false, hasPeak = false;
    const t0 = (typeof performance !== 'undefined' ? performance.now() : 0);
    const st = this._stride;
    const refl = this._refl, depthBuf = this._depth, abin = this._abin;
    const kindBuf = this._kind, maskBuf = this._mask, dcol = this._dcol, dalpha = this._dalpha;
    const NA = this._NA, halfInv = 1 / (2 * half);
    const elevOn = this.elevation !== false;
    // The stochastic slice-thickness estimate needs neighbouring samples to
    // average. When the sampler is coarsened on slow hardware each block carries a
    // SINGLE offset, so the reconstruction can no longer cancel and the effect
    // would read as blotching rather than partial-volume softening. Scale the
    // slab back with the stride: honest about what the sampling can support.
    const elevAmp = 1 / st;
    // hoist the plane normal out of the loop: repeated element reads off a plain
    // array are measurably slower than locals at this call frequency
    const eN = probe.normal, enx = eN[0], eny = eN[1], enz = eN[2];
    // Store the sampled scatterer/context of the st x st block anchored at (i,j)
    // into the reflectivity + context buffers. The B-mode grey is NOT formed here:
    // reflectivity is first PSF-convolved, then display-processed (TGC, artifacts,
    // log compression) and finally composited with the crisp Doppler overlay.
    const store = (i, j, e, depth, kind, mask, ab, dr, dg, db, da) => {
      for (let bj = j; bj < j + st && bj < SH; bj++) {
        for (let bi = i; bi < i + st && bi < SW; bi++) {
          const idx = bj * SW + bi;
          refl[idx] = e; depthBuf[idx] = depth; kindBuf[idx] = kind;
          maskBuf[idx] = mask; abin[idx] = ab;
          dalpha[idx] = da;
          if (da > 0) { const c3 = idx * 3; dcol[c3] = dr; dcol[c3 + 1] = dg; dcol[c3 + 2] = db; }
        }
      }
    };

    // bj/bi are SAMPLING-BLOCK indices, carried incrementally so the elevation
    // dither can be looked up without a per-sample division in the hot loop.
    for (let j = 0, bj = 0; j < SH; j += st, bj++) {
      for (let i = 0, bi = 0; i < SW; i += st, bi++) {
        const dx = i - cx;
        const dy = j - apexY;
        const rpix = Math.hypot(dx, dy);
        const theta = Math.atan2(dx, dy);
        if (Math.abs(theta) > half || dy <= 0 || rpix > SH * 0.94) {
          store(i, j, 0, 0, TISSUE.OUTSIDE, 0, 0, 0, 0, 0, 0);
          continue;
        }
        const depth = rpix / pxPerCm; // cm along the beam
        const bd = this.beamDir(probe, theta);
        let wx = P[0] + depth * bd[0];
        let wy = P[1] + depth * bd[1];
        let wz = P[2] + depth * bd[2];

        // --- elevational slice thickness (partial-volume averaging) ------------
        // A scan plane is not a mathematical surface: the beam has a finite
        // thickness OUT of plane, set by a fixed acoustic lens, so each pixel
        // reports the average over a slab, not a point. That is why real echo
        // shows structures fading in and out at plane edges and why thin
        // structures look softer than a geometric cut would suggest. The slab is
        // thinnest at the ELEVATION FOCUS — which, unlike the transmit focus, is
        // fixed by the lens and cannot be steered; that is a genuine constraint
        // operators work around by angling the probe, not a control they have.
        // Sampling is stochastic: each pixel takes ONE sample at a dithered
        // offset across the slab, and the PSF convolution downstream averages
        // neighbours that sampled different offsets — a Monte-Carlo estimate of
        // the partial-volume average at essentially no extra sampling cost.
        if (elevOn) {
          const hw = (ELEV_MIN + Math.abs(depth - ELEV_FOCUS) * ELEV_DIV) * elevAmp;
          // indexed by SAMPLING BLOCK, not pixel, so the pattern survives striding
          const off = hw * ELEV_DITHER[((bj & 1) << 1) | (bi & 1)];
          wx += enx * off; wy += eny * off; wz += enz * off;
        }

        const t = classify(wx, wy, wz, G, path);
        const tissue = t.echo; // capture before velocityAt reuses singletons
        const kind = t.tissue;

        // --- B-mode grey from echogenicity, speckle & depth gain ---
        // Correlated speckle sampled at a depth-dependent scale so the grain
        // coarsens slightly in the far field (deeper = larger speckle cells).
        const sc = 1 + depth * 0.06;
        // Depth-dependent offset (deterministic) so successive far-field bands
        // sample DIFFERENT regions of the noise buffer instead of magnifying and
        // repeating one enlarged top-left patch. The offsets wrap independently
        // in x and y so the grain spans the whole SW*SH buffer at every depth.
        const dOff = (depth * 11.0) | 0;
        const si = (((i / sc) | 0) + dOff) % SW;
        const sj = (((j / sc) | 0) + dOff * 3) % SH;
        const n = this.noise[((sj * SW + si) + shim) % this.noise.length];
        const nn = this.noise[((sj * SW + ((si + 37) % SW)) + shim) % this.noise.length];
        let e = tissue;
        if (kind === TISSUE.MYO) e = (0.16 + 0.20 * tissue) * (0.7 + 0.5 * n); // believable mid-grey, capped ~0.45, speckle kept
        else if (kind === TISSUE.VALVE) {
          // Normal leaflets are only moderately echogenic (~0.5-0.6); reserve the
          // bright, near-white band for calcified valves. The model flags calcified
          // leaflets with a higher raw echo (0.98 vs 0.85 normal) — key on that so
          // a healthy mitral/tricuspid plane no longer renders as a solid bright bar.
          const calcified = tissue > 0.92;
          e = calcified ? (0.72 + 0.30 * n) : (0.40 + 0.28 * n);
        }
        else if (kind === TISSUE.LV || kind === TISSUE.RV ||
                 kind === TISSUE.LA || kind === TISSUE.RA ||
                 kind === TISSUE.AORTA) e = 0.02 + 0.05 * n * n; // near-anechoic blood + faint noise floor
        else if (kind === TISSUE.PERICARDIUM) {
          e = 0.003; // echo-free effusion: distinctly dark, sharply bounded
          // inner (visceral) pericardial line where the effusion abuts myocardium:
          // probe one short step TOWARD the transducer for the myo interface.
          const nbIn = classify(wx - 0.2 * bd[0], wy - 0.2 * bd[1], wz - 0.2 * bd[2], G, path);
          if (nbIn.tissue === TISSUE.MYO) e = 0.85;
          else {
            // mirror: thin bright parietal-pericardium / lung line on the OUTER
            // side — probe one step AWAY from the transducer; a non-fluid,
            // non-myo neighbour there marks the fluid->outside interface.
            const nbOut = classify(wx + 0.2 * bd[0], wy + 0.2 * bd[1], wz + 0.2 * bd[2], G, path);
            if (nbOut.tissue !== TISSUE.PERICARDIUM && nbOut.tissue !== TISSUE.MYO) e = 0.8;
          }
        }
        else {
          // surrounding soft tissue / lung: a dim non-zero speckle floor so a
          // truly anechoic effusion reads as distinctly darker than the surround.
          e = Math.max(0.035 + 0.05 * nn, tissue * (0.35 + 0.8 * nn));
        }

        // reflectivity `e` is the scatterer field; TGC / artifacts / log
        // compression are applied AFTER the PSF convolution (see _psfBlur +
        // composite pass below). Only capture the crisp Doppler overlay here.
        let dR = 0, dG = 0, dB = 0, dA = 0;

        // --- colour Doppler overlay ---
        if (this.colorOn && depth < this.depthCm * 0.92) {
          const vel = velocityAt(wx, wy, wz, G, path);
          if (vel && vel.speed > 0.05) {
            // line-of-sight velocity: beamDir points away from probe
            const los = vel.vx * bd[0] + vel.vy * bd[1] + vel.vz * bd[2];
            const alos = Math.abs(los);
            // Steer the spectral gate to the TRUE peak-SPEED sample anywhere along
            // the beam — tracked independently of the wall filter / line-of-sight
            // projection so a jet running across the beam (small los, large speed,
            // e.g. severe MR) still fixes the gate on its real core, not on the
            // peak line-of-sight sample (which under-reads the true jet velocity).
            if (vel.speed > peakSpeed) {
              peakSpeed = vel.speed;
              peakFlow = vel.flow;
              pkx = wx; pky = wy; pkz = wz;
              pkLos = los; pkTurb = vel.turbulent; hasPeak = true;
            }
            // Depth-dependent EFFECTIVE Nyquist: PRF falls with depth so deeper
            // flow aliases earlier. this.nyquist stays the display baseline;
            // nyqEff drives the colour aliasing + wall-filter scaling only.
            const refDepth = this.depthCm * 0.5;
            const nyqEff = this.nyquist * clamp(refDepth / Math.max(depth, refDepth), 0.5, 1);
            // Wall filter scales with the effective Nyquist (PRF), so the colour
            // dead-zone shrinks with the velocity scale at depth.
            const wf = WALL_FILTER * nyqEff / this.nyquist;
            if (alos > wf) {           // wall filter: suppress near-zero flow
              // Colour-over-tissue priority keyed on the PRE-GAIN echogenicity e
              // (gain/TGC/log independent) rather than the fully-scaled brightness
              // b: only paint over near-anechoic blood, and changing the gain no
              // longer mutes the PISA/jet by pushing b past a brightness threshold.
              // Colour writes only over a genuine blood pool (LV/RV/LA/RA/AORTA)
              // AND where the pre-gain echo is truly anechoic, so the low-
              // echogenicity speckle tail of myocardium never takes colour.
              const isBlood = kind === TISSUE.LV || kind === TISSUE.RV ||
                              kind === TISSUE.LA || kind === TISSUE.RA ||
                              kind === TISSUE.AORTA;
              if (isBlood && e < 0.11) {
                // Position-seeded jitter dithers the aliasing wrap boundary so a
                // fast jet breaks into a TRUE mosaic instead of coherent hue-wrap
                // stripes. Two spatial-frequency terms (a low-frequency sequential
                // sample plus a high-frequency scrambled sample) sum toward
                // ~±0.5*nyqEff so the aliasing boundary shatters rather than bands.
                // Brightness in dopplerColor stays on the UNPERTURBED magnitude, so
                // this shifts only where the wrap lands, never how bright flow reads.
                const jLo = this.noise[(j * SW + i) % this.noise.length] - 0.5;
                const jHi = this.noise[((i * 131 + j * 57) & 0x7fffffff) % this.noise.length] - 0.5;
                const jitter = (jLo * 0.6 + jHi * 0.45) * nyqEff;
                const col = dopplerColor(los, nyqEff, vel.turbulent, jitter);
                // translucent mosaic: alpha grows with detectable velocity so
                // B-mode tissue stays visible under slow flow
                const mag = clamp((alos - wf) / (nyqEff - wf), 0, 1);
                dA = 0.40 + 0.34 * mag; // visible even for slow (normal) inflow
                dR = col[0]; dG = col[1]; dB = col[2];
              }
            }
          }
        }

        const ab = clamp((((theta + half) * halfInv * NA) | 0), 0, NA - 1);
        store(i, j, e, depth, kind, 1, ab, dR, dG, dB, dA);
      }
    }

    // adapt sampling stride only as a genuine slow-hardware fallback: coarsen
    // when sustained sampling cost is high, restore full resolution with headroom.
    if (t0) {
      const ms = performance.now() - t0;
      this._emaMs = this._emaMs * 0.85 + ms * 0.15;
      // Honest acquisition frame rate, measured from the real interval between
      // rendered frames (not the sampling cost). Scanners display this because
      // temporal resolution is a diagnostic property of the study, not a detail —
      // so it is reported here as achieved, never as an aspirational number.
      const tNow = performance.now();
      if (this._lastFrameT) {
        const iv = tNow - this._lastFrameT;
        if (iv > 0 && iv < 1000) this._emaFrameMs = this._emaFrameMs
          ? this._emaFrameMs * 0.9 + iv * 0.1 : iv;
      }
      this._lastFrameT = tNow;
      this.metrics.fps = this._emaFrameMs ? 1000 / this._emaFrameMs : null;
      // three tiers so a capable client renders the full-resolution buffer (stride 1)
      // while slow hardware coarsens the SAMPLING (2, then 3) without dropping frames.
      // Hysteresis band keeps it from oscillating between tiers. A pinned stride
      // (Detail control / tests) overrides the adaptation for a fixed resolution.
      if (this._pinStride) this._stride = this._pinStride;
      else if (this._emaMs > 30 && this._stride < 3) this._stride++;
      else if (this._emaMs < 14 && this._stride > 1) this._stride--;
    }

    // record the steered spectral gate for _spectral() (null => fall back to centre)
    this._peakGate = hasPeak ? {
      x: pkx, y: pky, z: pkz, los: pkLos, turbulent: pkTurb,
      speed: peakSpeed,               // true jet-core peak SPEED this frame (== metrics.peakVel)
      jet: pkTurb || peakSpeed > 1.8, // genuine jet: turbulent core, or clearly above inflow velocities
    } : null;

    // --- image formation: acoustic artifacts, PSF convolution, composite -------
    this._marchArtifacts(G, path, probe, half);   // shadowing / enhancement (polar)
    this._psfBlur();                               // depth/frequency separable PSF
    this._composite(data, cx, apexY, pxPerCm, half); // TGC + artifacts + log-compress + Doppler + gCNR

    this.bctx.putImageData(this.img, 0, 0);

    // paint to visible canvas
    const ctx = this.ctx;
    const W = this.canvas.width, H = this.canvas.height;
    ctx.fillStyle = '#05070c';
    ctx.fillRect(0, 0, W, H);
    ctx.imageSmoothingEnabled = true;
    ctx.imageSmoothingQuality = 'high';
    const scale = Math.min(W / SW, H / SH) * 0.98;
    const ox = (W - SW * scale) / 2;
    const oy = (H - SH * scale) / 2;
    ctx.drawImage(this.buf, ox, oy, SW * scale, SH * scale);

    const geo = { ox, oy, scale, cx, apexY, pxPerCm, half };
    // expose the current display geometry so an external measurement overlay can
    // convert canvas pixels ↔ centimetres. A straight-line displacement of `dpx`
    // display pixels is `dpx * cmPerPx` centimetres in the imaging plane (the sector
    // scales radial and lateral distance identically), so 2-D calipers are exact.
    this._geo = geo;
    this.cmPerPx = 1 / (scale * pxPerCm);
    this._overlays(ctx, geo);
    if (this.showLabels) this._labels(ctx, geo, G, path, probe);
    this._computeMetrics(G, path, peakSpeed, peakFlow);
    this._spectral(G, path, probe);
  }

  // March the sector in polar (angle x depth) space and accumulate two acoustic
  // artifact fields distal along each beam, gated by the tissue kinds already in
  // the model: (a) acoustic SHADOWING behind calcified/strongly-reflective valves
  // (attenuate everything deeper along the beam), and (b) posterior acoustic
  // ENHANCEMENT behind echo-free blood/effusion (brighten distal tissue). ~NA*NS
  // classify() calls — far fewer than the Cartesian sampling loop.
  _marchArtifacts(G, path, probe, half) {
    const NA = this._NA, NS = this._NS;
    const shadow = this._shadowGain, enh = this._enhGain;
    const P = probe.pos;
    const dStep = this.depthCm / NS;
    for (let a = 0; a < NA; a++) {
      const theta = -half + (a + 0.5) / NA * 2 * half;
      const bd = this.beamDir(probe, theta);
      const bx = bd[0], by = bd[1], bz = bd[2];
      let sh = 1;            // persistent distal attenuation (shadow)
      let fluidRun = 0;      // cm of contiguous anechoic column just traversed
      let enhAmt = 0;        // enhancement bonus captured on exiting a fluid column
      let sincePast = 1e9;   // cm since the fluid column ended
      const base = a * NS;
      for (let s = 0; s < NS; s++) {
        const depth = (s + 0.5) * dStep;
        const t = classify(P[0] + depth * bx, P[1] + depth * by, P[2] + depth * bz, G, path);
        const kind = t.tissue;
        // strong reflector: calcified valve (AS/MS) casts a distal shadow
        if (kind === TISSUE.VALVE && t.echo > 0.92) sh *= Math.exp(-2.4 * dStep);
        if (sh < 0.16) sh = 0.16; // shadow, not pure black
        const isFluid = kind === TISSUE.LV || kind === TISSUE.RV || kind === TISSUE.LA ||
                        kind === TISSUE.RA || kind === TISSUE.AORTA || kind === TISSUE.PERICARDIUM;
        if (isFluid) { fluidRun += dStep; sincePast = 0; }
        else {
          if (fluidRun > 0.6) enhAmt = Math.min(fluidRun * 0.24, 0.7); // brighten behind a real column
          fluidRun = 0; sincePast += dStep;
        }
        // enhancement applies to tissue behind the column, decaying with distance
        let en = 1;
        if (enhAmt > 0 && !isFluid) en = 1 + enhAmt * Math.exp(-sincePast / 2.2);
        shadow[base + s] = sh;
        enh[base + s] = en;
      }
    }
  }

  // Separable, spatially-varying Gaussian PSF convolution of the reflectivity
  // field: lateral (horizontal, across-beam) then axial (vertical, along-beam).
  // Both widths grow with depth (beam divergence) and shrink with frequency, and
  // lateral is always wider than axial — so the near field is crisp and the far
  // field blurs, lateral resolution is worse than axial, and the baked-in speckle
  // grain coarsens with depth exactly as real B-mode speckle tracks the PSF.
  _psfBlur() {
    const w = SW, h = SH, refl = this._refl, tmp = this._tmp, depthBuf = this._depth;
    const mask = this._mask, kern = this._kern;
    const freqScale = clamp(3.0 / this.freqMHz, 0.5, 2.0); // ref 3 MHz; higher freq => tighter PSF
    const [hLat, hAx] = this._harmonicScale();
    const latScale = freqScale * hLat, axScale = freqScale * hAx;
    const focus = this._focusCm();
    // The PSF radius is a function of DEPTH ALONE, so evaluating the two-way beam
    // (a sqrt and a divide) per pixel is pure waste — tabulate it once per frame
    // and let the hot loop do an array lookup instead.
    const LN = 256;
    if (!this._latLut) { this._latLut = new Uint8Array(LN); this._axLut = new Uint8Array(LN); }
    const latLut = this._latLut, axLut = this._axLut;
    const dMax = this.depthCm || 1;
    for (let i = 0; i < LN; i++) {
      const d = (i + 0.5) * dMax / LN;
      latLut[i] = this._latRadius(d, latScale, focus);
      axLut[i] = this._axRadius(d, axScale);
    }
    const lutK = LN / dMax;
    // lateral pass (across beams ~ horizontal): refl -> tmp
    for (let y = 0; y < h; y++) {
      const row = y * w;
      for (let x = 0; x < w; x++) {
        const idx = row + x;
        if (!mask[idx]) { tmp[idx] = 0; continue; }
        let li = (depthBuf[idx] * lutK) | 0; if (li >= LN) li = LN - 1;
        const r = latLut[li];
        if (r === 0) { tmp[idx] = refl[idx]; continue; }
        const k = kern[r]; let acc = 0, wsum = 0;
        for (let t = -r; t <= r; t++) {
          const xx = x + t;
          if (xx < 0 || xx >= w) continue;
          const j2 = row + xx;
          if (!mask[j2]) continue;
          const wv = k[t + r]; acc += refl[j2] * wv; wsum += wv;
        }
        tmp[idx] = wsum > 0 ? acc / wsum : refl[idx];
      }
    }
    // axial pass (along beam ~ vertical): tmp -> refl
    for (let y = 0; y < h; y++) {
      const row = y * w;
      for (let x = 0; x < w; x++) {
        const idx = row + x;
        if (!mask[idx]) { refl[idx] = 0; continue; }
        let ai = (depthBuf[idx] * lutK) | 0; if (ai >= LN) ai = LN - 1;
        const r = axLut[ai];
        if (r === 0) { refl[idx] = tmp[idx]; continue; }
        const k = kern[r]; let acc = 0, wsum = 0;
        for (let t = -r; t <= r; t++) {
          const yy = y + t;
          if (yy < 0 || yy >= h) continue;
          const j2 = yy * w + x;
          if (!mask[j2]) continue;
          const wv = k[t + r]; acc += tmp[j2] * wv; wsum += wv;
        }
        refl[idx] = wsum > 0 ? acc / wsum : tmp[idx];
      }
    }
  }

  // Form the displayed B-mode from the PSF-convolved reflectivity: depth
  // attenuation + TGC, acoustic artifacts (shadow/enhancement + near-field
  // reverberation), log compression, then composite the crisp Doppler overlay.
  // Simultaneously accumulate myocardium-vs-LV-blood grey samples for the live
  // image-quality readout (gCNR / CNR / speckle-SNR).
  _composite(data, cx, apexY, pxPerCm, half) {
    const refl = this._refl, depthBuf = this._depth, kindBuf = this._kind, mask = this._mask;
    const abin = this._abin, dcol = this._dcol, dalpha = this._dalpha;
    const shadow = this._shadowGain, enh = this._enhGain;
    const NS = this._NS, dInv = NS / this.depthCm;
    const gain = this.gain;
    // harmonic mode: steeper depth attenuation (the harmonic is attenuated at 2f),
    // a small overall gain-up as the scanner compensates, and near-total
    // suppression of near-field reverberation — those multiple internal echoes are
    // a FUNDAMENTAL-frequency artifact, and clearing them is harmonic imaging's
    // single most recognisable signature at the top of the sector.
    const pers = clamp(this.persistence || 0, 0, 0.85);
    if (pers > 0 && !this._prevGrey) this._prevGrey = new Float32Array(SW * SH);
    const prevGrey = this._prevGrey;
    const harm = this.harmonic;
    const atten = harm ? 0.030 : 0.023;
    const gainComp = harm ? 1.26 : 1.15;
    const revbGain = harm ? 0.16 : 1.0;
    // image-quality accumulators (target = myocardium, background = LV blood)
    const NB = 48; const hMyo = new Float32Array(NB), hLv = new Float32Array(NB);
    let nM = 0, sM = 0, sM2 = 0, nL = 0, sL = 0, sL2 = 0;
    for (let j = 0; j < SH; j++) {
      for (let i = 0; i < SW; i++) {
        const idx = j * SW + i;
        const k = idx * 4;
        if (!mask[idx]) { data[k] = 4; data[k + 1] = 6; data[k + 2] = 10; data[k + 3] = 255; continue; }
        const depth = depthBuf[idx];
        // depth attenuation (exp) partly compensated by a TGC ramp, + user gain.
        // Harmonic imaging pays a PENETRATION penalty: the harmonic is a
        // second-order effect, an order of magnitude weaker than the fundamental
        // and itself attenuated at twice the frequency, so the far field goes
        // softer. That trade-off (cleaner near/mid field, worse deep penetration)
        // is the honest reason harmonics are not simply "better" everywhere.
        const tgc = Math.exp(-depth * atten) * (1 + depth * 0.065);
        let b = refl[idx] * tgc * gain * gainComp;
        // acoustic artifacts (polar lookup): distal shadow / posterior enhancement
        let s = (depth * dInv) | 0; if (s >= NS) s = NS - 1;
        const abase = abin[idx] * NS + s;
        b *= shadow[abase] * enh[abase];
        // subtle near-field reverberation: faint repeating echoes close to the
        // transducer (multiple internal reflections), fading fast with depth.
        if (depth < 2.4 && revbGain > 0) {
          b += revbGain * 0.05 * Math.exp(-depth / 1.1) * (0.5 + 0.5 * Math.sin(depth * 9.0)) * (0.6 + 0.4 * refl[idx]);
        }
        b = Math.pow(clamp(b, 0, 1), 0.78); // log-ish compression
        // persistence: exponential blend with the previous frame's grey. Applied
        // AFTER compression so it averages what is displayed, as a scanner does.
        if (pers > 0) { b = b * (1 - pers) + prevGrey[idx] * pers; prevGrey[idx] = b; }
        // image-quality sampling on the B-mode grey (pre-Doppler)
        const kind = kindBuf[idx];
        if (kind === TISSUE.MYO) { const bi = b * NB | 0; hMyo[bi < NB ? bi : NB - 1]++; nM++; sM += b; sM2 += b * b; }
        else if (kind === TISSUE.LV) { const bi = b * NB | 0; hLv[bi < NB ? bi : NB - 1]++; nL++; sL += b; sL2 += b * b; }
        let R = b * 255, Gc = b * 255, Bc = b * 245;
        const a = dalpha[idx];
        if (a > 0) {
          const c3 = idx * 3;
          R = R * (1 - a) + dcol[c3] * a;
          Gc = Gc * (1 - a) + dcol[c3 + 1] * a;
          Bc = Bc * (1 - a) + dcol[c3 + 2] * a;
        }
        data[k] = R; data[k + 1] = Gc; data[k + 2] = Bc; data[k + 3] = 255;
      }
    }
    this._computeImageQuality(hMyo, hLv, NB, nM, sM, sM2, nL, sL, sL2);
  }

  // gCNR (generalized contrast-to-noise ratio, histogram-overlap form), CNR, and
  // speckle-SNR from the rendered B-mode — a live self-validation readout of image
  // quality between myocardium (target) and LV-cavity blood (background).
  _computeImageQuality(hMyo, hLv, NB, nM, sM, sM2, nL, sL, sL2) {
    const m = this.metrics;
    if (nM < 30 || nL < 30) { m.gcnr = null; m.cnr = null; m.speckleSNR = null; return; }
    // gCNR = 1 - OVL, where OVL is the overlap area of the two normalised histograms
    let ovl = 0;
    for (let b = 0; b < NB; b++) {
      const pm = hMyo[b] / nM, pl = hLv[b] / nL;
      ovl += pm < pl ? pm : pl;
    }
    m.gcnr = clamp(1 - ovl, 0, 1);
    const muM = sM / nM, muL = sL / nL;
    const vM = Math.max(0, sM2 / nM - muM * muM), vL = Math.max(0, sL2 / nL - muL * muL);
    m.cnr = Math.abs(muM - muL) / Math.sqrt(vM + vL + 1e-9);
    const sdM = Math.sqrt(vM);
    m.speckleSNR = sdM > 1e-6 ? muM / sdM : null; // fully-developed Rayleigh speckle ~1.91
  }

  _overlays(ctx, g) {
    const { ox, oy, scale, cx, apexY, pxPerCm, half } = g;
    ctx.save();
    ctx.strokeStyle = 'rgba(120,180,160,0.18)';
    ctx.fillStyle = 'rgba(150,220,190,0.65)';
    ctx.font = '10px monospace';
    ctx.lineWidth = 1;
    for (let d = 2; d <= this.depthCm; d += 2) {
      const rp = d * pxPerCm;
      const ex = cx + Math.sin(half) * rp, ey = apexY + Math.cos(half) * rp;
      ctx.beginPath();
      ctx.arc(ox + cx * scale, oy + apexY * scale, rp * scale, Math.PI / 2 - half, Math.PI / 2 + half);
      ctx.stroke();
      ctx.fillText(d + '', ox + (ex + 4) * scale, oy + ey * scale);
    }
    // Focal-zone marker on the depth scale — the caret every scanner puts beside
    // the depth ruler so the operator can see where the beam is narrowest.
    const fz = this._focusCm();
    if (fz > 0.5 && fz < this.depthCm) {
      const rp = fz * pxPerCm;
      const fx = ox + (cx + Math.sin(half) * rp) * scale;
      const fy = oy + (apexY + Math.cos(half) * rp) * scale;
      ctx.fillStyle = 'rgba(150,220,190,0.95)';
      ctx.beginPath();
      ctx.moveTo(fx + 2, fy);
      ctx.lineTo(fx + 9, fy - 4);
      ctx.lineTo(fx + 9, fy + 4);
      ctx.closePath();
      ctx.fill();
    }
    if (this.colorOn && this.showColorBox) {
      ctx.strokeStyle = 'rgba(90,200,255,0.5)';
      ctx.setLineDash([4, 3]);
      ctx.beginPath();
      ctx.moveTo(ox + cx * scale, oy + apexY * scale);
      const r2 = this.depthCm * 0.92 * pxPerCm;
      ctx.lineTo(ox + (cx + Math.sin(half * 0.92) * r2) * scale, oy + (apexY + Math.cos(half * 0.92) * r2) * scale);
      ctx.arc(ox + cx * scale, oy + apexY * scale, r2 * scale, Math.PI / 2 - half * 0.92, Math.PI / 2 + half * 0.92);
      ctx.lineTo(ox + cx * scale, oy + apexY * scale);
      ctx.stroke();
      ctx.setLineDash([]);
    }
    ctx.restore();
  }

  // Project each structure's 3D centroid onto the sector and label it, but only
  // for structures the imaging plane actually cuts through.
  _labels(ctx, g, G, path, probe) {
    const items = [
      ['LV', G.lv.c], ['RV', G.rv.c], ['LA', G.la.c], ['RA', G.ra.c], ['Ao', G.aorta.c],
      ['MV', VALVE_DEFS.mitral.c], ['TV', VALVE_DEFS.tricuspid.c], ['AV', VALVE_DEFS.aortic.c],
    ];
    // aortic-valve short axis: name the three cusps (R/L/N-coronary) at their centroids
    if (this.viewName === 'PSAX_AV') for (const c of AORTIC_CUSPS) items.push([c.name, c.p]);
    // LAA view: label the appendage body (mid-lobe of the appendage chain)
    if (this.viewName === 'MELAA' && G.A && G.A.la.aa && G.A.la.aa.length) {
      const s = G.A.la.aa[Math.floor(G.A.la.aa.length / 2)];
      items.push(['LAA', [(s.a[0] + s.b[0]) / 2, (s.a[1] + s.b[1]) / 2, (s.a[2] + s.b[2]) / 2]]);
    }
    // optional per-view whitelist: only label structures expected in this view
    const whitelist = this.viewName ? VIEW_WHITELIST[this.viewName] : null;
    ctx.save();
    ctx.font = '600 11px system-ui, sans-serif';
    ctx.textBaseline = 'middle';
    for (const [name, c] of items) {
      if (whitelist && !whitelist.includes(name)) continue;     // not expected in this view
      const rx = c[0] - probe.pos[0], ry = c[1] - probe.pos[1], rz = c[2] - probe.pos[2];
      const along = rx * probe.dir[0] + ry * probe.dir[1] + rz * probe.dir[2];
      const lat = rx * probe.lat[0] + ry * probe.lat[1] + rz * probe.lat[2];
      const nrm = rx * probe.normal[0] + ry * probe.normal[1] + rz * probe.normal[2];
      if (along <= 0.3 || along > this.depthCm) continue;      // out of depth
      if (Math.abs(Math.atan2(lat, along)) > this.sectorHalf * 0.96) continue; // out of sector
      // plane-cut test: a structure whitelisted for the active view is expected
      // in-plane, so allow a looser out-of-plane tolerance (chambers are large
      // and their centroid can sit a little off the exact scan plane).
      const nrmTol = (whitelist && whitelist.includes(name)) ? 1.7 : 0.8;
      if (Math.abs(nrm) > nrmTol) continue;                     // plane doesn't truly cut it
      const rpix = Math.hypot(along, lat) * g.pxPerCm;
      const theta = Math.atan2(lat, along);
      const bx = g.cx + Math.sin(theta) * rpix;
      const by = g.apexY + Math.cos(theta) * rpix;
      const px = g.ox + bx * g.scale, py = g.oy + by * g.scale;
      const w = ctx.measureText(name).width + 8;
      ctx.fillStyle = 'rgba(8,14,22,0.62)';
      ctx.fillRect(px - w / 2, py - 8, w, 16);
      ctx.fillStyle = 'rgba(120,230,255,0.95)';
      ctx.textAlign = 'center';
      ctx.fillText(name, px, py + 1);
    }
    ctx.restore();
  }

  _computeMetrics(G, path, peakSpeed, peakFlow) {
    const vol = (r) => (4 / 3) * Math.PI * r[0] * r[1] * r[2];
    // the two fixed-phase geometries only depend on `path`; memoise them and
    // recompute only when the path/pathology object identity changes.
    if (path !== this._lastMetricsPath) {
      this._ged = geometryAt(0.0, path);   // end-diastole
      this._ges = geometryAt(0.32, path);  // end-systole
      this._lastMetricsPath = path;
      this._peakVelMax = 0;                 // reset the peak-velocity tracker
    }
    const ged = this._ged, ges = this._ges;
    const edv = vol(ged.lv.r), esv = vol(ges.lv.r);
    this.metrics.lvidd = 2 * ged.lv.r[0];
    this.metrics.lvids = 2 * ges.lv.r[0];
    this.metrics.ef = edv > 0 ? ((edv - esv) / edv) * 100 : null;
    // MAPSE — mitral annular plane systolic excursion (mm): the LV base descends
    // toward the fixed apex in systole (apex-anchored contraction). Normal ≥ 10 mm;
    // a longitudinal-function marker reduced in DCM / infarction.
    const annED = ged.lv.c[1] + ged.lv.r[1], annES = ges.lv.c[1] + ges.lv.r[1];
    this.metrics.mapse = (annED - annES) * 10; // cm → mm
    this.metrics.peakVel = peakSpeed > 0.1 ? peakSpeed : null; // true jet-core speed magnitude
    this.metrics.peakGrad = peakSpeed > 0.1 ? 4 * peakSpeed * peakSpeed : null; // simplified Bernoulli
    // which flow compartment produced this frame's peak speed (teaching context)
    this.metrics.peakLabel = peakSpeed > 0.1 ? (FLOW_LABELS[peakFlow] || null) : null;
    this.metrics.effusion = path.effusion ? 1.6 : null;

    // --- severity grading (BSE/ASE valve thresholds), surfaced to the UI ---
    // Grade off the peak jet velocity seen since the pathology last changed: the
    // per-frame peak dips when the jet isn't firing, so the max over the cycle is
    // what the grade should reflect. Deterministic per case.
    if (peakSpeed > this._peakVelMax) this._peakVelMax = peakSpeed;
    const pv = this._peakVelMax;
    // regurgitations/shunts carry a severity grade (mild/moderate/severe); default
    // severe so single-severity presets are unchanged. AS/PS self-grade off velocity.
    const gr = (path.grade === 'mild' || path.grade === 'moderate' || path.grade === 'severe') ? path.grade : 'severe';
    let severity = null, severityLabel = null;
    if (path.aorticStenosis) {
      // grade off the model's continuity peak (what a CW sweep captures) rather than
      // the plane-dependent 2-D sample, so the band is robust to the imaging window.
      const av = Math.max(pv, (G.hemo && G.hemo.vAoPeak) || 0);
      // peak velocity bands (peak gradient 4v^2: <20 / 20-39 / 40-59 / >=60 mmHg)
      if (av >= 5.0) severity = 'very severe';
      else if (av >= 4.0) severity = 'severe';
      else if (av >= 3.0) severity = 'moderate';
      else if (av >= 2.5) severity = 'mild';
      if (severity) severityLabel = severity + ' AS';
    } else if (path.pulmonaryStenosis || path.pulmonicStenosis || path.ps) {
      if (pv > 4) severity = 'severe';
      else if (pv >= 3) severity = 'moderate';
      else if (pv > 0.1) severity = 'mild';
      if (severity) severityLabel = severity + ' PS';
    } else if (path.mr) {            // graded by regurgitant fraction in the model
      severity = gr; severityLabel = gr + ' MR';
    } else if (path.tr) {
      severity = gr; severityLabel = gr + ' TR';
    } else if (path.vsd) {          // modelled shunts are large by design
      severity = 'large shunt'; severityLabel = 'large VSD';
    } else if (path.asd) {
      severity = 'large shunt'; severityLabel = 'large ASD';
    } else if (path.mitralStenosis) {
      // graded by the modelled diastolic mean gradient / valve area
      severity = gr; severityLabel = gr + ' MS';
    }
    // normal / effusion / dcm: severity null (dcm severity comes from EF in the UI)
    this.metrics.severity = severity;
    this.metrics.severityLabel = severityLabel;
    // estimated PA systolic pressure from the peak TR jet: PASP = 4·V_TR² + RAP.
    // Use the model's TR jet velocity (what a CW sweep captures) rather than the
    // plane-dependent 2-D sample. RAP = 10 mmHg with a dilated/pressure-loaded RV, else 5.
    const vTR = (G.hemo && G.hemo.vTRPeak) || 0;
    if (path.tr && vTR > 1.5) {
      const rap = path.rvpo ? 10 : 5;
      this.metrics.pasp = 4 * vTR * vTR + rap;
    } else {
      this.metrics.pasp = null;
    }
  }

  _spectral(G, path, probe) {
    if (!this.sctx) return;
    const ctx = this.sctx;
    const W = this.spectral.width, H = this.spectral.height;
    const nyq = this.nyquist;
    const zero = H * 0.5;
    const amp = H * 0.46;

    // Fixed sample gate: this frame's steered peak-flow location if there is a
    // jet, else a central mid-beam gate. The beam direction at the gate gives
    // the line-of-sight projection.
    let gx, gy, gz;
    const gate = this._peakGate;
    if (gate) { gx = gate.x; gy = gate.y; gz = gate.z; }
    else {
      const bd0 = this.beamDir(probe, 0);
      const gd = Math.min(this.depthCm * 0.4, 7.5);
      gx = probe.pos[0] + gd * bd0[0]; gy = probe.pos[1] + gd * bd0[1]; gz = probe.pos[2] + gd * bd0[2];
    }
    const bx = gx - probe.pos[0], by = gy - probe.pos[1], bz = gz - probe.pos[2];
    const bl = Math.hypot(bx, by, bz) || 1;
    const bdx = bx / bl, bdy = by / bl, bdz = bz / bl; // unit gate beam direction

    // Sweep the WHOLE cardiac cycle every frame so the trace is always a full,
    // dense velocity-time waveform (E/A inflow, systolic jet envelope, …) that
    // resets cleanly per case — no stale history carried between pathologies.
    const step = 2;
    const cols = Math.ceil(W / step);
    const vs = this._sweepV || (this._sweepV = new Float32Array(2048));
    const tb = this._sweepT || (this._sweepT = new Uint8Array(2048));
    // A genuine jet (turbulent core / clearly above inflow) is displayed like CW:
    // the trace peak must read the frame's TRUE jet-core speed (the value feeding
    // metrics.peakVel), not the gated line-of-sight sample — a jet running across
    // the beam or straddling a thin vena contracta under-reads otherwise. Normal
    // laminar inflow keeps the honest gated line-of-sight (PW) sample.
    const isJet = !!(gate && gate.jet);
    const coreSpeed = gate ? gate.speed : 0;
    let peakAbs = 0;
    for (let c = 0; c < cols; c++) {
      const ph = c / cols;
      const Gp = geometryAt(ph, path);
      const vel = velocityAt(gx, gy, gz, Gp, path);
      let v = 0;
      if (vel) {
        const losb = vel.vx * bdx + vel.vy * bdy + vel.vz * bdz;
        // jet: angle-corrected true core speed, signed by beam projection (CW-like);
        // laminar inflow: honest line-of-sight PW sample.
        v = isJet ? (losb < 0 ? -1 : 1) * vel.speed : losb;
      }
      vs[c] = v; tb[c] = vel && vel.turbulent ? 1 : 0;
      const a = Math.abs(v); if (a > peakAbs) peakAbs = a;
    }
    // Scale the swept jet envelope so its peak equals the tracked true jet-core
    // speed: a severe MR/AS jet then shows ~5.5 / ~4.7 m/s on the trace, matching
    // the Measurements Peak V, instead of the gate straddling the vena contracta.
    if (isJet && coreSpeed > 0.1 && peakAbs > 1e-3) {
      const k = coreSpeed / peakAbs;
      for (let c = 0; c < cols; c++) vs[c] *= k;
      peakAbs = coreSpeed;
    }
    const aliased = peakAbs > nyq + 1e-6;
    const scaleV = Math.max(nyq, peakAbs) * 1.12; // auto-scale to fit the peak
    const mode = aliased ? 'CW' : 'PW';

    // background + calibrated gridlines + baseline
    ctx.fillStyle = '#04060a';
    ctx.fillRect(0, 0, W, H);
    ctx.strokeStyle = 'rgba(90,130,170,0.16)';
    ctx.lineWidth = 1;
    for (let q = 1; q <= 2; q++) {
      const dy = (q / 2) * amp;
      ctx.beginPath(); ctx.moveTo(0, zero - dy); ctx.lineTo(W, zero - dy); ctx.stroke();
      ctx.beginPath(); ctx.moveTo(0, zero + dy); ctx.lineTo(W, zero + dy); ctx.stroke();
    }
    ctx.strokeStyle = 'rgba(120,160,200,0.4)';
    ctx.beginPath(); ctx.moveTo(0, zero); ctx.lineTo(W, zero); ctx.stroke();

    // filled velocity envelope (toward probe = above baseline)
    for (let c = 0; c < cols; c++) {
      const v = vs[c];
      const y = clamp(zero - (-v / scaleV) * amp, 1, H - 1);
      const bright = clamp(Math.abs(v) / scaleV, 0, 1);
      if (bright < 0.02) continue;
      const x = c * step;
      // filled column from baseline to the envelope, brighter with velocity
      ctx.fillStyle = `rgba(${150 + 80 * bright},${210 + 30 * bright},${180},${0.22 + 0.5 * bright})`;
      ctx.fillRect(x, Math.min(zero, y), step, Math.abs(y - zero));
      // spectral broadening band for turbulent flow
      if (tb[c]) {
        const sp = amp * 0.32 * bright;
        ctx.fillStyle = `rgba(120,235,150,${0.10 + 0.25 * bright})`;
        ctx.fillRect(x, y - sp * 0.5, step, sp);
      }
    }
    // bright envelope contour
    ctx.strokeStyle = 'rgba(200,240,210,0.85)';
    ctx.lineWidth = 1.2;
    ctx.beginPath();
    for (let c = 0; c < cols; c++) {
      const y = clamp(zero - (-vs[c] / scaleV) * amp, 1, H - 1);
      const x = c * step;
      if (c === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
    }
    ctx.stroke();

    // moving time cursor at the current cardiac phase
    const cxp = (G.phase % 1) * W;
    ctx.strokeStyle = 'rgba(255,90,77,0.8)';
    ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(cxp, 0); ctx.lineTo(cxp, H); ctx.stroke();

    // labels: scale, mode (PW/CW), sweep + time
    ctx.fillStyle = 'rgba(150,200,240,0.7)';
    ctx.font = '9px monospace';
    ctx.textAlign = 'left';
    ctx.fillText(mode + '  m/s', 4, 11);
    ctx.textAlign = 'right';
    ctx.fillText('+' + scaleV.toFixed(1), W - 4, 11);
    ctx.fillText('-' + scaleV.toFixed(1), W - 4, H - 4);
    ctx.textAlign = 'left';
    ctx.fillStyle = 'rgba(150,200,240,0.5)';
    ctx.fillText('1 cardiac cycle', 4, H - 4);

    // peak-velocity caliper (true peak, always legible even when colour aliases)
    if (peakAbs > 0.12) {
      const py = zero - (peakAbs / scaleV) * amp;
      ctx.strokeStyle = 'rgba(255,210,90,0.7)';
      ctx.setLineDash([3, 3]);
      ctx.beginPath(); ctx.moveTo(0, py); ctx.lineTo(W, py); ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = 'rgba(255,225,140,0.95)';
      ctx.textAlign = 'right';
      ctx.font = '600 10px monospace';
      ctx.fillText('peak ' + peakAbs.toFixed(1) + ' m/s' + (aliased ? ' (aliased in colour)' : ''), W - 6, py - 3 < 8 ? py + 12 : py - 3);
    }
  }
}

// Map a line-of-sight velocity (m/s) to an RGB colour using BART convention.
// Negative los = towards probe = RED/orange. Positive = away = BLUE/cyan.
// Aliasing is decoupled from brightness: the WRAPPED value only picks the hue /
// red-vs-blue direction (so a fast jet still wraps into a mosaic), while the
// TRUE (pre-wrap) flow strength drives intensity so a fast jet stays bright
// instead of going black at every Nyquist multiple. Turbulence adds green.
function dopplerColor(los, nyq, turbulent, jitter) {
  const raw = Number.isFinite(los) ? los : 0;
  const jit = Number.isFinite(jitter) ? jitter : 0;
  // wrapped value chooses hue/direction (constant time, finite-safe). A small
  // position-seeded jitter perturbs los ONLY for the wrap decision, dithering the
  // aliasing boundary so hue-wrap stripes dissolve into a natural mosaic. The
  // brightness envelope below uses the UNPERTURBED strength, so the jitter shifts
  // only where the wrap lands, never how bright the flow reads.
  const wrapArg = raw + jit;
  const v = wrapArg - 2 * nyq * Math.round(wrapArg / (2 * nyq));
  const mag = clamp(Math.abs(v) / nyq, 0, 1);       // wrapped magnitude -> hue saturation
  const trueMag = clamp(Math.abs(raw) / nyq, 0, 1); // pre-wrap strength -> brightness
  const i = 0.45 + 0.55 * trueMag; // intensity envelope (bright for fast jets)
  let r, g, b;
  if (v < 0) {
    r = (150 + 105 * mag) * i; g = (20 + 210 * mag * mag) * i; b = 25 * i;
  } else {
    r = 25 * i; g = (30 + 200 * mag * mag) * i; b = (150 + 105 * mag) * i;
  }
  if (turbulent) {
    const t = 0.5;
    r = r * (1 - t) + 120 * t;
    g = g * (1 - t) + 235 * t;
    b = b * (1 - t) + 90 * t;
  }
  return [r, g, b];
}
