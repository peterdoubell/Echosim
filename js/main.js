// main.js — orchestrates the cardiac model, 3D anatomy view, 2D echo display,
// the moveable ultrasound plane, the cardiac-cycle clock and all UI controls.

import { geometryAt, hemoSummary, pvLoop, longitudinalStrain, regionalStrain, eaRatio } from './cardiac-model.js';
// BSE/ASE normal reference ranges — single source of truth for the
// Measurements thresholds and the displayed "normal band" text.
// (cardiac-model.js does not re-export REF, so import it directly.)
import { REF } from './anatomy.js';
import { Heart3D } from './heart3d.js';
import { EchoView } from './echo.js';
import { MeasureTool } from './measure.js';
import { vadd, vsub, vscale, vdot, vcross, norm, vrot } from './mathutils.js';
import { TTE_VIEWS, EXTRA_VIEWS, TEE_VIEWS, makeProbe } from './views.js';
import { Heartbeat, vcgToLeads, LEAD_NAMES } from './electrophysiology.js';

// ---------------------------------------------------------------------------
// Standard echocardiographic windows: derived from the anatomy landmarks in
// views.js, each displayed in the ASE/EACVI screen orientation (e.g. A4C with
// the LV on the right, PLAX with the aorta on the right, PSAX viewed from the
// apex with the RV on the left) at real clinical imaging depths.
// ---------------------------------------------------------------------------
const VIEWS = Object.fromEntries(Object.entries(TTE_VIEWS).map(([k, v]) => [k, v.probe]));
const viewDepth = (name) => (TTE_VIEWS[name] || EXTRA_VIEWS[name] || {}).depth || 16;
// friendly labels for the on-image view tag (raw keys may contain underscores)
const VIEW_LABEL = {
  PSAX_AV: 'PSAX-AV', MELAA: 'ME LAA', PSAX_MV: 'PSAX-MV', RVIT: 'RV inflow', SC_IVC: 'Subcostal IVC',
  SSN: 'Suprasternal', ME2C: 'ME 2C', MEAVSAX: 'ME AV SAX', MEBICAVAL: 'ME bicaval', MERVIO: 'ME RV in-out',
  DESCAO: 'Desc Ao SAX',
};

// Transesophageal (TEE) window set (views.js): a higher-frequency probe in the
// oesophagus directly behind the left atrium, so the atria sit in the near field.
const TEE_CAM = { ME4C: 'A4C', MELAX: 'PLAX', TGSAX: 'PSAX', MELAA: 'A2C', ME2C: 'A2C', MEAVSAX: 'PSAX_AV',
  MEBICAVAL: 'SUBCOSTAL', MERVIO: 'PSAX_AV', DESCAO: 'PSAX' }; // reuse a sensible 3D camera

const PATHOLOGY = {
  normal:    { flags: {}, view: 'PLAX', name: 'Normal heart',
    text: 'Four chambers contract and relax in sequence. Valves open and close cleanly; laminar inflow (<b>red</b>, toward the apex in diastole) and outflow (<b>away</b> in systole) stay below the Nyquist limit. <span class="cue">Baseline for comparison.</span>' },
  mr: { flags: { mr: true }, view: 'A4C', name: 'Mitral regurgitation',
    text: 'An incompetent mitral valve lets blood leak <b>backwards</b> from LV to LA during systole. <span class="cue">Look for a turbulent (mosaic) high-velocity jet firing into the left atrium in systole.</span>' },
  ms: { flags: { mitralStenosis: true }, view: 'A4C', name: 'Mitral stenosis',
    text: 'A narrowed, thickened mitral valve restricts LV filling. <span class="cue">Bright/calcified leaflets, a small orifice, and accelerated, prolonged diastolic inflow velocities.</span>' },
  as: { flags: { aorticStenosis: true, lvh: true }, view: 'PLAX', name: 'Aortic stenosis',
    text: 'A calcified aortic valve barely opens, so the LV must generate high pressure. <span class="cue">Restricted valve motion, a high-velocity systolic ejection jet, and a thickened (hypertrophied) LV wall.</span>' },
  dcm: { flags: { dilated: true }, view: 'A4C', name: 'Dilated cardiomyopathy',
    text: 'The left ventricle is enlarged and contracts poorly (low ejection fraction). <span class="cue">Big rounded LV with little wall thickening through the cycle and sluggish, low-velocity flow.</span>' },
  rwma: { flags: { rwma: 'septal' }, view: 'A4C', name: 'Regional wall-motion abnormality',
    text: 'After a myocardial infarction one territory stops contracting. <span class="cue">Watch the <b>septal</b> wall fail to thicken and move inward in systole while the other walls contract normally — regional, not global, dysfunction. Overall EF is reduced.</span>' },
  effusion: { flags: { effusion: true }, view: 'PLAX', name: 'Pericardial effusion',
    text: 'Fluid collects in the pericardial sac. <span class="cue">An echo-free (black) space surrounding the heart, typically first seen posteriorly/inferiorly.</span>' },
  vsd: { flags: { vsd: true }, view: 'A4C', name: 'Ventricular septal defect',
    text: 'A hole in the interventricular septum shunts blood LV → RV. <span class="cue">A turbulent systolic jet crossing the septum at high velocity (left-to-right shunt).</span>' },
  asd: { flags: { asd: true }, view: 'SUBCOSTAL', name: 'Atrial septal defect',
    text: 'A defect in the atrial septum shunts blood LA → RA. <span class="cue">Low-velocity flow crossing the atrial septum; best appreciated from the subcostal window.</span>' },
  tr: { flags: { tr: true }, view: 'A4C', name: 'Tricuspid regurgitation',
    text: 'The tricuspid valve leaks during systole, sending blood RV → RA. <span class="cue">A systolic jet into the right atrium; its peak velocity is used to estimate pulmonary pressures.</span>' },
  phtn: { flags: { rvpo: true, tr: true }, view: 'PSAX', name: 'Pulmonary hypertension (RV overload)',
    text: 'Chronically raised pulmonary pressures dilate and pressure-load the right ventricle, flattening the interventricular septum. <span class="cue">PSAX: the LV is <b>D-shaped</b> (the "D-sign") and the RV is dilated; a high-velocity tricuspid regurgitant jet estimates PA systolic pressure.</span>' },
};

// Per-pathology "what to look for" checklist + a guided-view hint.
const CHECKS = {
  normal:   { look: ['All four chambers', 'Valves opening & closing cleanly', 'Laminar inflow/outflow under Nyquist'], hint: 'Sweep through the cycle and watch red diastolic inflow flip to systolic ejection.' },
  mr:       { look: ['Mitral valve leaflets', 'Systolic mosaic jet into the LA', 'Dilated left atrium'], hint: 'Bring the mitral valve into the beam and pause in systole to see the jet fill the LA.' },
  ms:       { look: ['Thickened/bright mitral leaflets', 'Small orifice', 'Accelerated diastolic inflow'], hint: 'Watch the diastolic inflow velocity on the spectral trace — it stays high and prolonged.' },
  as:       { look: ['Restricted aortic valve', 'Thick (hypertrophied) LV wall', 'High-velocity systolic jet'], hint: 'Note the peak velocity & ΔP — severe AS is >4 m/s (>64 mmHg).' },
  dcm:      { look: ['Enlarged rounded LV', 'Little wall thickening', 'Low EF'], hint: 'Compare LVIDd and EF in Measurements against a normal heart.' },
  rwma:     { look: ['One akinetic wall segment', 'Normal thickening elsewhere', 'Reduced global EF'], hint: 'Pause in systole and watch the septum fail to thicken while the lateral wall does.' },
  effusion: { look: ['Echo-free rim around the heart', 'Largest posteriorly', 'Swinging heart if large'], hint: 'Look behind the LV posterior wall for the dependent black crescent.' },
  vsd:      { look: ['Interventricular septum', 'Systolic jet LV→RV', 'Turbulent mosaic'], hint: 'Centre the septum and pause in systole to catch the trans-septal jet.' },
  asd:      { look: ['Interatrial septum', 'Flow LA→RA', 'Dilated right heart'], hint: 'Use the subcostal window where the atrial septum lies across the beam.' },
  tr:       { look: ['Tricuspid valve', 'Systolic jet into the RA', 'Peak jet velocity → PA pressure'], hint: 'Peak TR velocity estimates pulmonary pressures — read it off the trace.' },
  phtn:     { look: ['Dilated RV', 'D-shaped (flattened) septum in PSAX', 'High-velocity TR jet'], hint: 'Switch to PSAX for the D-sign; the TR jet peak estimates PA systolic pressure.' },
};

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------
const state = {
  viewName: 'PLAX',
  baseProbe: VIEWS.PLAX(),
  path: {},
  phase: 0,
  playing: true,
  hr: 72,
  az: 0, tilt: 0, slide: 0,
  lastT: performance.now(),
};

// adjust the base probe by the manual azimuth/tilt/slide controls
function effectiveProbe(A) {
  let { pos, dir, lat, normal, target } = state.baseProbe;
  pos = pos.slice(); dir = dir.slice(); lat = lat.slice(); normal = normal.slice();
  // views that follow a moving landmark (PSAX-MV tracks the mitral annulus)
  const vw = TTE_VIEWS[state.viewName] || EXTRA_VIEWS[state.viewName] || TEE_VIEWS[state.viewName];
  const trk = A && vw && vw.track;
  if (trk) { const o = trk(A); pos = vadd(pos, o); target = vadd(target, o); }
  const Y = [0, 1, 0];
  if (state.az) {
    const a = state.az * Math.PI / 180;
    pos = vadd(target, vrot(vsub(pos, target), Y, a));
    dir = vrot(dir, Y, a); lat = vrot(lat, Y, a); normal = vrot(normal, Y, a);
  }
  if (state.tilt) {
    const a = state.tilt * Math.PI / 180;
    pos = vadd(target, vrot(vsub(pos, target), lat, a));
    dir = vrot(dir, lat, a); normal = vrot(normal, lat, a);
  }
  if (state.slide) {
    pos = vadd(pos, vscale(lat, state.slide * 0.12));
  }
  return { pos, dir: norm(dir), lat: norm(lat), normal: norm(normal), target };
}

// ---------------------------------------------------------------------------
// Setup
// ---------------------------------------------------------------------------
const heart3d = new Heart3D(document.getElementById('scene3d'));
const echo = new EchoView(document.getElementById('echo'), document.getElementById('spectral'));
// learner-placed caliper / Simpson-EF overlay, checked against model ground truth
const measure = new MeasureTool(document.getElementById('echoOverlay'), echo, {
  groundTruth: () => ({ ef: echo.metrics && echo.metrics.ef, lvidd: echo.metrics && echo.metrics.lvidd }),
  onResult: (r) => publishMeasure(r),
});
const ecgCanvas = document.getElementById('ecg');
const ecgCtx = ecgCanvas.getContext('2d');
// expose the core objects for debugging, headless tests and third-party extension
window.echosim = {
  echo, measure, get state() { return state; },
  // set an arbitrary probe (for TEE view calibration): makeProbe(target,dir,normal,standoff)
  setProbe: (t, d, n, s) => { state.viewName = 'CUSTOM'; state.baseProbe = makeProbe(t, d, n, s); state.az = state.tilt = state.slide = 0; },
  // select any named view, including the extra apical/parasternal ones (A5C, A3C, PSAX_MV)
  setView: (name) => setView(name),
};

// ---------------------------------------------------------------------------
// Electrophysiology: a VCG dipole → 12-lead engine is the master clock. It
// schedules beats (rhythm) and its ventricular-activation time drives the
// mechanical phase (electromechanical coupling), so the ECG and echo are one
// coupled system. We keep a rolling buffer of the VCG so any lead (and the full
// 12-lead grid) can be drawn from the same synchronized source.
// ---------------------------------------------------------------------------
const heart = new Heartbeat();
const ECG_FS = 250;                 // ECG sample rate (Hz) for the rolling buffer
const ECG_SECONDS = 6;
const vcgBuf = new Float32Array(ECG_FS * ECG_SECONDS * 3); // ring of [x,y,z]
let vcgHead = 0, vcgFilled = 0;
let ecgAccum = 0;                   // time accumulator to sample at ECG_FS
let ecgLead = 'II';                 // lead shown in the top-bar rhythm strip

function pushVcg(v) {
  vcgBuf[vcgHead * 3] = v[0]; vcgBuf[vcgHead * 3 + 1] = v[1]; vcgBuf[vcgHead * 3 + 2] = v[2];
  vcgHead = (vcgHead + 1) % (ECG_FS * ECG_SECONDS);
  if (vcgFilled < ECG_FS * ECG_SECONDS) vcgFilled++;
}
// oldest→newest ordered index helper
function vcgSample(i, out) {
  const cap = ECG_FS * ECG_SECONDS;
  const start = (vcgHead - vcgFilled + cap) % cap;
  const k = (start + i) % cap;
  out[0] = vcgBuf[k * 3]; out[1] = vcgBuf[k * 3 + 1]; out[2] = vcgBuf[k * 3 + 2];
  return out;
}

// ---------------------------------------------------------------------------
// 3D re-bake loading overlay (contract with heart3d: setPathology is async and
// fires onBakeStart / onBakeEnd around the marching-cubes re-bake).
// ---------------------------------------------------------------------------
const scene3dEl = document.getElementById('scene3d');
const rebuildOverlay = document.getElementById('rebuildOverlay');
heart3d.onBakeStart = () => {
  if (rebuildOverlay) rebuildOverlay.classList.remove('hidden');
  if (scene3dEl) { scene3dEl.setAttribute('aria-busy', 'true'); scene3dEl.style.cursor = 'wait'; }
};
heart3d.onBakeEnd = () => {
  if (rebuildOverlay) rebuildOverlay.classList.add('hidden');
  if (scene3dEl) { scene3dEl.setAttribute('aria-busy', 'false'); scene3dEl.style.cursor = ''; }
};

// Keep the 2D canvas backing-stores matched to their CSS size × DPR so the
// B-mode and spectral traces stay crisp on large / high-DPI displays. The
// render loop reads canvas.width/height each frame, so resizing is safe.
function fitCanvas(cv) {
  if (!cv) return;
  const w = cv.clientWidth, h = cv.clientHeight;
  if (!w || !h) return; // hidden / zero-size — skip
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  const bw = Math.round(w * dpr), bh = Math.round(h * dpr);
  if (cv.width !== bw || cv.height !== bh) { cv.width = bw; cv.height = bh; }
}
function resizeCanvases() {
  fitCanvas(document.getElementById('echo'));
  // keep the measurement overlay's backing store identical to #echo's, so a
  // caliper pixel maps through the same echo.cmPerPx to centimetres.
  fitCanvas(document.getElementById('echoOverlay'));
  fitCanvas(document.getElementById('spectral'));
  // #ecg is left at its fixed backing-store size.
}
let _resizeTimer = null;
window.addEventListener('resize', () => {
  clearTimeout(_resizeTimer);
  _resizeTimer = setTimeout(resizeCanvases, 150);
});

function phaseName(p) {
  if (p < 0.05) return 'Isovolumic contraction';
  if (p < 0.34) return 'Systole · ejection';
  if (p < 0.42) return 'Isovolumic relaxation';
  if (p < 0.55) return 'Early diastole · E wave';
  if (p < 0.84) return 'Diastasis';
  return 'Atrial kick · A wave';
}

// Draw one lead from the rolling VCG buffer onto a canvas as a scrolling strip.
// mmPerMv scales amplitude; secShown sets the time window across the width.
const _vs = [0, 0, 0];
function drawLead(ctx, W, H, leadName, secShown, opts = {}) {
  ctx.clearRect(0, 0, W, H);
  if (opts.grid) {
    // faint ECG graph paper (major squares ~0.2 s)
    ctx.strokeStyle = 'rgba(255,120,120,0.10)'; ctx.lineWidth = 1;
    const px = W / secShown, step = 0.2 * px;
    ctx.beginPath();
    for (let x = 0; x <= W; x += step) { ctx.moveTo(x, 0); ctx.lineTo(x, H); }
    for (let y = 0; y <= H; y += H / 6) { ctx.moveTo(0, y); ctx.lineTo(W, y); }
    ctx.stroke();
  }
  const nShow = Math.min(vcgFilled, Math.floor(secShown * ECG_FS));
  const base = H * 0.6, gain = (H * 0.32) / 1.2; // ~1.2 mV to 0.32 H
  ctx.strokeStyle = opts.color || '#35d0a0';
  ctx.lineWidth = opts.lw || 1.5;
  ctx.beginPath();
  const startIdx = vcgFilled - nShow;
  for (let s = 0; s < nShow; s++) {
    vcgSample(startIdx + s, _vs);
    const val = vcgToLeads(_vs)[leadName];
    const x = (s / (nShow - 1 || 1)) * W;
    const y = base - val * gain;
    if (s === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
  }
  ctx.stroke();
  if (opts.label) {
    ctx.fillStyle = 'rgba(200,230,255,0.85)'; ctx.font = 'bold 10px system-ui';
    ctx.fillText(leadName, 4, 11);
  }
}

function drawEcg() {
  drawLead(ecgCtx, ecgCanvas.width, ecgCanvas.height, ecgLead, ECG_SECONDS * 0.62, { label: true });
}

// Full 12-lead grid (3 rows × 4 cols) drawn into the #twelveLead modal canvas.
function drawTwelveLead() {
  const cv = document.getElementById('twelveLead');
  if (!cv || cv.offsetParent === null) return; // hidden
  const ctx = cv.getContext('2d');
  const W = cv.width, H = cv.height;
  ctx.fillStyle = '#0a0f18'; ctx.fillRect(0, 0, W, H);
  const cols = 4, rows = 3, cw = W / cols, ch = H / rows, secShown = 2.4;
  ctx.save();
  for (let r = 0; r < rows; r++) for (let c = 0; c < cols; c++) {
    const name = LEAD_NAMES[r * cols + c];
    ctx.save();
    ctx.translate(c * cw, r * ch);
    ctx.beginPath(); ctx.rect(0, 0, cw, ch); ctx.clip();
    drawLead(ctx, cw, ch, name, secShown, { grid: true, label: true, color: '#39e0a0', lw: 1.2 });
    ctx.strokeStyle = 'rgba(120,150,190,0.18)'; ctx.strokeRect(0.5, 0.5, cw - 1, ch - 1);
    ctx.restore();
  }
  ctx.restore();
}

// ---------------------------------------------------------------------------
// Main loop
// ---------------------------------------------------------------------------
let lastBannerT = 0;   // throttle for the probe-banner frame-rate readout
function loop(now) {
  const dt = Math.min((now - state.lastT) / 1000, 0.05);
  state.lastT = now;

  if (state.playing) {
    // Advance the electrophysiology engine at a fixed ECG sample rate and buffer
    // the VCG; the mechanical phase is derived from ventricular depolarisation
    // (electromechanical coupling), so an irregular rhythm (AF, PVC) makes the
    // heart beat — and the echo move — irregularly.
    const step = 1 / ECG_FS;
    ecgAccum += dt;
    let guard = 0;
    while (ecgAccum >= step && guard++ < 60) {
      heart.advance(step);
      heart.leads();        // updates the current VCG
      pushVcg(heart.vcg);
      ecgAccum -= step;
    }
    state.phase = heart.mechanical().phase;
    const scrub = document.getElementById('scrub');
    scrub.value = Math.round(state.phase * 1000);
    scrub.setAttribute('aria-valuetext', Math.round(state.phase * 100) + '% — ' + phaseName(state.phase));
  }

  const G = geometryAt(state.phase, state.path);
  const probe = effectiveProbe(G.A);

  echo.depthCm = state.depthCm || 17;
  heart3d.updateProbe(probe, echo.sectorHalf, echo.depthCm);
  heart3d.update(G, state.path, probe, dt || 0.016);
  echo.render(G, state.path, probe, now * 0.001);

  drawEcg();
  drawTwelveLead();
  // refresh the scanner banner's frame-rate readout a couple of times a second —
  // often enough to be live, rarely enough not to flicker or cost layout work
  if (now - lastBannerT > 500) { lastBannerT = now; updateProbeInfo(!!TEE_VIEWS[state.viewName]); }
  document.getElementById('phaseName').textContent = phaseName(state.phase);
  publishMetrics();
  drawPvLoop(G);
  publishHemo(G);
  publishImageQuality();
  measure.draw();
  drawBullseye();

  requestAnimationFrame(loop);
}

// Peak velocity / gradient upper limits (not part of the REF chamber ranges).
const VEL_MAX = 2.5;   // m/s
const GRAD_MAX = 25;   // mmHg

// Render the displayed "normal band" text straight from REF so the shown range
// always matches the thresholds used to flag abnormal values (one source).
function applyRefText() {
  const cm = (r) => r[0].toFixed(1) + '–' + r[1].toFixed(1) + ' cm';
  setText('refLvidd', 'normal ' + cm(REF.lviddNormal));
  setText('refLvids', 'normal ' + cm(REF.lvidsNormal));
  setText('refEf', 'normal ≥' + REF.efNormal[0] + '%');
}

// mirror echo.metrics into the Measurements panel
const fmt = (v, d, u) => (v == null ? '—' : v.toFixed(d) + u);
function publishMetrics() {
  const m = echo.metrics || {};
  setText('mLvidd', fmt(m.lvidd, 1, ' cm'));
  setText('mLvids', fmt(m.lvids, 1, ' cm'));
  setText('mEf', fmt(m.ef, 0, ' %'));
  setText('mMapse', fmt(m.mapse, 0, ' mm'));
  setText('mVel', fmt(m.peakVel, 1, ' m/s'));
  // show the jet/flow label (e.g. 'MR jet') next to Peak V when the echo engine provides one
  setText('mVelLabel', m.peakLabel ? '(' + m.peakLabel + ')' : '');
  setText('mGrad', fmt(m.peakGrad, 0, ' mmHg'));
  setText('mEff', m.effusion == null ? '—' : m.effusion.toFixed(1) + ' cm');
  // flag every value that falls outside its normal reference range, so the
  // learner sees exactly which measurement is abnormal (not just the EF).
  const warn = (id, bad) => { const el = document.getElementById(id); if (el) el.classList.toggle('warn', !!bad); };
  // thresholds driven from the REF normal ranges so there's a single source
  warn('mLvidd', m.lvidd != null && m.lvidd > REF.lviddNormal[1]);
  warn('mLvids', m.lvids != null && m.lvids > REF.lvidsNormal[1]);
  warn('mEf', m.ef != null && m.ef < REF.efNormal[0]);
  warn('mMapse', m.mapse != null && m.mapse < 10);
  warn('mVel', m.peakVel != null && m.peakVel > VEL_MAX);
  warn('mGrad', m.peakGrad != null && m.peakGrad > GRAD_MAX);
  warn('mEff', m.effusion != null);
  const anyAbn = (m.ef != null && m.ef < REF.efNormal[0]) || m.effusion != null ||
    (m.peakVel != null && m.peakVel > VEL_MAX) || (m.lvidd != null && m.lvidd > REF.lviddNormal[1]);
  const flag = document.getElementById('measFlag');
  if (flag) flag.textContent = anyAbn ? '⚠ abnormal' : '';

  // --- BSE/ASE severity grading badges ---
  // valve/shunt severity comes straight from the echo engine, beside Peak V
  if (m.severity) setBadge('mVel', 'mVelBadge', m.severityLabel || m.severity, sevLevel(m.severity));
  else setBadge('mVel', 'mVelBadge', '', 'normal');
  // DCM: derive an EF impairment band (BSE ≥55 normal, 45–54 mild, 36–44 moderate, ≤35 severe)
  if (m.ef != null && state.path && state.path.dilated) {
    const band = m.ef >= 55 ? ['normal', 'normal'] : m.ef >= 45 ? ['mild', 'mild']
      : m.ef >= 36 ? ['moderate', 'moderate'] : ['severe', 'severe'];
    setBadge('mEf', 'mEfBadge', band[0], band[1]);
  } else setBadge('mEf', 'mEfBadge', '', 'normal');
}
// map an echo-engine severity string to a badge colour class
function sevLevel(sev) {
  if (/very severe|severe|large/.test(sev)) return 'severe';
  if (/moderate/.test(sev)) return 'moderate';
  if (/mild|restrictive|shunt/.test(sev)) return 'mild';
  return 'normal';
}
// create/update/remove a coloured severity badge as a sibling of a metric value
function setBadge(anchorId, badgeId, text, level) {
  const anchor = document.getElementById(anchorId);
  if (!anchor) return;
  let badge = document.getElementById(badgeId);
  if (!text) { if (badge) badge.remove(); return; }
  if (!badge) {
    badge = document.createElement('span');
    badge.id = badgeId;
    anchor.insertAdjacentElement('afterend', badge);
  }
  badge.textContent = text;
  badge.className = 'sev-badge sev-' + level;
}
function setText(id, t) { const el = document.getElementById(id); if (el) el.textContent = t; }

// ---------------------------------------------------------------------------
// Haemodynamics: LV pressure–volume loop + numeric readout.
// The loop shape depends only on the pathology, so it is memoised per pathKey;
// only the moving operating-point marker is redrawn each frame.
// ---------------------------------------------------------------------------
const pvCanvas = document.getElementById('pvLoop');
const pvCtx = pvCanvas ? pvCanvas.getContext('2d') : null;
let _pv = { key: null, pts: null, b: null };

function buildPvCache(path, key) {
  const pts = pvLoop(path, 160);
  let vmin = 1e9, vmax = -1e9, pmin = 1e9, pmax = -1e9;
  for (const p of pts) {
    if (p.V < vmin) vmin = p.V; if (p.V > vmax) vmax = p.V;
    if (p.P < pmin) pmin = p.P; if (p.P > pmax) pmax = p.P;
  }
  const vpad = (vmax - vmin) * 0.14 + 3, ppad = (pmax - pmin) * 0.12 + 4;
  _pv = { key, pts, b: { vmin: vmin - vpad, vmax: vmax + vpad, pmin: Math.min(0, pmin - ppad), pmax: pmax + ppad } };
}

function drawPvLoop(G) {
  if (!pvCtx) return;
  if (_pv.key !== state.pathKey || !_pv.pts) buildPvCache(state.path, state.pathKey);
  const b = _pv.b, W = pvCanvas.width, H = pvCanvas.height;
  const padL = 34, padB = 22, padT = 8, padR = 8;
  const x = (V) => padL + (V - b.vmin) / (b.vmax - b.vmin) * (W - padL - padR);
  const y = (P) => H - padB - (P - b.pmin) / (b.pmax - b.pmin) * (H - padT - padB);
  const css = getComputedStyle(document.documentElement);
  const ink = css.getPropertyValue('--ink').trim() || '#e8eef6';
  const dim = css.getPropertyValue('--muted').trim() || '#8b98a8';
  pvCtx.clearRect(0, 0, W, H);
  // axes
  pvCtx.strokeStyle = 'rgba(139,152,168,0.35)'; pvCtx.lineWidth = 1;
  pvCtx.beginPath(); pvCtx.moveTo(padL, padT); pvCtx.lineTo(padL, H - padB); pvCtx.lineTo(W - padR, H - padB); pvCtx.stroke();
  pvCtx.fillStyle = dim; pvCtx.font = '9px system-ui, sans-serif';
  pvCtx.textAlign = 'center'; pvCtx.fillText('LV volume (mL)', (padL + W) / 2, H - 4);
  pvCtx.save(); pvCtx.translate(9, (padT + H) / 2); pvCtx.rotate(-Math.PI / 2);
  pvCtx.fillText('LV pressure (mmHg)', 0, 0); pvCtx.restore();
  // the loop
  pvCtx.strokeStyle = '#c0392b'; pvCtx.lineWidth = 1.8; pvCtx.beginPath();
  _pv.pts.forEach((p, i) => { const px = x(p.V), py = y(p.P); i ? pvCtx.lineTo(px, py) : pvCtx.moveTo(px, py); });
  pvCtx.closePath(); pvCtx.stroke();
  // operating-point marker at the current phase (from the live per-frame snapshot)
  const hd = G && G.hemo;
  if (hd && hd.Vlv != null && hd.Plv != null) {
    pvCtx.fillStyle = '#39c0e8';
    pvCtx.beginPath(); pvCtx.arc(x(hd.Vlv), y(hd.Plv), 3.4, 0, Math.PI * 2); pvCtx.fill();
    pvCtx.strokeStyle = 'rgba(255,255,255,0.8)'; pvCtx.lineWidth = 1; pvCtx.stroke();
  }
}

// numeric haemodynamics readout (EDV/ESV/SV/EF/pressures/gradient/regurg fraction)
let _hemoKey = null, _hemoSum = null;
function publishHemo(G) {
  if (state.pathKey !== _hemoKey || !_hemoSum) { _hemoSum = hemoSummary(state.path); _hemoKey = state.pathKey; }
  const s = _hemoSum;
  const n = (v, d, u) => (v == null || Number.isNaN(v)) ? '—' : v.toFixed(d) + (u || '');
  setText('hEdv', n(s.EDV, 0, ' mL'));
  setText('hEsv', n(s.ESV, 0, ' mL'));
  setText('hSv', n(s.forwardSV != null ? s.forwardSV : s.SV, 0, ' mL'));
  setText('hEf', n(s.EF, 0, ' %'));
  setText('hPlv', n(s.PlvSys, 0, ' mmHg'));
  setText('hPao', s.PaoSys != null ? (Math.round(s.PaoSys) + '/' + Math.round(s.PaoDia)) : '—');
  // The peak instantaneous LV–aortic gradient is only a clinically meaningful
  // headline in aortic stenosis; across a normal valve it is trivial. Show the
  // modelled value only for AS so the readout can't be mistaken for a stenosis.
  setText('hGrad', state.path && state.path.aorticStenosis ? n(s.gradient, 0, ' mmHg') : 'trivial');
  // estimated PA systolic pressure from the TR jet (echo metric), when a TR jet exists
  const pasp = echo.metrics && echo.metrics.pasp;
  setText('hPasp', pasp != null ? Math.round(pasp) + ' mmHg' : '—');
  const pEl = document.getElementById('hPasp');
  if (pEl) pEl.classList.toggle('warn', pasp != null && pasp > 35); // PH threshold ~>35
  // mitral inflow E/A (diastolic function) — mitral-stenosis inflow is gated differently
  const ea = eaRatio(state.path);
  setText('hEa', ea != null ? ea.toFixed(2) : '—');
  const eEl = document.getElementById('hEa');
  if (eEl) eEl.classList.toggle('warn', ea != null && ea < 0.8); // impaired relaxation
  setText('hRf', (s.regurgFraction != null && s.regurgFraction > 0.01) ? Math.round(s.regurgFraction * 100) + ' %' : '—');
  // global longitudinal strain (memoised per pathology in the model); flag when
  // impaired (ASE: normal more negative than about −18%).
  const gls = longitudinalStrain(state.path);
  setText('hGls', gls == null ? '—' : gls.toFixed(1) + ' %');
  const gEl = document.getElementById('hGls');
  if (gEl) gEl.classList.toggle('warn', gls != null && gls > -16);
}

// live image-quality self-validation readout (gCNR / CNR / speckle SNR) from echo
function publishImageQuality() {
  const m = echo.metrics || {};
  const n = (v, d) => (v == null || Number.isNaN(v)) ? '—' : v.toFixed(d);
  setText('qGcnr', n(m.gcnr, 2));
  setText('qCnr', n(m.cnr, 1));
  setText('qSnr', n(m.speckleSNR, 2));
  const g = document.getElementById('qGcnr');
  if (g) g.classList.toggle('warn', m.gcnr != null && m.gcnr < 0.6);
}

// ---------------------------------------------------------------------------
// Regional-strain bullseye (AHA 17-segment). Depends only on the pathology, so
// it is redrawn when the pathology key changes, not every frame.
// ---------------------------------------------------------------------------
const bullCanvas = document.getElementById('bullseye');
const bullCtx = bullCanvas ? bullCanvas.getContext('2d') : null;
// diverging strain colormap: deep red = normal (≈ −20 %), through pink, to blue (≈ 0)
function strainColor(s) {
  const t = Math.max(0, Math.min(1, (-s) / 22)); // 0 (blue) … 1 (deep red)
  const stops = [
    [0.00, [43, 108, 176]],   // blue  (0 %, akinetic/dyskinetic)
    [0.30, [150, 170, 210]],  // pale
    [0.45, [233, 150, 150]],  // pink
    [0.70, [192, 57, 43]],    // red
    [1.00, [122, 21, 18]],    // deep red (normal shortening)
  ];
  let a = stops[0], b = stops[stops.length - 1];
  for (let i = 0; i < stops.length - 1; i++) { if (t >= stops[i][0] && t <= stops[i + 1][0]) { a = stops[i]; b = stops[i + 1]; break; } }
  const f = (t - a[0]) / (b[0] - a[0] || 1);
  const c = a[1].map((v, i) => Math.round(v + (b[1][i] - v) * f));
  return `rgb(${c[0]},${c[1]},${c[2]})`;
}
let _bullKey = null;
function drawBullseye() {
  if (!bullCtx) return;
  if (state.pathKey === _bullKey) return; // only redraw on pathology change
  _bullKey = state.pathKey;
  const { segments, gls } = regionalStrain(state.path);
  const W = bullCanvas.width, H = bullCanvas.height, cx = W / 2, cy = H / 2, R = Math.min(W, H) / 2 - 6;
  bullCtx.clearRect(0, 0, W, H);
  // ring radii: apex centre, then apical(4), mid(6), basal(6)
  const rings = [[0.74, 1.0, 0], [0.48, 0.74, 1], [0.22, 0.48, 2]]; // [inner,outer,level]
  const wedge = (a0, a1, r0, r1, col) => {
    bullCtx.beginPath();
    bullCtx.arc(cx, cy, r1 * R, a0, a1);
    bullCtx.arc(cx, cy, r0 * R, a1, a0, true);
    bullCtx.closePath();
    bullCtx.fillStyle = col; bullCtx.fill();
    bullCtx.strokeStyle = 'rgba(8,11,17,0.9)'; bullCtx.lineWidth = 1.5; bullCtx.stroke();
  };
  for (const [ri, ro, level] of rings) {
    const segs = segments.filter((s) => s.level === level);
    const span = (Math.PI * 2) / segs.length;
    for (const s of segs) {
      const ca = -s.angle;                 // model→canvas (anterior up, lateral right)
      wedge(ca - span / 2, ca + span / 2, ri, ro, strainColor(s.strain));
    }
  }
  // apex centre disc
  const apex = segments.find((s) => s.level === 3);
  bullCtx.beginPath(); bullCtx.arc(cx, cy, 0.22 * R, 0, Math.PI * 2);
  bullCtx.fillStyle = strainColor(apex ? apex.strain : -20); bullCtx.fill();
  bullCtx.strokeStyle = 'rgba(8,11,17,0.9)'; bullCtx.lineWidth = 1.5; bullCtx.stroke();
  // centre GLS label
  bullCtx.fillStyle = '#fff'; bullCtx.font = 'bold 15px system-ui, sans-serif';
  bullCtx.textAlign = 'center'; bullCtx.textBaseline = 'middle';
  bullCtx.fillText(Math.round(gls) + '%', cx, cy);
  // orientation cue
  bullCtx.fillStyle = 'rgba(200,215,230,0.75)'; bullCtx.font = '9px system-ui, sans-serif';
  bullCtx.fillText('ANT', cx, cy - R + 6); bullCtx.fillText('INF', cx, cy + R - 6);
  bullCtx.textAlign = 'left'; bullCtx.fillText('SEPT', 3, cy); bullCtx.textAlign = 'right'; bullCtx.fillText('LAT', W - 3, cy);
  // readouts
  setText('sGls', Math.round(gls) + ' %');
  const worst = segments.reduce((m, s) => (s.strain > m.strain ? s : m), segments[0]);
  const abn = worst.strain > -12;
  setText('sPattern', !abn ? 'normal / uniform' : (state.path.dilated ? 'global impairment' : worst.name + ' loss'));
  const gEl = document.getElementById('sGls'); if (gEl) gEl.classList.toggle('warn', gls > -16);
}

// learner-measurement readout + live accuracy vs the model's own ground truth
function publishMeasure(r) {
  setText('measDist', r.distCm != null ? r.distCm.toFixed(1) + ' cm' : '—');
  const tag = document.getElementById('measSimpTag');
  const efEl = document.getElementById('measSimpEf');
  const errEl = document.getElementById('measErr');
  const accRow = document.getElementById('measAcc');
  if (r.ef != null) {
    efEl.textContent = Math.round(r.ef) + ' %';
    if (tag) tag.textContent = `(EDV ${Math.round(r.edv)} / ESV ${Math.round(r.esv)} mL)`;
    // compare learner EF to the model ground-truth EF (echo.metrics.ef)
    const truth = echo.metrics && echo.metrics.ef;
    if (truth != null && errEl) {
      const d = r.ef - truth;
      errEl.textContent = (d >= 0 ? '+' : '') + d.toFixed(0) + ' % (truth ' + Math.round(truth) + ')';
      // ASE test–retest for biplane EF is ~±5–6 %; flag good vs off
      const good = Math.abs(d) <= 6;
      errEl.classList.toggle('good', good);
      errEl.classList.toggle('warn', !good);
      if (accRow) accRow.style.display = '';
    }
  } else {
    efEl.textContent = r.which === 'ED only' ? 'trace ES…' : r.which === 'ES only' ? 'trace ED…' : '—';
    if (tag) tag.textContent = '';
    if (errEl) { errEl.textContent = '—'; errEl.classList.remove('good', 'warn'); }
  }
}

// ---------------------------------------------------------------------------
// UI wiring
// ---------------------------------------------------------------------------
function setView(name) {
  state.viewName = name;
  if (echo) echo.viewName = name; // let the echo label whitelist know the current view
  const tee = TEE_VIEWS[name];
  if (tee) {
    state.baseProbe = tee.probe();
    state.depthCm = tee.depth; echo.depthCm = tee.depth;
    echo.freqMHz = 5.5;                 // TEE probes run higher-frequency → finer PSF
    echo.nearFieldCm = 0.3;             // only the oesophageal wall lies in front
    // small high-frequency TEE aperture: a thin slab focused in the near field
    echo.elevFocus = 4.0; echo.elevMin = 0.05; echo.elevDiv = 0.015;
    setRange('depth', tee.depth, 'depthOut', String(tee.depth), tee.depth + ' cm');
    heart3d.setCameraForView(TEE_CAM[name] || 'A4C');
  } else {
    state.baseProbe = (VIEWS[name] || EXTRA_VIEWS[name].probe)();
    echo.freqMHz = 2.7;                 // standard transthoracic probe
    echo.elevFocus = 8.0; echo.elevMin = 0.09; echo.elevDiv = 0.025;
    // chest wall (parasternal / apical) or abdominal wall (subcostal) under the probe
    echo.nearFieldCm = name === 'SUBCOSTAL' ? 1.6 : /^A\dC$/.test(name) ? 1.4 : 2.0;
    const d = viewDepth(name);          // fill the FOV per view
    state.depthCm = d; echo.depthCm = d;
    setRange('depth', d, 'depthOut', String(d), d + ' cm');
    heart3d.setCameraForView(name);
  }
  state.az = 0; state.tilt = 0; state.slide = 0;
  setRange('az', 0, 'azOut', '0°', '0°');
  setRange('tilt', 0, 'tiltOut', '0°', '0°');
  setRange('slide', 0, 'slideOut', '0', '0 cm');
  document.getElementById('viewTag').textContent = VIEW_LABEL[name] || name;
  document.querySelectorAll('[data-view]').forEach((b) => {
    const on = b.dataset.view === name;
    b.classList.toggle('active', on);
    b.setAttribute('aria-pressed', on ? 'true' : 'false');
  });
  updateProbeInfo(tee);
  updatePlaneFlag();
}

// Scanner-style probe banner. Real machines annotate the active imaging mode on
// screen, so the harmonic state belongs here rather than only in the console.
function updateProbeInfo(tee) {
  const info = document.getElementById('probeInfo');
  if (!info) return;
  const base = tee ? '5.5 MHz · TEE' : '2.5 MHz · sector';
  const fps = echo.metrics && echo.metrics.fps;
  const rate = fps ? ' · ' + Math.round(fps) + ' Hz' : '';
  info.textContent = (echo.harmonic ? base + ' · THI' : base) + rate;
}

document.querySelectorAll('[data-view]').forEach((b) =>
  b.addEventListener('click', () => setView(b.dataset.view)));

// TTE / TEE modality switch — swaps the standard-view button row and jumps to that
// modality's default window.
document.querySelectorAll('#modalityCtl .seg-btn').forEach((b) => b.addEventListener('click', () => {
  if (challenge.active) return;
  const mod = b.dataset.modality;
  document.querySelectorAll('#modalityCtl .seg-btn').forEach((x) => {
    const on = x.dataset.modality === mod; x.classList.toggle('active', on); x.setAttribute('aria-pressed', String(on));
  });
  const tee = mod === 'TEE';
  document.getElementById('viewBtns').classList.toggle('hidden', tee);
  document.getElementById('teeBtns').classList.toggle('hidden', !tee);
  // restore the standard TTE per-view depth (setView will re-apply on the next view)
  if (!tee) { const d = viewDepth(state.viewName); state.depthCm = d; echo.depthCm = d; setRange('depth', d, 'depthOut', String(d), d + ' cm'); }
  setView(tee ? 'ME4C' : 'PLAX');
}));

function renderTeaching(key) {
  if (challenge.active) return; // challenge mode hides the answer
  const p = PATHOLOGY[key], c = CHECKS[key] || { look: [], hint: '' };
  // render each item as a tickable checkbox the learner can mark off (client-side only)
  const list = c.look.map((x) =>
    '<li><label class="cl-item"><input type="checkbox" />' + x + '</label></li>').join('');
  document.getElementById('teaching').innerHTML =
    '<b>' + p.name + '</b><br>' + p.text +
    '<div class="checklist"><span class="cl-title">What to look for</span><ul>' + list + '</ul>' +
    '<span class="cl-hint">💡 ' + c.hint + '</span></div>';
}

// pathologies whose severity can be graded mild/moderate/severe
const GRADEABLE = new Set(['as', 'mr', 'ms', 'tr']);
const severityCtl = document.getElementById('severityCtl');
const gradeBtns = severityCtl ? [...severityCtl.querySelectorAll('.seg-btn')] : [];
function setGradeButtons(grade) {
  gradeBtns.forEach((x) => { const on = x.dataset.grade === grade; x.classList.toggle('active', on); x.setAttribute('aria-pressed', String(on)); });
}
// show/hide + reset the severity control for the current pathology key
function syncSeverityCtl(key) {
  const on = GRADEABLE.has(key);
  if (severityCtl) severityCtl.classList.toggle('hidden', !on);
}

document.getElementById('pathology').addEventListener('change', (e) => {
  if (challenge.active) return;
  state.pathKey = e.target.value;
  state.path = { ...PATHOLOGY[e.target.value].flags, grade: 'severe' };
  setGradeButtons('severe');
  syncSeverityCtl(e.target.value);
  heart3d.setPathology(state.path);
  heart.setPathology(state.path);   // ECG morphology follows the pathology (LVH, low-voltage, …)
  renderTeaching(e.target.value);
  setView(PATHOLOGY[e.target.value].view);
});

// severity grade — reapply a fresh path object so the echo re-measures (it memoises
// by path identity) and the graded orifice feeds the continuity-Doppler + PV loop.
gradeBtns.forEach((b) => b.addEventListener('click', () => {
  if (challenge.active) return;
  const grade = b.dataset.grade;
  state.path = { ...state.path, grade };
  setGradeButtons(grade);
  heart.setPathology(state.path);
  announce(`${grade} severity`);
}));

const playBtn = document.getElementById('playBtn');
playBtn.addEventListener('click', () => {
  state.playing = !state.playing;
  playBtn.textContent = state.playing ? '⏸ Pause' : '▶ Play';
});

// --- caliper / Simpson-EF measurement tools ---
const measBtns = document.querySelectorAll('.meas-tools .seg-btn');
measBtns.forEach((b) => b.addEventListener('click', () => {
  const active = measure.setTool(b.dataset.tool);
  measBtns.forEach((x) => { const on = x.dataset.tool === active; x.classList.toggle('active', on); x.setAttribute('aria-pressed', String(on)); });
  // tracing an endocardial border is only meaningful on a still frame — auto-freeze
  if (active === 'simpsonED' || active === 'simpsonES') {
    state.playing = false; playBtn.textContent = '▶ Play';
  }
}));
document.getElementById('measClear').addEventListener('click', () => {
  measure.clear();
  measBtns.forEach((x) => { x.classList.remove('active'); x.setAttribute('aria-pressed', 'false'); });
  measure.setTool('off');
});
document.getElementById('measFreeze').addEventListener('click', () => {
  state.playing = !state.playing;
  playBtn.textContent = state.playing ? '⏸ Pause' : '▶ Play';
});

document.getElementById('hr').addEventListener('input', (e) => {
  state.hr = +e.target.value;
  heart.setHR(state.hr);            // drives the rhythm engine's R-R interval
  document.getElementById('hrOut').textContent = e.target.value;
  document.getElementById('hrReadout').textContent = e.target.value + ' bpm';
  e.target.setAttribute('aria-valuetext', e.target.value + ' bpm');
});
// Rhythm selector (sinus / AF / brady / tachy / AV block / PVC / LBBB / RBBB / STEMI / hyperK)
const rhythmSel = document.getElementById('rhythm');
if (rhythmSel) rhythmSel.addEventListener('change', (e) => {
  heart.setRhythm(e.target.value);
  state.playing = true; playBtn.textContent = '⏸ Pause';
});
// Lead selector for the top rhythm strip
const leadSel = document.getElementById('ecgLead');
if (leadSel) leadSel.addEventListener('change', (e) => { ecgLead = e.target.value; });
// 12-lead modal open/close
const twelveModal = document.getElementById('twelveModal');
const twelveBtn = document.getElementById('twelveBtn');
if (twelveBtn) twelveBtn.addEventListener('click', () => {
  twelveModal.classList.remove('hidden');
  const cv = document.getElementById('twelveLead');
  cv.width = cv.clientWidth; cv.height = cv.clientHeight;
});
const twelveClose = document.getElementById('twelveClose');
if (twelveClose) twelveClose.addEventListener('click', () => twelveModal.classList.add('hidden'));
if (twelveModal) twelveModal.addEventListener('click', (e) => { if (e.target === twelveModal) twelveModal.classList.add('hidden'); });
document.getElementById('scrub').addEventListener('input', (e) => {
  state.playing = false;
  playBtn.textContent = '▶ Play';
  state.phase = +e.target.value / 1000;
  e.target.setAttribute('aria-valuetext', Math.round(+e.target.value / 10) + '%');
});

// set a range slider's value, its <output> text and its aria-valuetext together
function setRange(id, value, outId, outText, valuetext) {
  const el = document.getElementById(id);
  if (el) { el.value = value; if (valuetext != null) el.setAttribute('aria-valuetext', valuetext); }
  if (outId) { const o = document.getElementById(outId); if (o) o.textContent = outText; }
}
const bind = (id, outId, fn, suffix = '', unit) => {
  document.getElementById(id).addEventListener('input', (e) => {
    const v = +e.target.value;
    fn(v);
    const disp = (suffix === '°' ? v : e.target.value) + suffix;
    if (outId) document.getElementById(outId).textContent = disp;
    // screen readers announce the value with its unit
    e.target.setAttribute('aria-valuetext', unit != null ? e.target.value + unit : disp);
  });
};
function updatePlaneFlag() {
  const modified = !!(state.az || state.tilt || state.slide);
  const f = document.getElementById('planeFlag');
  if (f) f.textContent = modified ? '• modified' : '';
  // the reset affordance reflects whether the plane has actually been moved
  const reset = document.getElementById('resetPlane');
  if (reset) { reset.disabled = !modified; reset.classList.toggle('muted', !modified); }
}
bind('az', 'azOut', (v) => { state.az = v; updatePlaneFlag(); }, '°', '°');
bind('tilt', 'tiltOut', (v) => { state.tilt = v; updatePlaneFlag(); }, '°', '°');
bind('slide', 'slideOut', (v) => { state.slide = v; updatePlaneFlag(); }, '', ' cm');
bind('depth', 'depthOut', (v) => (state.depthCm = v), '', ' cm');
bind('gain', 'gainOut', (v) => (echo.gain = v));
bind('nyq', 'nyqOut', (v) => (echo.nyquist = v), '', ' m/s');
bind('pers', 'persOut', (v) => (echo.persistence = v));
// Focus: 0 means "auto" (the renderer tracks mid-field), any other value pins the
// transmit focal depth in cm. Handled directly rather than via bind() so the zero
// position can read as "auto" instead of "0 cm".
document.getElementById('focus').addEventListener('input', (e) => {
  const v = +e.target.value;
  echo.focusCm = v > 0 ? v : null;
  const txt = v > 0 ? v.toFixed(1) + ' cm' : 'auto';
  document.getElementById('focusOut').textContent = txt;
  e.target.setAttribute('aria-valuetext', txt);
});
// Detail: 0 = Auto (FPS-adaptive), 1/2/3 = pinned sampling stride (resolution).
document.getElementById('detail').addEventListener('change', (e) => {
  const v = +e.target.value;
  echo._pinStride = v || null;      // 0 → null → adaptive
  if (v) echo._stride = v;          // apply immediately
});

document.getElementById('colorOn').addEventListener('change', (e) => (echo.colorOn = e.target.checked));
document.getElementById('elevOn').addEventListener('change', (e) => {
  echo.elevation = e.target.checked;
  announce(e.target.checked ? 'Slice thickness on' : 'Slice thickness off (infinitely thin plane)');
});
document.getElementById('harmonicOn').addEventListener('change', (e) => {
  echo.harmonic = e.target.checked;
  updateProbeInfo(!!TEE_VIEWS[state.viewName]);
  announce(e.target.checked ? 'Harmonic imaging on' : 'Harmonic imaging off (fundamental)');
});
document.getElementById('labelsOn').addEventListener('change', (e) => (echo.showLabels = e.target.checked));
document.getElementById('flowOn').addEventListener('change', (e) => (heart3d.showFlow = e.target.checked));
document.getElementById('transp').addEventListener('change', (e) => (heart3d.transparent = e.target.checked));
document.getElementById('fibersOn').addEventListener('change', (e) => heart3d.showFibers(e.target.checked));

document.getElementById('resetPlane').addEventListener('click', () => { setView(state.viewName); });

// Reset all — restore gain, nyquist, depth, HR, plane offsets and toggles to defaults
const DEFAULTS = { hr: 72, depth: 15, gain: 1.0, nyq: 0.62 };
function resetAll() {
  // heart rate
  state.hr = DEFAULTS.hr;
  setRange('hr', DEFAULTS.hr, 'hrOut', String(DEFAULTS.hr), DEFAULTS.hr + ' bpm');
  document.getElementById('hrReadout').textContent = DEFAULTS.hr + ' bpm';
  // imaging — depth is view-appropriate so the heart fills the field of view
  const dep = viewDepth(state.viewName);
  state.depthCm = dep; echo.depthCm = dep;
  setRange('depth', dep, 'depthOut', String(dep), dep + ' cm');
  echo.gain = DEFAULTS.gain;
  setRange('gain', DEFAULTS.gain, 'gainOut', DEFAULTS.gain.toFixed(1), DEFAULTS.gain.toFixed(1));
  echo.nyquist = DEFAULTS.nyq;
  setRange('nyq', DEFAULTS.nyq, 'nyqOut', DEFAULTS.nyq.toFixed(2), DEFAULTS.nyq + ' m/s');
  // toggles back on
  // harmonic imaging is off by default (fundamental), unlike the toggles below
  const harm = document.getElementById('harmonicOn');
  if (harm) { harm.checked = false; echo.harmonic = false; }
  echo.focusCm = null;   // back to auto (mid-field) focus
  setRange('focus', 0, 'focusOut', 'auto', 'auto');
  echo.persistence = 0;
  setRange('pers', 0, 'persOut', '0', '0');
  const toggles = [['elevOn', (v) => (echo.elevation = v)],
    ['colorOn', (v) => (echo.colorOn = v)], ['labelsOn', (v) => (echo.showLabels = v)],
    ['flowOn', (v) => (heart3d.showFlow = v)], ['transp', (v) => (heart3d.transparent = v)]];
  toggles.forEach(([id, fn]) => { const c = document.getElementById(id); if (c) { c.checked = true; fn(true); } });
  // plane offsets (az/tilt/slide) reset via setView, which also clears the flag
  setView(state.viewName);
  announce('Controls reset to defaults');
}
// announce a transient message through the polite live region for screen readers
function announce(msg) {
  const el = document.getElementById('liveStatus');
  if (el) el.textContent = msg;
}
document.getElementById('resetAll').addEventListener('click', resetAll);

// ---- help / onboarding overlay ----
const helpOverlay = document.getElementById('helpOverlay');
const helpBtn = document.getElementById('helpBtn');
const helpCard = helpOverlay.querySelector('.overlay-card');
let helpOpener = null;
const helpOpen = () => !helpOverlay.classList.contains('hidden');
const showHelp = () => {
  helpOpener = document.activeElement;
  helpOverlay.classList.remove('hidden');
  // move focus into the dialog for keyboard users
  const first = document.getElementById('helpGot');
  if (first) first.focus();
};
const hideHelp = () => {
  helpOverlay.classList.add('hidden');
  // restore focus to whatever opened the dialog (usually the help button)
  const back = helpOpener || helpBtn;
  if (back && typeof back.focus === 'function') back.focus();
  helpOpener = null;
};
helpBtn.addEventListener('click', showHelp);
document.getElementById('helpClose').addEventListener('click', hideHelp);
document.getElementById('helpGot').addEventListener('click', () => { hideHelp(); try { localStorage.setItem('echosim_seen', '1'); } catch (e) {} });
// click on the backdrop (but not the card) dismisses
helpOverlay.addEventListener('click', (e) => { if (e.target === helpOverlay) hideHelp(); });
// Escape closes; Tab is trapped within the card while open
helpOverlay.addEventListener('keydown', (e) => {
  if (!helpOpen()) return;
  if (e.key === 'Escape') { e.preventDefault(); hideHelp(); return; }
  if (e.key === 'Tab') {
    const focusables = helpCard.querySelectorAll('button, [href], input, [tabindex]:not([tabindex="-1"])');
    if (!focusables.length) return;
    const first = focusables[0], last = focusables[focusables.length - 1];
    if (e.shiftKey && document.activeElement === first) { e.preventDefault(); last.focus(); }
    else if (!e.shiftKey && document.activeElement === last) { e.preventDefault(); first.focus(); }
  }
});
try { if (!localStorage.getItem('echosim_seen')) showHelp(); } catch (e) { showHelp(); }

// ---- keyboard shortcuts ----
window.addEventListener('keydown', (e) => {
  if (/INPUT|SELECT|TEXTAREA/.test(document.activeElement?.tagName || '')) return;
  if (helpOpen()) return; // don't let shortcuts fire behind the open help modal
  if (e.key === ' ') { e.preventDefault(); playBtn.click(); }
  else if (e.key === 'ArrowLeft') { state.playing = false; playBtn.textContent = '▶ Play'; state.phase = (state.phase - 0.02 + 1) % 1; }
  else if (e.key === 'ArrowRight') { state.playing = false; playBtn.textContent = '▶ Play'; state.phase = (state.phase + 0.02) % 1; }
  else if (e.key >= '1' && e.key <= '5') { const v = Object.keys(VIEWS)[+e.key - 1]; if (v) setView(v); }
  else if (e.key === 'r' || e.key === 'R') { setView(state.viewName); }
  else if (e.key === 'l' || e.key === 'L') { const c = document.getElementById('labelsOn'); c.checked = !c.checked; echo.showLabels = c.checked; }
  else if (e.key === 'c' || e.key === 'C') { document.getElementById('challengeBtn').click(); }
});

// ---- challenge / quiz mode ----
// One-line "signature" of each finding, used to contrast a wrong guess with the
// correct answer in challenge feedback.
const SIGNATURE = {
  normal:   'laminar flow with clean valves and no jet',
  mr:       'a SYSTOLIC regurgitant jet INTO the left atrium',
  ms:       'restricted, accelerated DIASTOLIC mitral inflow',
  as:       'a SYSTOLIC ejection jet OUT through a restricted aortic valve',
  dcm:      'a dilated, poorly-contracting LV with a low EF',
  rwma:     'a single akinetic wall segment with normal thickening elsewhere',
  effusion: 'an echo-free fluid rim surrounding the heart',
  vsd:      'a systolic shunt jet crossing the VENTRICULAR septum',
  asd:      'a low-velocity shunt across the ATRIAL septum',
  tr:       'a systolic regurgitant jet INTO the right atrium',
  phtn:     'a dilated pressure-loaded RV flattening the septum (D-sign) with a TR jet',
};
const challenge = { active: false, answer: null, score: 0, total: 0, difficulty: 'easy', selectedView: null };
const CHOICES = Object.keys(PATHOLOGY);
const VIEW_KEYS = Object.keys(VIEWS);
const challengeBtn = document.getElementById('challengeBtn');

// show/hide the view-name hint tag on the Echo Display header
function setViewTagHidden(hidden) {
  const t = document.getElementById('viewTag');
  if (t) t.style.visibility = hidden ? 'hidden' : '';
}

// difficulty segmented control
document.querySelectorAll('#difficultyCtl .seg-btn').forEach((b) => b.addEventListener('click', () => {
  challenge.difficulty = b.dataset.diff;
  document.querySelectorAll('#difficultyCtl .seg-btn').forEach((x) => {
    const on = x === b;
    x.classList.toggle('active', on);
    x.setAttribute('aria-pressed', on ? 'true' : 'false');
  });
  if (challenge.active) newCase(); // re-deal the current case under the new tier
}));
function newCase() {
  const diff = challenge.difficulty;
  const key = CHOICES[Math.floor(pseudoRand() * CHOICES.length)];
  challenge.answer = key;
  challenge.selectedView = null;
  state.pathKey = key;
  state.path = { ...PATHOLOGY[key].flags };
  heart3d.setPathology(state.path);
  setView(PATHOLOGY[key].view);

  // Easy = view preset + finding pre-centred (view tag visible).
  // Medium/Hard = hide the view-name hint and start the plane slightly off-axis
  // so the learner must adjust to bring the finding into the beam.
  setViewTagHidden(diff !== 'easy');
  // Hard = blank the (disabled) pathology dropdown so no label leaks the answer.
  document.getElementById('pathology').selectedIndex = (diff === 'hard') ? -1 : 0;
  if (diff !== 'easy') {
    const az = Math.round((pseudoRand() * 2 - 1) * 22);
    const tilt = Math.round((pseudoRand() * 2 - 1) * 12);
    state.az = az; state.tilt = tilt;
    setRange('az', az, 'azOut', az + '°', az + '°');
    setRange('tilt', tilt, 'tiltOut', tilt + '°', tilt + '°');
    updatePlaneFlag();
  }

  // build 4 finding options incl. the answer
  const opts = new Set([key]);
  while (opts.size < 4) opts.add(CHOICES[Math.floor(pseudoRand() * CHOICES.length)]);
  const shuffled = [...opts].sort(() => pseudoRand() - 0.5);
  const findingBtns = shuffled.map((k) => `<button class="quiz-opt" data-k="${k}">${PATHOLOGY[k].name}</button>`).join('');

  // Hard also requires naming the view
  let viewBlock = '';
  if (diff === 'hard') {
    const viewBtns = VIEW_KEYS.map((v) => `<button class="quiz-view" data-v="${v}">${v}</button>`).join('');
    viewBlock = `<div class="quiz-sub">1 · Standard view</div><div class="quiz-views">${viewBtns}</div>` +
      `<div class="quiz-sub">2 · Finding</div>`;
  }
  const prompt = diff === 'hard' ? 'Name the view and the finding' : 'Identify the finding';
  document.getElementById('teaching').innerHTML =
    `<div class="quiz"><b>${prompt}</b><div class="quiz-score">Score ${challenge.score}/${challenge.total}</div>` +
    viewBlock + `<div class="quiz-opts">${findingBtns}</div><div id="quizFeedback" role="status" aria-live="polite"></div></div>`;
  document.querySelectorAll('.quiz-opt').forEach((b) => b.addEventListener('click', () => answerCase(b.dataset.k)));
  document.querySelectorAll('.quiz-view').forEach((b) => b.addEventListener('click', () => {
    challenge.selectedView = b.dataset.v;
    document.querySelectorAll('.quiz-view').forEach((x) => x.classList.toggle('sel', x === b));
  }));
}
function answerCase(k) {
  const hard = challenge.difficulty === 'hard';
  // Hard: a view must be named first (both answers required)
  if (hard && !challenge.selectedView) {
    const fb = document.getElementById('quizFeedback');
    if (fb) fb.innerHTML = `<div class="quiz-fb no">Pick the standard view first, then the finding.</div>`;
    return;
  }
  challenge.total++;
  const correctFinding = k === challenge.answer;
  const correctViewName = PATHOLOGY[challenge.answer].view;
  const correctView = !hard || challenge.selectedView === correctViewName;
  const correct = correctFinding && correctView;
  if (correct) challenge.score++;
  const p = PATHOLOGY[challenge.answer];
  // on a wrong finding, contrast the chosen pathology with the correct one
  let discriminator = '';
  if (!correctFinding && k && k !== challenge.answer) {
    const chosen = PATHOLOGY[k];
    discriminator = `<div class="quiz-why">Why not ${chosen.name}? ${chosen.name} is ${SIGNATURE[k]}, whereas ${p.name} is ${SIGNATURE[challenge.answer]}.</div>`;
  }
  let viewLine = '';
  if (hard) {
    viewLine = correctView
      ? `<div class="quiz-why">View: ✓ ${correctViewName}.</div>`
      : `<div class="quiz-why">View: ✗ you chose ${challenge.selectedView || '—'}; this is <b>${correctViewName}</b>.</div>`;
  }
  document.getElementById('quizFeedback').innerHTML =
    `<div class="quiz-fb ${correct ? 'ok' : 'no'}">${correct ? '✓ Correct' : '✗ Not quite'} — <b>${p.name}</b><br>${p.text}${viewLine}${discriminator}</div>` +
    `<button id="nextCase" class="wide">Next case →</button>`;
  document.querySelector('.quiz-score').textContent = `Score ${challenge.score}/${challenge.total}`;
  document.getElementById('nextCase').addEventListener('click', newCase);
  document.querySelectorAll('.quiz-opt').forEach((b) => { b.disabled = true; if (b.dataset.k === challenge.answer) b.classList.add('correct'); });
  document.querySelectorAll('.quiz-view').forEach((b) => { b.disabled = true; if (b.dataset.v === correctViewName) b.classList.add('correct'); });
}
// tiny deterministic-ish PRNG seeded by phase+total (Math.random unavailable in some sandboxes)
let _seed = 20260707;
function pseudoRand() { _seed = (_seed * 1664525 + 1013904223) & 0x7fffffff; return _seed / 0x7fffffff; }
challengeBtn.addEventListener('click', () => {
  challenge.active = !challenge.active;
  challengeBtn.classList.toggle('active', challenge.active);
  challengeBtn.textContent = challenge.active ? '✕ Exit challenge' : '🎯 Challenge mode';
  document.getElementById('pathology').disabled = challenge.active;
  // hide the severity control during a challenge — its state would leak the lesion
  if (severityCtl) severityCtl.classList.add('hidden');
  if (challenge.active) { challenge.score = 0; challenge.total = 0; _seed = 20260707 + Math.floor(state.phase * 99999); newCase(); }
  else { state.pathKey = 'normal'; state.path = {}; heart3d.setPathology({}); document.getElementById('pathology').value = 'normal'; setViewTagHidden(false); renderTeaching('normal'); setView('PLAX'); }
});

// initialise teaching text + kick off
state.pathKey = 'normal';
state.depthCm = viewDepth('PLAX');
echo.depthCm = state.depthCm;
echo.showLabels = true;
applyRefText();
renderTeaching('normal');
setView('PLAX');
resizeCanvases();

// ---------------------------------------------------------------------------
// Collapsible console panels — the console grew long as analysis panels were
// added; let each section fold so the learner focuses on what they need. State
// persists per-panel in localStorage; secondary/analysis panels start collapsed.
// ---------------------------------------------------------------------------
(function initCollapsiblePanels() {
  const KEY = 'echosim_panels_v1';
  let saved = {};
  try { saved = JSON.parse(localStorage.getItem(KEY) || '{}'); } catch (e) {}
  const defaultCollapsed = new Set(['Caliper tools', 'Image quality', 'Imaging', 'Training']);
  document.querySelectorAll('.console .ctl-group').forEach((grp) => {
    const h3 = grp.querySelector('h3');
    if (!h3) return;
    const key = ((h3.firstChild && h3.firstChild.textContent) || '').trim();
    h3.setAttribute('role', 'button');
    h3.setAttribute('tabindex', '0');
    const collapsed = key in saved ? saved[key] : defaultCollapsed.has(key);
    grp.classList.toggle('collapsed', collapsed);
    h3.setAttribute('aria-expanded', String(!collapsed));
    const toggle = () => {
      const now = grp.classList.toggle('collapsed');
      h3.setAttribute('aria-expanded', String(!now));
      saved[key] = now;
      try { localStorage.setItem(KEY, JSON.stringify(saved)); } catch (e) {}
    };
    h3.addEventListener('click', toggle);
    h3.addEventListener('keydown', (e) => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); toggle(); } });
  });
})();

requestAnimationFrame((t) => { state.lastT = t; loop(t); });
