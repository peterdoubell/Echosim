// Headless B-mode renderer for audits: runs the real EchoView pipeline (sampling,
// specular interfaces, artifacts, PSF, compositing) in node against stub
// canvases, and returns the grey image together with the per-pixel tissue kind,
// depth and beam index, so image-physics properties (which reflector is the
// brightest, how a pleural line fades) can be checked numerically.
import { TTE_VIEWS, EXTRA_VIEWS, TEE_VIEWS } from '../js/views.js';

const noop = () => {};
const stubCtx = () => new Proxy({
  createImageData: (w, h) => ({ width: w, height: h, data: new Uint8ClampedArray(w * h * 4) }),
  measureText: () => ({ width: 0 }),
  getImageData: (x, y, w, h) => ({ width: w, height: h, data: new Uint8ClampedArray(w * h * 4) }),
  createLinearGradient: () => ({ addColorStop: noop }),
}, { get: (t, k) => (k in t ? t[k] : noop), set: () => true });
const stubCanvas = () => ({ width: 800, height: 800, getContext: () => stubCtx(), style: {} });
if (typeof globalThis.document === 'undefined') {
  globalThis.document = { createElement: () => stubCanvas() };
}
const { EchoView } = await import('../js/echo.js');
const { geometryAt } = await import('../js/cardiac-model.js');

const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
let EV = null;
// Render one standard view (name as in views.js) at `phase` for pathology `path`
// with the same per-view settings main.js applies; colour Doppler off.
export function renderBmode(name, phase = 0, path = {}) {
  if (!EV) { EV = new EchoView(stubCanvas(), null); EV._pinStride = 1; EV._stride = 1; }
  const e = EV;
  const tee = TEE_VIEWS[name], vw = tee || TTE_VIEWS[name] || EXTRA_VIEWS[name];
  const base = tee ? tee.probe() : (TTE_VIEWS[name] ? TTE_VIEWS[name].probe() : EXTRA_VIEWS[name].probe());
  e.viewName = name; e.colorOn = false; e.showLabels = false; e.showColorBox = false;
  if (tee) { e.depthCm = tee.depth; e.freqMHz = 5.5; e.nearFieldCm = 0.3; e.elevFocus = 4.0; e.elevMin = 0.05; e.elevDiv = 0.015; }
  else {
    e.freqMHz = 2.7; e.elevFocus = 8.0; e.elevMin = 0.09; e.elevDiv = 0.025;
    e.nearFieldCm = name === 'SUBCOSTAL' ? 1.6 : /^A\dC$/.test(name) ? 1.4 : 2.0;
    e.depthCm = vw.depth || 16;
  }
  const G = geometryAt(phase, path);
  let probe = { ...base };
  if (vw.track) { const o = vw.track(G.A); probe = { ...base, pos: add(base.pos, o), target: base.target && add(base.target, o) }; }
  e.render(G, path, probe, 0);
  const SW = e.img.width, SH = e.img.height, d = e.img.data;
  const grey = new Float32Array(SW * SH);
  for (let i = 0; i < SW * SH; i++) grey[i] = d[i * 4];
  const apexY = SH * 0.04, pxPerCm = (SH * 0.92) / e.depthCm;
  return { SW, SH, grey, spk: e._spk, kind: e._kind, depth: e._depth, mask: e._mask, abin: e._abin, NA: e._NA,
    lungD: e._lungDepth, cx: SW * 0.5, apexY, pxPerCm, half: e.sectorHalf, G, probe };
}
// sample the grey (and kind) along a beam at angle th (rad, + = screen right)
export function beamSample(R, th, depth) {
  const x = R.cx + Math.sin(th) * depth * R.pxPerCm, y = R.apexY + Math.cos(th) * depth * R.pxPerCm;
  const i = Math.round(x), j = Math.round(y);
  if (i < 0 || j < 0 || i >= R.SW || j >= R.SH) return null;
  const idx = j * R.SW + i;
  return { grey: R.grey[idx], kind: R.kind[idx], idx };
}
// orientation (deg, screen frame, 0 = horizontal, + = clockwise on screen) of
// the major axis of the speckle autocorrelation in a patch centred on beam angle
// th at depth d, and the beam-perpendicular orientation there
export function speckleAxis(R, th, d, half = 12) {
  const cx = Math.round(R.cx + Math.sin(th) * d * R.pxPerCm), cy = Math.round(R.apexY + Math.cos(th) * d * R.pxPerCm);
  const f = R.spk, W = R.SW;
  let m = 0, n = 0;
  for (let y = cy - half; y <= cy + half; y++) for (let x = cx - half; x <= cx + half; x++) { m += f[y * W + x]; n++; }
  m /= n;
  const C = (dx, dy) => {
    let s = 0, k = 0;
    for (let y = cy - half; y <= cy + half; y++) for (let x = cx - half; x <= cx + half; x++) { s += (f[y * W + x] - m) * (f[(y + dy) * W + x + dx] - m); k++; }
    return s / k;
  };
  const c0 = C(0, 0);
  let sxx = 0, syy = 0, sxy = 0;
  for (let dy = -4; dy <= 4; dy++) for (let dx = -6; dx <= 6; dx++) {
    const c = C(dx, dy) / c0; if (c < 0.2) continue;
    sxx += c * dx * dx; syy += c * dy * dy; sxy += c * dx * dy;
  }
  const ang = 0.5 * Math.atan2(2 * sxy, sxx - syy) * 180 / Math.PI;
  const perp = Math.atan2(-Math.sin(th), Math.cos(th)) * 180 / Math.PI;
  return { ang, perp };
}
