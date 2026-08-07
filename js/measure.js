// measure.js — learner-placed measurement tools drawn on a transparent overlay
// canvas above the B-mode sector (B8: quantification & measurement accuracy).
//
// Tools:
//   • caliper   — click two points, read the true distance in cm.
//   • Simpson    — trace the LV endocardium at end-diastole and end-systole; the
//                 module computes EDV/ESV/EF by the method of discs (Simpson
//                 single-plane) and reports it against the model's own ground
//                 truth (the same §A2 measurement-accuracy check as the harness).
//
// Coordinate model: both #echo and #echoOverlay are 440×500 rasters shown with
// `object-fit: contain`, so a client point maps to the shared 440-space canvas
// coordinate after removing the letterbox. A displacement of `d` canvas pixels is
// `d × echo.cmPerPx` centimetres (the sector scales radial and lateral identically),
// so linear distances and disc areas are exact.

const ND = 20; // Simpson discs

export class MeasureTool {
  constructor(overlay, echo, opts = {}) {
    this.cv = overlay;
    this.ctx = overlay.getContext('2d');
    this.echo = echo;
    this.tool = 'off';          // 'off' | 'caliper' | 'simpsonED' | 'simpsonES'
    this.caliper = null;        // { a:{x,y}, b:{x,y} } in canvas space
    this._pendingA = null;      // first caliper click awaiting the second
    this.traceED = null;        // array of {x,y}
    this.traceES = null;
    this._trace = null;         // trace currently being drawn
    this.hover = null;          // live cursor point
    this.result = { distCm: null, edv: null, esv: null, ef: null, which: null };
    this.onResult = opts.onResult || (() => {});
    this.groundTruth = opts.groundTruth || (() => ({}));  // () => { ef, lvidd }

    this._onDown = this._down.bind(this);
    this._onMove = this._move.bind(this);
    this._onDbl = this._dbl.bind(this);
    this._onLeave = () => { this.hover = null; };
    overlay.addEventListener('pointerdown', this._onDown);
    overlay.addEventListener('pointermove', this._onMove);
    overlay.addEventListener('dblclick', this._onDbl);
    overlay.addEventListener('pointerleave', this._onLeave);
  }

  setTool(t) {
    this.tool = (this.tool === t) ? 'off' : t;
    this._pendingA = null;
    if (this.tool === 'simpsonED') this._trace = this.traceED = [];
    else if (this.tool === 'simpsonES') this._trace = this.traceES = [];
    else this._trace = null;
    this.cv.classList.toggle('measuring', this.tool !== 'off');
    return this.tool;
  }

  clear() {
    this.caliper = this._pendingA = this._trace = null;
    this.traceED = this.traceES = null;
    this.result = { distCm: null, edv: null, esv: null, ef: null, which: null };
    this.onResult(this.result);
  }

  // client (mouse) → shared 440×500 canvas coordinate, undoing object-fit:contain
  _toCanvas(e) {
    const r = this.cv.getBoundingClientRect();
    const cw = this.cv.width, ch = this.cv.height;
    const s = Math.min(r.width / cw, r.height / ch);
    const padX = (r.width - cw * s) / 2, padY = (r.height - ch * s) / 2;
    const x = (e.clientX - r.left - padX) / s;
    const y = (e.clientY - r.top - padY) / s;
    return { x, y, inside: x >= 0 && x <= cw && y >= 0 && y <= ch };
  }

  _down(e) {
    if (this.tool === 'off') return;
    const p = this._toCanvas(e);
    if (!p.inside) return;
    e.preventDefault();
    if (this.tool === 'caliper') {
      if (!this._pendingA) { this._pendingA = { x: p.x, y: p.y }; }
      else {
        this.caliper = { a: this._pendingA, b: { x: p.x, y: p.y } };
        this._pendingA = null;
        this._emitCaliper();
      }
    } else if (this._trace) {
      this._trace.push({ x: p.x, y: p.y });
    }
  }

  _move(e) {
    if (this.tool === 'off') { this.hover = null; return; }
    const p = this._toCanvas(e);
    this.hover = p.inside ? { x: p.x, y: p.y } : null;
  }

  _dbl(e) {
    if (this.tool !== 'simpsonED' && this.tool !== 'simpsonES') return;
    e.preventDefault();
    // finalize the current trace (needs ≥6 points to be a usable contour)
    if (this._trace && this._trace.length >= 6) {
      this._compute();
      this.setTool('off');
    }
  }

  _emitCaliper() {
    const c = this.caliper;
    const d = Math.hypot(c.a.x - c.b.x, c.a.y - c.b.y) * this.echo.cmPerPx;
    this.result.distCm = d;
    this.onResult(this.result);
  }

  // Simpson method of discs on a traced endocardial contour. Clinical convention:
  // the trace runs from one mitral hinge, around the apex, to the other hinge, so the
  // two endpoints define the base; the long axis is apex → base-midpoint and the
  // cavity is closed by the base chord. Returns cavity volume in mL, or null.
  _discVolume(trace) {
    if (!trace || trace.length < 6) return null;
    const cm = this.echo.cmPerPx;
    const P = trace.map((q) => ({ x: q.x * cm, y: q.y * cm })); // → cm space
    const E0 = P[0], E1 = P[P.length - 1];
    const bm = { x: (E0.x + E1.x) / 2, y: (E0.y + E1.y) / 2 }; // base midpoint (mitral plane)
    // apex = contour point farthest from the base midpoint
    let apex = P[0], L = -1;
    for (const q of P) { const d = Math.hypot(q.x - bm.x, q.y - bm.y); if (d > L) { L = d; apex = q; } }
    if (L < 1) return null;
    const ux = (apex.x - bm.x) / L, uy = (apex.y - bm.y) / L;  // long-axis unit (base→apex)
    const px = -uy, py = ux;                                    // perpendicular unit
    // project every contour point onto (long from base midpoint, perp)
    const proj = P.map((q) => ({ t: (q.x - bm.x) * ux + (q.y - bm.y) * uy, s: (q.x - bm.x) * px + (q.y - bm.y) * py }));
    const h = L / ND;
    let vol = 0;
    for (let i = 0; i < ND; i++) {
      const tc = (i + 0.5) * h;
      let smin = Infinity, smax = -Infinity;
      // perpendicular extent of the contour where it crosses this slab centre
      for (let k = 0; k < proj.length; k++) {
        const a = proj[k], b = proj[(k + 1) % proj.length]; // treat as closed
        if ((a.t - tc) * (b.t - tc) <= 0 && a.t !== b.t) {
          const f = (tc - a.t) / (b.t - a.t);
          const s = a.s + f * (b.s - a.s);
          if (s < smin) smin = s; if (s > smax) smax = s;
        }
      }
      if (smax > smin) { const d = smax - smin; vol += Math.PI * (d * d) / 4 * h; }
    }
    return vol > 0 ? vol : null; // cm³ = mL
  }

  _compute() {
    if (this.traceED) { const v = this._discVolume(this.traceED); if (v != null) this.result.edv = v; }
    if (this.traceES) { const v = this._discVolume(this.traceES); if (v != null) this.result.esv = v; }
    const { edv, esv } = this.result;
    if (edv != null && esv != null && edv > 0) {
      this.result.ef = (edv - esv) / edv * 100;
      this.result.which = 'EF';
    } else {
      this.result.which = edv != null ? 'ED only' : (esv != null ? 'ES only' : null);
    }
    this.onResult(this.result);
  }

  // redraw the overlay each animation frame
  draw() {
    const ctx = this.ctx, W = this.cv.width, H = this.cv.height;
    ctx.clearRect(0, 0, W, H);
    // finalized caliper
    if (this.caliper) this._line(this.caliper.a, this.caliper.b, '#ffd24a', true);
    // pending caliper rubber-band
    if (this.tool === 'caliper' && this._pendingA) {
      this._dot(this._pendingA, '#ffd24a');
      if (this.hover) this._line(this._pendingA, this.hover, 'rgba(255,210,74,0.6)', false);
    }
    // traces
    this._drawTrace(this.traceED, '#57d9a3');
    this._drawTrace(this.traceES, '#5aa9e8');
    if (this._trace && this.hover && (this.tool === 'simpsonED' || this.tool === 'simpsonES')) {
      const last = this._trace[this._trace.length - 1];
      if (last) this._line(last, this.hover, 'rgba(255,255,255,0.4)', false);
    }
  }

  _drawTrace(tr, col) {
    if (!tr || tr.length === 0) return;
    const ctx = this.ctx;
    ctx.strokeStyle = col; ctx.lineWidth = 1.6; ctx.beginPath();
    tr.forEach((q, i) => (i ? ctx.lineTo(q.x, q.y) : ctx.moveTo(q.x, q.y)));
    ctx.stroke();
    ctx.fillStyle = col;
    for (const q of tr) { ctx.beginPath(); ctx.arc(q.x, q.y, 2, 0, Math.PI * 2); ctx.fill(); }
  }

  _line(a, b, col, ticks) {
    const ctx = this.ctx;
    ctx.strokeStyle = col; ctx.lineWidth = 1.6;
    ctx.beginPath(); ctx.moveTo(a.x, a.y); ctx.lineTo(b.x, b.y); ctx.stroke();
    this._dot(a, col); this._dot(b, col);
    if (ticks) {
      const d = Math.hypot(a.x - b.x, a.y - b.y) * this.echo.cmPerPx;
      ctx.fillStyle = col; ctx.font = 'bold 12px monospace'; ctx.textAlign = 'left';
      ctx.fillText(d.toFixed(1) + ' cm', (a.x + b.x) / 2 + 6, (a.y + b.y) / 2 - 4);
    }
  }
  _dot(p, col) { const ctx = this.ctx; ctx.fillStyle = col; ctx.beginPath(); ctx.arc(p.x, p.y, 3, 0, Math.PI * 2); ctx.fill(); }
}
