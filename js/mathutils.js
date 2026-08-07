// mathutils.js
// A tiny, dependency-free collection of scalar easing helpers and plain-array
// vector maths shared across the EchoSim modules (cardiac-model.js, echo.js,
// main.js). Everything is pure: no allocation of hidden state, no imports.
//
// Vectors are plain JS arrays [x, y, z]. The vec helpers return fresh arrays and
// never mutate their inputs, so callers can treat the maths as value-semantics.

// clamp v into [a, b]
export const clamp = (v, a, b) => (v < a ? a : v > b ? b : v);

// classic Hermite smoothstep: 0 below a, 1 above b, smooth in between
export const smoothstep = (a, b, x) => {
  const t = clamp((x - a) / (b - a), 0, 1);
  return t * t * (3 - 2 * t);
};

// a smooth 0->1->0 bump over [a,b]
export const pulse = (a, b, x) => {
  if (x <= a || x >= b) return 0;
  const t = (x - a) / (b - a);
  return Math.sin(Math.PI * t);
};

// Euclidean length of a vector.
export const vlen = (v) => Math.hypot(...v);

// unit vector (safe on the zero vector: returns it unchanged rather than NaN).
export function unit(v) {
  const l = vlen(v) || 1;
  return [v[0] / l, v[1] / l, v[2] / l];
}
// `norm` is a conventional alias for `unit`.
export const norm = unit;

// --- plain-array vector helpers (3-component) ---
export const vadd = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
export const vsub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
export const vscale = (v, s) => [v[0] * s, v[1] * s, v[2] * s];
export const vdot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
export const vcross = (a, b) => [
  a[1] * b[2] - a[2] * b[1],
  a[2] * b[0] - a[0] * b[2],
  a[0] * b[1] - a[1] * b[0],
];

// Rotate vector v about `axis` (need not be unit) by `ang` radians using
// Rodrigues' rotation formula: v*cos + (k x v)*sin + k*(k.v)*(1-cos).
export function vrot(v, axis, ang) {
  const k = unit(axis);
  const c = Math.cos(ang), s = Math.sin(ang);
  const kv = vdot(k, v);
  const cr = vcross(k, v);
  return [
    v[0] * c + cr[0] * s + k[0] * kv * (1 - c),
    v[1] * c + cr[1] * s + k[1] * kv * (1 - c),
    v[2] * c + cr[2] * s + k[2] * kv * (1 - c),
  ];
}
