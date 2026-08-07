// sdf.js — signed-distance-field primitives and smooth CSG operators used to
// assemble an anatomically-structured heart. All functions are pure scalar
// maths (no allocation) so the field can be sampled tens of thousands of times
// per frame for the 2D echo cross-section and on a voxel grid for the 3D
// marching-cubes surface. Distances are in centimetres; negative = inside.

// smooth minimum (rounded union). k = blend radius in cm.
export function smin(a, b, k) {
  const h = Math.max(k - Math.abs(a - b), 0) / k;
  return Math.min(a, b) - h * h * h * k * (1 / 6);
}
// smooth maximum (rounded intersection).
export function smax(a, b, k) {
  const h = Math.max(k - Math.abs(a - b), 0) / k;
  return Math.max(a, b) + h * h * h * k * (1 / 6);
}
// smooth subtraction: carve b out of a.
export function ssub(a, b, k) {
  return smax(a, -b, k);
}

export function sdSphere(px, py, pz, cx, cy, cz, r) {
  const dx = px - cx, dy = py - cy, dz = pz - cz;
  return Math.sqrt(dx * dx + dy * dy + dz * dz) - r;
}

// Approximate ellipsoid SDF (Inigo Quilez bound). Good enough for blending.
export function sdEllipsoid(px, py, pz, cx, cy, cz, rx, ry, rz) {
  const dx = (px - cx), dy = (py - cy), dz = (pz - cz);
  const k1 = Math.sqrt((dx * dx) / (rx * rx) + (dy * dy) / (ry * ry) + (dz * dz) / (rz * rz));
  if (k1 === 0) return -Math.min(rx, ry, rz);
  const k2 = Math.sqrt((dx * dx) / (rx * rx * rx * rx) + (dy * dy) / (ry * ry * ry * ry) + (dz * dz) / (rz * rz * rz * rz));
  return (k1 * (k1 - 1.0)) / k2;
}

// Capsule: distance to the segment a-b, minus radius r (optionally tapered r1->r2).
export function sdCapsule(px, py, pz, ax, ay, az, bx, by, bz, r1, r2) {
  const pax = px - ax, pay = py - ay, paz = pz - az;
  const bax = bx - ax, bay = by - ay, baz = bz - az;
  const baLen2 = bax * bax + bay * bay + baz * baz || 1e-6;
  let h = (pax * bax + pay * bay + paz * baz) / baLen2;
  h = h < 0 ? 0 : h > 1 ? 1 : h;
  const dx = pax - bax * h, dy = pay - bay * h, dz = paz - baz * h;
  const r = r2 === undefined ? r1 : r1 + (r2 - r1) * h;
  return Math.sqrt(dx * dx + dy * dy + dz * dz) - r;
}

// Round cone (tapered capsule with spherical caps of radius r1 at a, r2 at b).
export function sdRoundCone(px, py, pz, ax, ay, az, bx, by, bz, r1, r2) {
  const bax = bx - ax, bay = by - ay, baz = bz - az;
  const l2 = bax * bax + bay * bay + baz * baz || 1e-6;
  const pax = px - ax, pay = py - ay, paz = pz - az;
  const paba = (pax * bax + pay * bay + paz * baz) / l2;
  // approximate: interpolate radius by projection then straight distance
  let h = paba; h = h < 0 ? 0 : h > 1 ? 1 : h;
  const cx = ax + bax * h, cy = ay + bay * h, cz = az + baz * h;
  const r = r1 + (r2 - r1) * h;
  const dx = px - cx, dy = py - cy, dz = pz - cz;
  return Math.sqrt(dx * dx + dy * dy + dz * dz) - r;
}
