// marching.js — build a triangle mesh from a signed-distance field using Naive
// Surface Nets (Lysenko). Chosen over classic marching cubes because it needs
// no 256-row triangle table and yields smoother, watertight organic surfaces —
// ideal for the blended-SDF heart. Returns INDEXED position + normal arrays plus
// an index buffer, ready for a THREE.BufferGeometry (setIndex).

const cubeEdges = new Int32Array(24);
const edgeTable = new Int32Array(256);
(function init() {
  let k = 0;
  for (let i = 0; i < 8; ++i) {
    for (let j = 1; j <= 4; j <<= 1) {
      const p = i ^ j;
      if (i <= p) { cubeEdges[k++] = i; cubeEdges[k++] = p; }
    }
  }
  for (let i = 0; i < 256; ++i) {
    let em = 0;
    for (let j = 0; j < 24; j += 2) {
      const a = !!(i & (1 << cubeEdges[j]));
      const b = !!(i & (1 << cubeEdges[j + 1]));
      em |= a !== b ? (1 << (j >> 1)) : 0;
    }
    edgeTable[i] = em;
  }
})();

// Sample the field on a grid and extract the zero-isosurface.
// fieldFn(x,y,z) -> signed distance (negative inside).
// bbox = {min:[x,y,z], max:[x,y,z]}, cell = voxel size (cm).
// Returns { positions: Float32Array, normals: Float32Array, indices: Uint32Array |
// Uint16Array } (indexed triangles), or null if the surface is empty.
export function buildMesh(fieldFn, bbox, cell) {
  const dims = [
    Math.max(2, Math.ceil((bbox.max[0] - bbox.min[0]) / cell) + 1),
    Math.max(2, Math.ceil((bbox.max[1] - bbox.min[1]) / cell) + 1),
    Math.max(2, Math.ceil((bbox.max[2] - bbox.min[2]) / cell) + 1),
  ];
  const [nx, ny, nz] = dims;
  const data = new Float32Array(nx * ny * nz);
  let p = 0;
  for (let z = 0; z < nz; ++z) {
    const wz = bbox.min[2] + z * cell;
    for (let y = 0; y < ny; ++y) {
      const wy = bbox.min[1] + y * cell;
      for (let x = 0; x < nx; ++x, ++p) {
        data[p] = fieldFn(bbox.min[0] + x * cell, wy, wz);
      }
    }
  }

  const verts = [];        // grid-space vertices
  const quads = [];        // index quads
  const R = [1, nx + 1, (nx + 1) * (ny + 1)];
  const buffer = new Int32Array(R[2] * 2);
  const grid = new Float32Array(8);
  let bufNo = 1;
  let n = 0;
  const xpos = [0, 0, 0];

  for (xpos[2] = 0; xpos[2] < nz - 1; ++xpos[2], n += nx, bufNo ^= 1, R[2] = -R[2]) {
    let m = 1 + (nx + 1) * (1 + bufNo * (ny + 1));
    for (xpos[1] = 0; xpos[1] < ny - 1; ++xpos[1], ++n, m += 2) {
      for (xpos[0] = 0; xpos[0] < nx - 1; ++xpos[0], ++n, ++m) {
        let mask = 0, g = 0, idx = n;
        for (let k = 0; k < 2; ++k, idx += nx * (ny - 2)) {
          for (let j = 0; j < 2; ++j, idx += nx - 2) {
            for (let i = 0; i < 2; ++i, ++g, ++idx) {
              const val = data[idx];
              grid[g] = val;
              mask |= (val < 0) ? (1 << g) : 0;
            }
          }
        }
        if (mask === 0 || mask === 255) continue;
        const em = edgeTable[mask];
        const v = [0, 0, 0];
        let eCount = 0;
        for (let i = 0; i < 12; ++i) {
          if (!(em & (1 << i))) continue;
          ++eCount;
          const e0 = cubeEdges[i << 1], e1 = cubeEdges[(i << 1) + 1];
          const g0 = grid[e0], g1 = grid[e1];
          let t = g0 - g1;
          if (Math.abs(t) > 1e-6) t = g0 / t; else t = 0.5;
          for (let j = 0, kk = 1; j < 3; ++j, kk <<= 1) {
            const a = e0 & kk, b = e1 & kk;
            if (a !== b) v[j] += a ? 1.0 - t : t;
            else v[j] += a ? 1.0 : 0;
          }
        }
        const s = 1.0 / eCount;
        v[0] = xpos[0] + s * v[0];
        v[1] = xpos[1] + s * v[1];
        v[2] = xpos[2] + s * v[2];
        buffer[m] = verts.length / 3;
        verts.push(v[0], v[1], v[2]);
        for (let i = 0; i < 3; ++i) {
          if (!(em & (1 << i))) continue;
          const iu = (i + 1) % 3, iv = (i + 2) % 3;
          if (xpos[iu] === 0 || xpos[iv] === 0) continue;
          const du = R[iu], dv = R[iv];
          if (mask & 1) {
            quads.push(buffer[m], buffer[m - du], buffer[m - du - dv], buffer[m - dv]);
          } else {
            quads.push(buffer[m], buffer[m - dv], buffer[m - du - dv], buffer[m - du]);
          }
        }
      }
    }
  }

  if (verts.length === 0 || quads.length === 0) return null;

  // INDEXED output: each surface-net vertex is unique, so map it to world space
  // and evaluate its gradient normal exactly ONCE (not once per triangle-corner).
  // This roughly thirds both the vertex count and the normal-field evaluations.
  const nVerts = verts.length / 3;
  const positions = new Float32Array(nVerts * 3);
  const normals = new Float32Array(nVerts * 3);
  const eps = cell * 0.5;
  for (let i = 0; i < nVerts; ++i) {
    const wx = bbox.min[0] + verts[i * 3] * cell;
    const wy = bbox.min[1] + verts[i * 3 + 1] * cell;
    const wz = bbox.min[2] + verts[i * 3 + 2] * cell;
    positions[i * 3] = wx; positions[i * 3 + 1] = wy; positions[i * 3 + 2] = wz;
    let nxv = fieldFn(wx + eps, wy, wz) - fieldFn(wx - eps, wy, wz);
    let nyv = fieldFn(wx, wy + eps, wz) - fieldFn(wx, wy - eps, wz);
    let nzv = fieldFn(wx, wy, wz + eps) - fieldFn(wx, wy, wz - eps);
    const l = Math.hypot(nxv, nyv, nzv) || 1;
    normals[i * 3] = nxv / l; normals[i * 3 + 1] = nyv / l; normals[i * 3 + 2] = nzv / l;
  }

  // expand quads -> two triangles each, preserving the original winding
  const triQuads = quads.length / 4;
  const IndexArray = nVerts > 65535 ? Uint32Array : Uint16Array;
  const indices = new IndexArray(triQuads * 6);
  let o = 0;
  for (let q = 0; q < quads.length; q += 4) {
    const a = quads[q], b = quads[q + 1], c = quads[q + 2], d = quads[q + 3];
    indices[o++] = a; indices[o++] = b; indices[o++] = c;
    indices[o++] = a; indices[o++] = c; indices[o++] = d;
  }
  return { positions, normals, indices };
}
