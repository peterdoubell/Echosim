// CPU-heavy 3D surface baking runs away from the live echo/ECG UI thread.
// Transfer the typed arrays directly; Three.js objects stay on the UI thread.
import { geometryAt } from './cardiac-model.js';
import { myoDist, lumenDist } from './anatomy.js';
import { buildMesh } from './marching.js';

self.onmessage = ({ data: { path, startPhase, phaseCount, bbox, cell } }) => {
  for (let phase = startPhase; phase < phaseCount; phase++) {
    const A = geometryAt(phase / phaseCount, path).A;
    const fields = {
      myo: (x, y, z) => myoDist(x, y, z, A),
      LV: (x, y, z) => lumenDist(x, y, z, A, 'LV'),
      RV: (x, y, z) => lumenDist(x, y, z, A, 'RV'),
      LA: (x, y, z) => lumenDist(x, y, z, A, 'LA'),
      RA: (x, y, z) => lumenDist(x, y, z, A, 'RA'),
      VES: (x, y, z) => Math.min(lumenDist(x, y, z, A, 'AOROOT'), lumenDist(x, y, z, A, 'PA')),
    };
    const meshes = {}, transfers = [];
    for (const [name, field] of Object.entries(fields)) {
      const mesh = meshes[name] = buildMesh(field, bbox, cell);
      if (mesh) transfers.push(mesh.positions.buffer, mesh.normals.buffer, mesh.indices.buffer);
    }
    self.postMessage({ phase, meshes }, transfers);
  }
};
