// Bakes the 3D panel's thoracic context (thorax3d.js) off the UI thread and
// transfers the typed arrays back; Three.js objects stay on the UI thread.
import { geometryAt } from './cardiac-model.js';
import { buildMesh } from './marching.js';
import { buildContext } from './thorax3d.js';

self.onmessage = ({ data: { heartBox } }) => {
  const meshes = buildContext(geometryAt(0, {}).A, buildMesh, heartBox), transfers = [];
  for (const m of Object.values(meshes)) if (m) transfers.push(m.positions.buffer, m.normals.buffer, m.indices.buffer);
  self.postMessage({ meshes }, transfers);
};
