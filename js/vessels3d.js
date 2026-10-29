// vessels3d.js — the great vessels, coronaries and coronary sinus of the 3D
// panel, as distance fields for the marching-cubes mesher (marching.js). Kept
// free of Three.js so the node audits can mesh and count them.
//
// Colour follows the panel's legend: oxygenated (aorta, its arch and head and
// neck branches, the descending aorta) red, the coronary arteries a lighter red;
// deoxygenated (pulmonary trunk and branches, coronary sinus, both venae cavae
// and the left brachiocephalic vein) venous blue.

import { lumenDist } from './anatomy.js';

export const ARTERIAL = 0xc0392b;      // = the LV / 'Oxygenated / left' legend swatch
export const CORONARY = 0xff7f66;      // lighter arterial red (small epicardial arteries)
export const VENOUS = 0x2b6cb0;        // = the RV / 'Deoxygenated / right' legend swatch
export const VENOUS_DEEP = 0x203f79;   // cavae and coronary sinus: the RA's deeper venous tint

// The heart box (heart3d BBOX) holds the moving structures; the arch, its
// branches, the descending aorta and the cavae run beyond it, so the static
// vessels are meshed once in a taller box (body superior is mostly heart +z).
export const VESSEL_BOX = { min: [-8.4, -5.0, -10.5], max: [5.6, 9.6, 11.4] };
export const VESSEL_CELL = 0.4;

// moving: re-meshed per cardiac-phase keyframe with the chambers (heart box);
// static: meshed once from the end-diastolic anatomy (VESSEL_BOX). The small
// coronaries and the CS are drawn a little fatter than their lumens so the
// heart-box voxel grid resolves them.
export const VESSELS = [
  { key: 'aorticRoot', color: ARTERIAL, moving: true, fn: (x, y, z, A) => lumenDist(x, y, z, A, 'AOROOT') },
  { key: 'PA', color: VENOUS, moving: true, fn: (x, y, z, A) => lumenDist(x, y, z, A, 'PA') },
  { key: 'coronary', color: CORONARY, moving: true, fn: (x, y, z, A) => lumenDist(x, y, z, A, 'CORONARY') - 0.16 },
  { key: 'CS', color: VENOUS_DEEP, moving: true, fn: (x, y, z, A) => lumenDist(x, y, z, A, 'CS') - 0.08 },
  { key: 'arch', color: ARTERIAL, moving: false, fn: (x, y, z, A) => lumenDist(x, y, z, A, 'ARCH') },
  { key: 'desc', color: ARTERIAL, moving: false, fn: (x, y, z, A) => lumenDist(x, y, z, A, 'DESC') },
  { key: 'SVC', color: VENOUS_DEEP, moving: false, fn: (x, y, z, A) => Math.min(lumenDist(x, y, z, A, 'SVC'), lumenDist(x, y, z, A, 'BCV')) },
  { key: 'IVC', color: VENOUS_DEEP, moving: false, fn: (x, y, z, A) => lumenDist(x, y, z, A, 'IVC') },
];
