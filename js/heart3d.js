// heart3d.js
// The interactive 3D "anatomy" view: a continuously beating schematic heart,
// a moveable ultrasound transducer with its imaging fan (the plane the echo
// view samples), and Doppler-coloured flow particles. Built on Three.js.

import * as THREE from 'three';
import { OrbitControls } from '../vendor/OrbitControls.js';
import { geometryAt, velocityAt, VALVE_DEFS } from './cardiac-model.js';
import { myoDist, lumenDist, LM } from './anatomy.js';
import { buildMesh } from './marching.js';

const MUSCLE = 0x8a3330;
// Chamber tints. Ventricles read bright; atria are a touch deeper/darker and
// slightly hue-shifted so the four overlapping translucent pools separate into
// distinct volumes instead of one muddy blob. BART convention preserved:
// red = left heart, blue = right heart.
const LV_COLOR = 0xc0392b; // bright arterial red
const LA_COLOR = 0x8e2036; // deeper crimson (darker + shifted toward burgundy)
const RV_COLOR = 0x2b6cb0; // bright venous blue
const RA_COLOR = 0x203f79; // deeper indigo (darker + shifted)

// Voxel volume + resolution for the marching-cubes heart surface, and how many
// cardiac-phase keyframes to pre-bake (swapped per frame for the beating loop).
// The box frames the heart at adult size (~12.5 x 15 x 11 cm with the root and
// PA) and clips the cavae; the cell is set so a phase bakes in about the same
// time as the old compact schematic did.
const BBOX = { min: [-7.4, -9.4, -4.6], max: [5.2, 6.0, 6.6] };
const CELL = 0.36;
const NPHASE = 12;

export class Heart3D {
  constructor(container) {
    this.container = container;
    const w = container.clientWidth, h = container.clientHeight;

    this.renderer = new THREE.WebGLRenderer({ antialias: true, alpha: false });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    this.renderer.setSize(w, h);
    this.renderer.localClippingEnabled = true;
    // Cinematic tone + correct colour management: filmic roll-off keeps the
    // bright blood pools from clipping and gives the muscle real gradation.
    this.renderer.toneMapping = THREE.ACESFilmicToneMapping;
    this.renderer.toneMappingExposure = 1.1;
    this.renderer.outputColorSpace = THREE.SRGBColorSpace;
    container.appendChild(this.renderer.domElement);

    this.scene = new THREE.Scene();
    this.scene.background = new THREE.Color(0x0a0e16);

    this.camera = new THREE.PerspectiveCamera(42, w / h, 0.1, 200);
    this.camera.position.set(11, 4, 13);

    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = true;
    this.controls.dampingFactor = 0.08;
    this.controls.target.set(0, 0.2, 0);

    // lights: low fill so form is carved by a strong warm key and a cool
    // back/rim light that separates the heart from the dark background.
    this.scene.add(new THREE.AmbientLight(0xffffff, 0.25));
    const key = new THREE.DirectionalLight(0xfff1e0, 2.1);
    key.position.set(7, 11, 9);
    this.scene.add(key);
    const fill = new THREE.DirectionalLight(0xffffff, 0.35);
    fill.position.set(4, -2, 8);
    this.scene.add(fill);
    const rim = new THREE.DirectionalLight(0x7fa8ff, 1.3);
    rim.position.set(-9, 2, -8);
    this.scene.add(rim);

    // subtle ground grid for spatial reference
    const grid = new THREE.GridHelper(40, 40, 0x1f3550, 0x14202f);
    grid.position.y = BBOX.min[1] - 1.2;                // just below the apex
    this.scene.add(grid);

    this.heartGroup = new THREE.Group();
    this.scene.add(this.heartGroup);

    this._buildHeart();
    this._buildProbe();
    this._buildFlow();

    window.addEventListener('resize', () => this._resize());
    this.showFlow = true;
    this.transparent = true;
  }

  _buildHeart() {
    const mkMat = (color, opacity, extra = {}) => new THREE.MeshStandardMaterial({
      color, transparent: true, opacity, roughness: 0.55, metalness: 0.0,
      side: THREE.FrontSide, ...extra,
    });
    // myocardium: a single-pass FrontSide shell (no stacked double-sided
    // translucency) with a subtle Fresnel rim so the form and chamber walls
    // read cleanly while blood pools stay legible through it.
    this.myoMesh = new THREE.Mesh(new THREE.BufferGeometry(), mkMat(MUSCLE, 0.4, {
      depthWrite: false, roughness: 0.62,
      emissive: new THREE.Color(0x2a0f0d), emissiveIntensity: 0.25,
    }));
    this._addFresnelRim(this.myoMesh.material, 0xff8a72, 2.4, 0.9);
    this.myoMesh.renderOrder = 4;
    // blood pools + great vessels (retuned for the sRGB / filmic pipeline).
    // Ventricles are brighter/less transparent; atria deeper and slightly denser
    // so paired chambers separate by value as well as hue.
    this.lvMesh = new THREE.Mesh(new THREE.BufferGeometry(), mkMat(LV_COLOR, 0.72, { depthWrite: true, roughness: 0.42 }));
    this.rvMesh = new THREE.Mesh(new THREE.BufferGeometry(), mkMat(RV_COLOR, 0.72, { depthWrite: true, roughness: 0.42 }));
    this.laMesh = new THREE.Mesh(new THREE.BufferGeometry(), mkMat(LA_COLOR, 0.62, { depthWrite: true, roughness: 0.42 }));
    this.raMesh = new THREE.Mesh(new THREE.BufferGeometry(), mkMat(RA_COLOR, 0.62, { depthWrite: true, roughness: 0.42 }));
    this.vesMesh = new THREE.Mesh(new THREE.BufferGeometry(), mkMat(0xb85c40, 0.66, { depthWrite: true, roughness: 0.5 }));
    // Per-chamber Fresnel rim: each pool gets a thin glowing silhouette in its
    // own hue, so the depth-sorted translucent volumes read as separate bodies
    // (a shader-side outline that needs no extra per-phase geometry).
    this._addFresnelRim(this.lvMesh.material, 0xff7a63, 2.2, 0.55);
    this._addFresnelRim(this.laMesh.material, 0xcf5a72, 2.2, 0.55);
    this._addFresnelRim(this.rvMesh.material, 0x6fb2ff, 2.2, 0.55);
    this._addFresnelRim(this.raMesh.material, 0x5f83d6, 2.2, 0.55);
    this.bloodMeshes = [this.lvMesh, this.rvMesh, this.laMesh, this.raMesh, this.vesMesh];
    this.bloodMeshes.forEach((m) => { m.renderOrder = 2; });
    [this.myoMesh, ...this.bloodMeshes].forEach((m) => this.heartGroup.add(m));

    // Rule-based (Streeter) myofiber overlay: helical fibres wound around the LV
    // long axis whose helix angle sweeps transmurally from ~+60° (endocardium) to
    // ~−60° (epicardium). Off by default; built lazily from the LV geometry.
    this.fiberLines = new THREE.LineSegments(
      new THREE.BufferGeometry(),
      new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, opacity: 0.9, depthWrite: false }),
    );
    this.fiberLines.visible = false;
    this.fiberLines.renderOrder = 6;
    this.heartGroup.add(this.fiberLines);
    this._fiberKey = null;

    // Subvalvular apparatus: papillary muscles (LV ×2 + RV anterior), the
    // moderator band, and the chordae-tendineae fans up to the mitral & tricuspid
    // leaflet rings — sourced from the shared anatomy so the 3D matches the echo.
    this.subValv = null;
    this._buildSubvalvular(geometryAt(0, {}).A);

    // valves as slim annular rings (open/close cue). Geometry from shared VALVE_DEFS.
    this.valveMeshes = {};
    const valveColor = { mitral: 0xe8d8b0, tricuspid: 0xe8d8b0, aortic: 0xf0e4c0, pulmonic: 0xf0e4c0 };
    for (const name in VALVE_DEFS) {
      const d = VALVE_DEFS[name];
      // thin annulus: inner radius near the outer so the ring reads as a rim
      const geo = new THREE.RingGeometry(d.r * 0.82, d.r, 40, 1);
      const mat = new THREE.MeshStandardMaterial({ color: valveColor[name], side: THREE.DoubleSide, transparent: true, opacity: 0.6, roughness: 0.5, emissive: 0x332a10, emissiveIntensity: 0.25, depthWrite: false });
      const m = new THREE.Mesh(geo, mat);
      const nrm = new THREE.Vector3(d.n[0], d.n[1], d.n[2]).normalize();
      m.position.set(d.c[0], d.c[1], d.c[2]);
      m.quaternion.setFromUnitVectors(new THREE.Vector3(0, 0, 1), nrm);
      m.userData.normal = nrm; // world-space annulus normal for edge-on fade
      m.renderOrder = 5;
      this.heartGroup.add(m);
      this.valveMeshes[name] = m;
    }

    // Async-bake cache (see setPathology). Only phase 0 of the normal heart is
    // baked synchronously so the first frame has geometry to bind; the other
    // phases are deferred over the same time-sliced path (see _bakeInitial).
    this._cache = new Map();     // pathKey -> frames[] (LRU by insertion order)
    this._cacheMax = 4;
    this._bakeToken = 0;         // bumped on every setPathology to supersede bakes
    this._overlayShown = false;  // guards onBakeStart/onBakeEnd balance
    // Frame-budget-adaptive time slicing: start at 2 phases/tick and let the
    // bake loop steer this between 1 and _maxPhasesPerTick from measured cost.
    this._phasesPerTick = 2;
    this._maxPhasesPerTick = 4;
    this._tickEMA = null;        // EMA of per-phase bake cost (ms)
    this._bakeInitial();
  }

  // Initial normal-heart bake. Phase 0 is baked synchronously so the very first
  // frame can bind real geometry; the remaining phases are deferred to the next
  // animation frame and time-sliced, with onBakeStart/onBakeEnd bracketing that
  // deferred portion. Deferring to rAF also ensures onBakeStart fires after the
  // caller has wired up its callbacks.
  _bakeInitial() {
    const path = {};
    const key = this._pathKey(path);
    const token = ++this._bakeToken;
    const frames = [this._bakePhase(0, path)];
    this._bindFrames(frames); // first paint binds phase 0 immediately
    requestAnimationFrame(() => {
      if (token !== this._bakeToken) return; // superseded before we began
      this._showOverlay();
      this._runBake(path, token, frames, 1).then((done) => {
        if (!done) return; // a newer setPathology owns the overlay + binding
        this._cacheStore(key, frames);
        this._bindFrames(frames);
        this._hideOverlay();
      });
    });
  }

  // Add a view-dependent Fresnel rim to a standard material via a tiny
  // onBeforeCompile patch (grazing angles glow, so the silhouette reads).
  _addFresnelRim(mat, color, power, strength) {
    mat.onBeforeCompile = (shader) => {
      shader.uniforms.uRimColor = { value: new THREE.Color(color) };
      shader.uniforms.uRimPower = { value: power };
      shader.uniforms.uRimStrength = { value: strength };
      shader.fragmentShader = 'uniform vec3 uRimColor;\nuniform float uRimPower;\nuniform float uRimStrength;\n' + shader.fragmentShader;
      shader.fragmentShader = shader.fragmentShader.replace(
        '#include <emissivemap_fragment>',
        '#include <emissivemap_fragment>\n' +
        '  float rimF = pow(1.0 - abs(dot(normalize(normal), normalize(vViewPosition))), uRimPower);\n' +
        '  totalEmissiveRadiance += uRimColor * rimF * uRimStrength;'
      );
    };
    mat.needsUpdate = true;
  }

  // Bake the six structure geometries for a single cardiac-phase keyframe.
  _bakePhase(p, path) {
    const A = geometryAt(p / NPHASE, path).A;
    return {
      myo: this._geomFrom((x, y, z) => myoDist(x, y, z, A)),
      LV: this._geomFrom((x, y, z) => lumenDist(x, y, z, A, 'LV')),
      RV: this._geomFrom((x, y, z) => lumenDist(x, y, z, A, 'RV')),
      LA: this._geomFrom((x, y, z) => lumenDist(x, y, z, A, 'LA')),
      RA: this._geomFrom((x, y, z) => lumenDist(x, y, z, A, 'RA')),
      VES: this._geomFrom((x, y, z) => Math.min(lumenDist(x, y, z, A, 'AOROOT'), lumenDist(x, y, z, A, 'PA'))),
    };
  }

  // Steer _phasesPerTick from the measured cost of the phases just baked.
  // A short EMA smooths the per-phase estimate; an 0.8x/1.2x dead-band around
  // the frame budget gives hysteresis so the rate doesn't oscillate. Weak
  // hardware settles toward 1 phase/tick so a re-bake never drops a frame.
  _adaptSlice(ms, phasesDone) {
    if (phasesDone <= 0) return;
    const perPhase = ms / phasesDone;
    const a = 0.3;
    this._tickEMA = this._tickEMA == null ? perPhase : this._tickEMA * (1 - a) + perPhase * a;
    const BUDGET = 8; // ms of bake work we allow per animation frame
    const cur = this._phasesPerTick;
    if (this._tickEMA * (cur + 1) < BUDGET * 0.8 && cur < this._maxPhasesPerTick) {
      this._phasesPerTick = cur + 1;
    } else if (this._tickEMA * cur > BUDGET * 1.2 && cur > 1) {
      this._phasesPerTick = cur - 1;
    }
  }

  // Time-sliced bake driver shared by the initial bake and setPathology. Bakes
  // phases [startP, NPHASE) into `frames` across rAF ticks, adapting the slice
  // size each tick. Resolves true when the full set is baked, or false if a
  // newer bake (token bump) superseded this one mid-flight.
  _runBake(path, token, frames, startP) {
    return new Promise((resolve) => {
      let p = startP;
      const step = () => {
        if (token !== this._bakeToken) { resolve(false); return; }
        const t0 = performance.now();
        const end = Math.min(NPHASE, p + this._phasesPerTick);
        let done = 0;
        for (; p < end; p++) { frames.push(this._bakePhase(p, path)); done++; }
        this._adaptSlice(performance.now() - t0, done);
        if (p < NPHASE) requestAnimationFrame(step);
        else resolve(true);
      };
      requestAnimationFrame(step);
    });
  }

  _bindFrames(frames) {
    this.frames = frames;
    this._boundFrame = null; // force a re-bind of the per-phase geometry next update
  }

  static _disposeFrames(frames) {
    for (const f of frames) for (const k in f) if (f[k]) f[k].dispose();
  }

  // Canonical cache key: JSON of the sorted truthy [flag, value] pairs, so two
  // pathologies with the same active flags share a bake and never collide.
  _pathKey(path) {
    const entries = Object.keys(path || {})
      .filter((k) => path[k])
      .sort()
      .map((k) => [k, path[k]]);
    return JSON.stringify(entries);
  }

  _showOverlay() {
    if (!this._overlayShown) { this._overlayShown = true; this.onBakeStart && this.onBakeStart(); }
  }

  _hideOverlay() {
    if (this._overlayShown) { this._overlayShown = false; this.onBakeEnd && this.onBakeEnd(); }
  }

  // Insert a freshly baked frame set into the LRU cache, evicting the oldest
  // entries (never the currently bound set) if we exceed the cap.
  _cacheStore(key, frames) {
    this._cache.set(key, frames);
    while (this._cache.size > this._cacheMax) {
      const oldestKey = this._cache.keys().next().value;
      const old = this._cache.get(oldestKey);
      this._cache.delete(oldestKey);
      if (old !== this.frames) Heart3D._disposeFrames(old);
    }
  }

  // Rebind the beating loop for a new pathology. NON-BLOCKING and cached:
  //   - cache HIT: bind the stored frames immediately.
  //   - cache MISS: bake the 12 phases time-sliced across rAF ticks so the main
  //     thread never freezes, then bind and cache the result.
  // Fires this.onBakeStart()/this.onBakeEnd() around a real bake and returns a
  // Promise that resolves once the frames are bound (or the bake is superseded).
  setPathology(path) {
    path = path || {};
    this._buildSubvalvular(geometryAt(0, path).A); // papillaries/chordae track dilatation
    const key = this._pathKey(path);
    const token = ++this._bakeToken; // supersede any in-flight bake

    const cached = this._cache.get(key);
    if (cached) {
      // refresh LRU recency, bind immediately, drop any overlay a prior bake left up
      this._cache.delete(key); this._cache.set(key, cached);
      this._bindFrames(cached);
      this._hideOverlay();
      return Promise.resolve();
    }

    this._showOverlay();
    const frames = [];
    return this._runBake(path, token, frames, 0).then((done) => {
      if (!done) {
        // a newer setPathology superseded us: discard partial work and let the
        // newer call own the overlay lifecycle.
        Heart3D._disposeFrames(frames);
        return;
      }
      this._cacheStore(key, frames);
      this._bindFrames(frames);
      this._hideOverlay();
    });
  }

  _geomFrom(fn) {
    const g = new THREE.BufferGeometry();
    const mesh = buildMesh(fn, BBOX, CELL);
    if (mesh) {
      g.setAttribute('position', new THREE.BufferAttribute(mesh.positions, 3));
      g.setAttribute('normal', new THREE.BufferAttribute(mesh.normals, 3));
      g.setIndex(new THREE.BufferAttribute(mesh.indices, 1));
      g.computeBoundingSphere();
    }
    return g;
  }

  _buildProbe() {
    this.probeGroup = new THREE.Group();
    // transducer body
    const body = new THREE.Mesh(
      new THREE.CylinderGeometry(0.5, 0.7, 1.4, 20),
      new THREE.MeshStandardMaterial({ color: 0x222831, roughness: 0.5, metalness: 0.3 })
    );
    body.position.y = 0.7;
    this.probeGroup.add(body);
    const face = new THREE.Mesh(
      new THREE.CylinderGeometry(0.5, 0.5, 0.12, 20),
      new THREE.MeshStandardMaterial({ color: 0x39c0e8, emissive: 0x0a4a5a, emissiveIntensity: 0.6, roughness: 0.3 })
    );
    this.probeGroup.add(face);
    this.scene.add(this.probeGroup);

    // imaging fan (triangle sector) + edges.
    // Geometry is rebuilt only when the probe actually moves (see updateProbe);
    // buffers are pre-allocated once here and mutated in place afterwards.
    this.fanSeg = 22;
    const seg = this.fanSeg;
    this.fanMat = new THREE.MeshBasicMaterial({ color: 0x39c0e8, transparent: true, opacity: 0.14, side: THREE.DoubleSide, depthWrite: false });
    this.fanGeo = new THREE.BufferGeometry();
    // one triangle per segment: seg * 3 verts * 3 floats
    this.fanVertArr = new Float32Array(seg * 9);
    this.fanGeo.setAttribute('position', new THREE.BufferAttribute(this.fanVertArr, 3));
    this.fanMesh = new THREE.Mesh(this.fanGeo, this.fanMat);
    this.fanMesh.renderOrder = 6;
    this.scene.add(this.fanMesh);

    this.fanEdgeMat = new THREE.LineBasicMaterial({ color: 0x66e0ff, transparent: true, opacity: 0.7 });
    this.fanEdgeGeo = new THREE.BufferGeometry();
    // 2 side edges + seg arc edges, each a 2-vertex line segment (6 floats)
    this.fanEdgeArr = new Float32Array((2 + seg) * 6);
    this.fanEdgeGeo.setAttribute('position', new THREE.BufferAttribute(this.fanEdgeArr, 3));
    this.fanEdges = new THREE.LineSegments(this.fanEdgeGeo, this.fanEdgeMat);
    this.fanEdges.renderOrder = 7;
    this.scene.add(this.fanEdges);

    // scratch + dirty-flag cache for the fan rebuild
    this.fanRim = new Float32Array((seg + 1) * 3);
    this._fanBd = new THREE.Vector3();
    this._fanTmp = new THREE.Vector3();
    this._lastFan = null; // {pos:[3], dir:[3], lat:[3], half, depth}
  }

  _fanUnchanged(pos, dir, lat, half, depth) {
    const c = this._lastFan;
    if (!c) return false;
    const EPS = 1e-4;
    const near = (a, b) =>
      Math.abs(a[0] - b[0]) < EPS && Math.abs(a[1] - b[1]) < EPS && Math.abs(a[2] - b[2]) < EPS;
    return (
      Math.abs(c.half - half) < EPS &&
      Math.abs(c.depth - depth) < EPS &&
      near(c.pos, [pos.x, pos.y, pos.z]) &&
      near(c.dir, [dir.x, dir.y, dir.z]) &&
      near(c.lat, [lat.x, lat.y, lat.z])
    );
  }

  // Small radial-gradient alpha sprite so each flow point is a soft round
  // streamer rather than a hard square.
  _makeFlowSprite() {
    const s = 64;
    const cv = document.createElement('canvas');
    cv.width = cv.height = s;
    const ctx = cv.getContext('2d');
    const g = ctx.createRadialGradient(s / 2, s / 2, 0, s / 2, s / 2, s / 2);
    g.addColorStop(0.0, 'rgba(255,255,255,1)');
    g.addColorStop(0.35, 'rgba(255,255,255,0.75)');
    g.addColorStop(1.0, 'rgba(255,255,255,0)');
    ctx.fillStyle = g;
    ctx.fillRect(0, 0, s, s);
    const tex = new THREE.CanvasTexture(cv);
    tex.colorSpace = THREE.SRGBColorSpace;
    tex.needsUpdate = true;
    return tex;
  }

  _buildFlow() {
    this.flowCount = 900;
    const geo = new THREE.BufferGeometry();
    this.flowPos = new Float32Array(this.flowCount * 3);
    this.flowCol = new Float32Array(this.flowCount * 3);
    this.flowAlpha = new Float32Array(this.flowCount); // per-particle life fade
    this.flowState = []; // {p:[x,y,z], age}
    for (let i = 0; i < this.flowCount; i++) {
      this.flowState.push({ x: 0, y: 0, z: 0, age: Math.random() * 60 });
      this.flowPos[i * 3 + 1] = -100; // hidden initially
      this.flowAlpha[i] = 0;
    }
    geo.setAttribute('position', new THREE.BufferAttribute(this.flowPos, 3));
    geo.setAttribute('color', new THREE.BufferAttribute(this.flowCol, 3));
    geo.setAttribute('alpha', new THREE.BufferAttribute(this.flowAlpha, 1));
    const sprite = this._makeFlowSprite();
    const mat = new THREE.PointsMaterial({
      size: 0.34, vertexColors: true, transparent: true, opacity: 0.95,
      depthWrite: false, map: sprite, alphaMap: sprite, blending: THREE.NormalBlending,
      sizeAttenuation: true,
    });
    // Feed the per-particle `alpha` attribute into the fragment alpha so
    // particles fade smoothly at the end of their life (not popping out).
    mat.onBeforeCompile = (shader) => {
      shader.vertexShader = 'attribute float alpha;\nvarying float vAlpha;\n' + shader.vertexShader;
      shader.vertexShader = shader.vertexShader.replace(
        '#include <begin_vertex>',
        '#include <begin_vertex>\n  vAlpha = alpha;'
      );
      shader.fragmentShader = 'varying float vAlpha;\n' + shader.fragmentShader;
      shader.fragmentShader = shader.fragmentShader.replace(
        '#include <color_fragment>',
        '#include <color_fragment>\n  diffuseColor.a *= vAlpha;'
      );
    };
    this.flowPoints = new THREE.Points(geo, mat);
    this.scene.add(this.flowPoints);
  }

  updateProbe(probe, sectorHalf, depthCm) {
    // orient/position transducer group
    const pos = new THREE.Vector3(...probe.pos);
    const dir = new THREE.Vector3(...probe.dir).normalize();
    this.probeGroup.position.copy(pos);
    // probe body points opposite the beam (sticking out of the body)
    this.probeGroup.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir.clone().negate());

    const half = sectorHalf, seg = this.fanSeg, depth = depthCm;
    const lat = new THREE.Vector3(...probe.lat).normalize();

    // Dirty-flag: the fan geometry only depends on apex/dir/lat/half/depth.
    // When the probe is static (the common case) skip the whole rebuild.
    if (this._fanUnchanged(pos, dir, lat, half, depth)) return;
    this._lastFan = {
      pos: [pos.x, pos.y, pos.z], dir: [dir.x, dir.y, dir.z],
      lat: [lat.x, lat.y, lat.z], half, depth,
    };

    // Rim points along the far arc, computed once and reused for both the
    // filled triangles and the edge outline.
    const rim = this.fanRim, bd = this._fanBd;
    for (let s = 0; s <= seg; s++) {
      const th = -half + (2 * half) * (s / seg);
      bd.copy(dir).multiplyScalar(Math.cos(th)).addScaledVector(lat, Math.sin(th));
      rim[s * 3] = pos.x + bd.x * depth;
      rim[s * 3 + 1] = pos.y + bd.y * depth;
      rim[s * 3 + 2] = pos.z + bd.z * depth;
    }

    // Filled sector: one triangle (apex, rim[s-1], rim[s]) per segment.
    const va = this.fanVertArr;
    let vi = 0;
    for (let s = 1; s <= seg; s++) {
      va[vi++] = pos.x; va[vi++] = pos.y; va[vi++] = pos.z;
      va[vi++] = rim[(s - 1) * 3]; va[vi++] = rim[(s - 1) * 3 + 1]; va[vi++] = rim[(s - 1) * 3 + 2];
      va[vi++] = rim[s * 3]; va[vi++] = rim[s * 3 + 1]; va[vi++] = rim[s * 3 + 2];
    }

    // Edge outline: two side rays (apex->near/far edge) + the far arc.
    // Flat fan uses MeshBasicMaterial (normals ignored), so no computeVertexNormals.
    const ea = this.fanEdgeArr;
    let ei = 0;
    // left side edge
    ea[ei++] = pos.x; ea[ei++] = pos.y; ea[ei++] = pos.z;
    ea[ei++] = rim[0]; ea[ei++] = rim[1]; ea[ei++] = rim[2];
    // right side edge
    ea[ei++] = pos.x; ea[ei++] = pos.y; ea[ei++] = pos.z;
    ea[ei++] = rim[seg * 3]; ea[ei++] = rim[seg * 3 + 1]; ea[ei++] = rim[seg * 3 + 2];
    // far arc
    for (let s = 0; s < seg; s++) {
      ea[ei++] = rim[s * 3]; ea[ei++] = rim[s * 3 + 1]; ea[ei++] = rim[s * 3 + 2];
      ea[ei++] = rim[(s + 1) * 3]; ea[ei++] = rim[(s + 1) * 3 + 1]; ea[ei++] = rim[(s + 1) * 3 + 2];
    }

    this.fanGeo.attributes.position.needsUpdate = true;
    this.fanGeo.computeBoundingSphere();
    this.fanEdgeGeo.attributes.position.needsUpdate = true;
    this.fanEdgeGeo.computeBoundingSphere();
  }

  // A tapered muscular strut (papillary muscle / moderator band) between two
  // points, base slightly wider than the tip.
  _strut(a, b, r, mat) {
    const V = THREE.Vector3;
    const pa = new V(a[0], a[1], a[2]), pb = new V(b[0], b[1], b[2]);
    const dir = new V().subVectors(pb, pa);
    const len = dir.length() || 1e-3;
    const geo = new THREE.CylinderGeometry(r * 0.7, r, len, 12, 1);
    const mesh = new THREE.Mesh(geo, mat);
    mesh.position.copy(pa).addScaledVector(dir, 0.5);
    mesh.quaternion.setFromUnitVectors(new V(0, 1, 0), dir.clone().normalize());
    return mesh;
  }

  // Append a chordal fan (tip → a ring of leaflet-edge points) into `pos`.
  _fan(tip, ringC, t1, t2, rho, NC, pos) {
    const V = THREE.Vector3;
    for (let i = 0; i < NC; i++) {
      const th = (i / NC) * 2 * Math.PI;
      const r = ringC.clone().addScaledVector(t1, rho * Math.cos(th)).addScaledVector(t2, rho * Math.sin(th));
      pos.push(tip.x, tip.y, tip.z, r.x, r.y, r.z);
    }
  }

  // Build the subvalvular apparatus from the shared anatomy A: the LV papillary
  // muscles + RV anterior papillary (fused to the moderator band), and the
  // chordae-tendineae fans up to the mitral & tricuspid leaflet free-edge rings.
  // Rebuilt on pathology change so it tracks LV dilatation.
  _buildSubvalvular(A) {
    const V = THREE.Vector3;
    if (this.subValv) {                                   // dispose the old group
      this.heartGroup.remove(this.subValv);
      this.subValv.traverse((o) => { if (o.geometry) o.geometry.dispose(); });
    }
    const grp = new THREE.Group();
    const papMat = new THREE.MeshStandardMaterial({ color: 0xb0504a, roughness: 0.75, transparent: true, opacity: 0.92, depthWrite: false });
    const chordaMat = new THREE.LineBasicMaterial({ color: 0xd8cfbb, transparent: true, opacity: 0.4, depthWrite: false });
    const ringBasis = (n) => {
      let t1 = new V().crossVectors(n, new V(0, 0, 1));
      if (t1.lengthSq() < 1e-4) t1 = new V().crossVectors(n, new V(1, 0, 0));
      t1.normalize();
      return { t1, t2: new V().crossVectors(n, t1).normalize() };
    };
    const pos = [];
    // mitral: both LV papillary muscles → mitral leaflet ring (hanging into the LV)
    const mv = VALVE_DEFS.mitral, mn = new V(mv.n[0], mv.n[1], mv.n[2]).normalize();
    const mB = ringBasis(mn), mRing = new V(mv.c[0], mv.c[1], mv.c[2]).addScaledVector(mn, -0.55);
    for (const p of A.pap) {
      grp.add(this._strut(p.a, p.b, p.r, papMat));
      this._fan(new V(p.b[0], p.b[1], p.b[2]), mRing, mB.t1, mB.t2, 0.95, 9, pos);
    }
    // RV anterior papillary + moderator band → tricuspid leaflet ring
    if (A.rvPap) {
      grp.add(this._strut(A.rvPap.a, A.rvPap.b, A.rvPap.r, papMat));
      if (A.mod) grp.add(this._strut(A.mod.a, A.mod.b, A.mod.r, papMat));
      const tv = VALVE_DEFS.tricuspid, tn = new V(tv.n[0], tv.n[1], tv.n[2]).normalize();
      const tB = ringBasis(tn), tRing = new V(tv.c[0], tv.c[1], tv.c[2]).addScaledVector(tn, -0.5);
      this._fan(new V(A.rvPap.b[0], A.rvPap.b[1], A.rvPap.b[2]), tRing, tB.t1, tB.t2, 1.0, 9, pos);
    }
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
    const chordae = new THREE.LineSegments(g, chordaMat);
    chordae.renderOrder = 4;
    grp.add(chordae);
    grp.renderOrder = 4;
    this.subValv = grp;
    this.heartGroup.add(grp);
  }

  // Toggle the myofiber overlay; builds it lazily on first show.
  showFibers(on) {
    this.fiberLines.visible = !!on;
    if (on && this._fiberKey === null) this._buildFibers(geometryAt(0.5, this._fiberPath || {}));
  }

  // Build the helical fibre field from the LV ellipsoid proxy. Two transmural
  // shells (sub-endo, sub-epi) are drawn as short segments whose direction is the
  // circumferential vector rotated by the transmural helix angle — the counter-wound
  // helices that give the LV its wringing (torsion) contraction. Coloured endo→epi.
  _buildFibers(G) {
    const lv = G.lv, c = lv.c, r = lv.r, wall = lv.wall;
    const V = THREE.Vector3;
    const pos = [], col = [];
    const NU = 30, NV = 9, LEN = 0.42;               // azimuth / polar samples, segment length
    const shells = [{ d: 0.22, ang: 55 }, { d: 0.82, ang: -55 }]; // endo +, epi −
    const endoCol = new THREE.Color(0xff7a4d), epiCol = new THREE.Color(0x4d9dff);
    for (const sh of shells) {
      const R = [r[0] + wall * sh.d, r[1] + wall * sh.d, r[2] + wall * sh.d];
      const a = sh.ang * Math.PI / 180;
      const tone = sh.d < 0.5 ? endoCol : epiCol;
      for (let iv = 1; iv < NV; iv++) {
        const v = (0.16 + 0.74 * iv / NV) * Math.PI;  // avoid apex/base poles
        for (let iu = 0; iu < NU; iu++) {
          const u = (iu / NU) * 2 * Math.PI;
          const sv = Math.sin(v), cv = Math.cos(v);
          const p = new V(c[0] + R[0] * sv * Math.cos(u), c[1] + R[1] * cv, c[2] + R[2] * sv * Math.sin(u));
          // outward normal = ellipsoid gradient
          const n = new V((p.x - c[0]) / (R[0] * R[0]), (p.y - c[1]) / (R[1] * R[1]), (p.z - c[2]) / (R[2] * R[2])).normalize();
          // longitudinal = long axis (+y) projected into the tangent plane
          const lng = new V(0, 1, 0).addScaledVector(n, -n.y).normalize();
          const circ = new V().crossVectors(n, lng).normalize();
          // fibre = circumferential rotated toward longitudinal by the helix angle
          const fib = circ.clone().multiplyScalar(Math.cos(a)).addScaledVector(lng, Math.sin(a)).normalize();
          const p0 = p.clone().addScaledVector(fib, -LEN / 2);
          const p1 = p.clone().addScaledVector(fib, LEN / 2);
          pos.push(p0.x, p0.y, p0.z, p1.x, p1.y, p1.z);
          col.push(tone.r, tone.g, tone.b, tone.r, tone.g, tone.b);
        }
      }
    }
    const g = this.fiberLines.geometry;
    g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
    g.setAttribute('color', new THREE.Float32BufferAttribute(col, 3));
    g.computeBoundingSphere();
    this._fiberKey = 'built';
  }

  update(G, path, probe, dt) {
    // beating loop: bind the pre-baked marching-cubes geometry for this phase.
    // During the deferred startup bake not every phase exists yet, so fall back
    // to phase 0 for any phase that hasn't been baked. Track the bound frame
    // object (not the index) so the fallback doesn't rebind geometry needlessly.
    if (this.frames && this.frames.length) {
      const idx = ((Math.round((G.phase % 1) * NPHASE) % NPHASE) + NPHASE) % NPHASE;
      const f = this.frames[idx] || this.frames[0];
      if (f && f !== this._boundFrame) {
        this.myoMesh.geometry = f.myo;
        this.lvMesh.geometry = f.LV; this.rvMesh.geometry = f.RV;
        this.laMesh.geometry = f.LA; this.raMesh.geometry = f.RA;
        this.vesMesh.geometry = f.VES;
        this._boundFrame = f;
      }
    }
    // keep the fibre overlay sized to the current pathology's LV (rebuild on change)
    if (this.fiberLines.visible && path !== this._fiberPath) {
      this._fiberPath = path;
      this._buildFibers(geometryAt(0.5, path || {}));
    }

    // transparent walls: see chambers through the muscle; solid: opaque muscle
    this.myoMesh.material.opacity = this.transparent ? 0.42 : 0.98;
    this.myoMesh.material.depthWrite = !this.transparent;
    const bloodVis = this.transparent;
    this.bloodMeshes.forEach((m) => (m.visible = bloodVis));

    // valves: cue opening/closing via opacity (a shut valve reads as a bright
    // annulus; an open valve fades as the leaflets swing to the wall). Also fade
    // toward transparent when the ring is viewed edge-on (its normal ~⊥ view),
    // where a flat annulus would otherwise degenerate to an ugly bright line.
    const camPos = this.camera.position;
    for (const name in this.valveMeshes) {
      const m = this.valveMeshes[name];
      const open = G.valves[name];
      const base = 0.2 + 0.4 * (1 - open);
      const n = m.userData.normal;
      const vx = camPos.x - m.position.x, vy = camPos.y - m.position.y, vz = camPos.z - m.position.z;
      const vl = Math.hypot(vx, vy, vz) || 1;
      const facing = Math.abs((n.x * vx + n.y * vy + n.z * vz) / vl); // 1 face-on, 0 edge-on
      const edge = THREE.MathUtils.smoothstep(facing, 0.05, 0.4);
      m.material.opacity = base * (0.1 + 0.9 * edge);
    }

    // flow particles
    if (this.showFlow) {
      this.flowPoints.visible = true;
      this._updateFlow(G, path, probe, dt);
    } else {
      this.flowPoints.visible = false;
    }

    this.controls.update();
    this.renderer.render(this.scene, this.camera);
  }

  _seedFlow(st) {
    // seed a particle at an inlet region depending on cycle handled generically:
    // pick a random chamber-ish location; velocity field will carry it.
    const mid = (a, b, t) => [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t];
    const spots = [
      mid(LM.M, LM.LA_ROOF, 0.4), mid(LM.M, LM.apex, 0.35),           // LA -> LV column
      mid(LM.T, LM.RA_ROOF, 0.4), mid(LM.T, [-2.6, -5.5, 1.2], 0.35), // RA -> RV
      LM.LVOT0, mid(LM.A, LM.AO_STJ, 0.6),                             // LVOT / aorta
    ];
    const s = spots[(Math.random() * spots.length) | 0];
    st.x = s[0] + (Math.random() - 0.5) * 1.6;
    st.y = s[1] + (Math.random() - 0.5) * 1.6;
    st.z = s[2] + (Math.random() - 0.5) * 1.2;
    st.age = 0;
  }

  _updateFlow(G, path, probe, dt) {
    const probePos = probe.pos;
    // Frame-time budget: the per-particle velocityAt() sampling dominates cost,
    // so on slow frames we shrink the active particle set and grow it back when
    // frames are cheap. Degrades gracefully on weak hardware without stutter.
    if (this.flowActive == null) this.flowActive = this.flowCount;
    if (dt > 0.040) this.flowActive = Math.max(200, this.flowActive - 40);
    else if (dt < 0.022) this.flowActive = Math.min(this.flowCount, this.flowActive + 25);
    for (let i = 0; i < this.flowCount; i++) {
      if (i >= this.flowActive) {
        // parked: skip the expensive velocity sampling, just keep it hidden
        if (this.flowPos[i * 3 + 1] !== -100) this.flowPos[i * 3 + 1] = -100;
        this.flowAlpha[i] = 0;
        continue;
      }
      const st = this.flowState[i];
      const vel = velocityAt(st.x, st.y, st.z, G, path);
      // Cull: with flow now confined to jets/tracts, most of the heart has no
      // meaningful velocity. Hide stagnant/out-of-flow particles (no floating
      // cluster near the probe) and reseed them so they can re-enter a stream.
      const outOfBounds = Math.abs(st.x) > 7 || st.y > 6 || st.y < -7 || Math.abs(st.z) > 6;
      if (!vel || vel.speed < 0.07 || st.age > 55 || outOfBounds) {
        this.flowPos[i * 3 + 1] = -100; // hidden
        this.flowAlpha[i] = 0;
        st.age += 2;
        if (st.age > 55 + (i % 20)) this._seedFlow(st);
        continue;
      }
      if (vel) {
        const sp = Math.min(vel.speed, 5);
        st.x += vel.vx * dt * 3.2 / Math.max(vel.speed, 0.001) * sp;
        st.y += vel.vy * dt * 3.2 / Math.max(vel.speed, 0.001) * sp;
        st.z += vel.vz * dt * 3.2 / Math.max(vel.speed, 0.001) * sp;
        st.age += dt * 12;
        this.flowPos[i * 3] = st.x;
        this.flowPos[i * 3 + 1] = st.y;
        this.flowPos[i * 3 + 2] = st.z;
        // colour by line-of-sight relative to probe (Doppler convention)
        const bx = st.x - probePos[0], by = st.y - probePos[1], bz = st.z - probePos[2];
        const bl = Math.hypot(bx, by, bz) || 1;
        const los = (vel.vx * bx + vel.vy * by + vel.vz * bz) / bl; // + = away
        const mag = Math.min(Math.abs(los) / 1.2, 1);
        if (los < 0) { // toward probe -> red
          this.flowCol[i * 3] = 0.8 + 0.2 * mag; this.flowCol[i * 3 + 1] = 0.25 * mag; this.flowCol[i * 3 + 2] = 0.15;
        } else {
          this.flowCol[i * 3] = 0.12; this.flowCol[i * 3 + 1] = 0.35 + 0.4 * mag; this.flowCol[i * 3 + 2] = 0.85;
        }
        if (vel.turbulent) { this.flowCol[i * 3 + 1] = 0.85; }
        // life fade: ramp up briefly after seeding, fade out over the last
        // third of life so streamers dissolve instead of blinking off.
        const life = st.age / 55;
        let a = life > 0.66 ? Math.max(0, (1 - life) / 0.34) : 1;
        a *= Math.min(1, 0.25 + st.age * 0.35);
        this.flowAlpha[i] = a;
      }
    }
    this.flowPoints.geometry.attributes.position.needsUpdate = true;
    this.flowPoints.geometry.attributes.color.needsUpdate = true;
    this.flowPoints.geometry.attributes.alpha.needsUpdate = true;
  }

  setCameraForView(name) {
    // Per-view direction, but distance is derived from the heart's bounding
    // sphere so every view is framed consistently (never oversized/clipped).
    const dirs = {
      PLAX: [12, 2, 6], A4C: [2, -3, 15], PSAX: [3, 12, 4], PSAX_AV: [3, 10, 5], SUBCOSTAL: [4, -8, 11], A2C: [13, 1, 3],
    };
    const dir = new THREE.Vector3(...(dirs[name] || [11, 4, 13])).normalize();
    const box = new THREE.Box3().setFromObject(this.heartGroup);
    if (box.isEmpty()) { // no geometry yet — fall back to a fixed frame
      this.camera.position.copy(dir.multiplyScalar(18));
      this.controls.target.set(0, 0.2, 0);
      return;
    }
    const sph = box.getBoundingSphere(new THREE.Sphere());
    const fov = THREE.MathUtils.degToRad(this.camera.fov);
    const margin = 1.25;
    const dist = (sph.radius / Math.sin(fov / 2)) * margin;
    this.controls.target.copy(sph.center);
    this.camera.position.copy(sph.center).addScaledVector(dir, dist);
    this.camera.lookAt(sph.center);
  }

  _resize() {
    const w = this.container.clientWidth, h = this.container.clientHeight;
    if (!w || !h) return;
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
    this.renderer.setSize(w, h);
  }
}
