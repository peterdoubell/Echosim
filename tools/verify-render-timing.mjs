// Regression checks for render-independent cardiac playback. Run with Node;
// no browser or timing-sensitive wall-clock assertions are required.
import assert from 'node:assert/strict';
import { Heartbeat } from '../js/electrophysiology.js';
import { PlaybackClock } from '../js/playback-clock.js';

function sample(frames, hr = 72, rhythm = 'sinus') {
  const heart = new Heartbeat();
  heart.setHR(hr);
  heart.setRhythm(rhythm);
  const history = [];
  const clock = new PlaybackClock(heart, (vcg) => {
    history.push(vcg.slice());
    if (history.length > 1500) history.shift();
  });
  for (const dt of frames) clock.advance(dt);
  return { heart, clock, history };
}

for (const hr of [50, 72, 120, 160]) {
  for (const rhythm of ['sinus', 'afib', 'pvc']) {
    const fast = sample(Array(1200).fill(1 / 60), hr, rhythm);
    // 80 ms per frame previously advanced only 50 ms: 20 seconds of playback
    // fell 7.5 seconds behind the displayed heart rate.
    const slow = sample(Array(250).fill(0.08), hr, rhythm);
    assert.ok(Math.abs(slow.heart._t - 20) < 1e-8, `${hr} ${rhythm}: no elapsed time discarded`);
    assert.ok(Math.abs(slow.heart.mechanical().phase - fast.heart.mechanical().phase) < 1e-8,
      `${hr} ${rhythm}: mechanical phase independent of render cadence`);
    assert.deepEqual(slow.history, fast.history, `${hr} ${rhythm}: same ECG samples at 60 and 12.5 fps`);
  }
}

{
  const normal = sample(Array(3000).fill(0.004));
  const stall = sample([12]);
  assert.equal(stall.history.length, 1500, 'long stall only samples visible six-second ECG history');
  assert.ok(Math.abs(stall.heart._t - 12) < 1e-8, 'long stall retains elapsed time');
  assert.ok(Math.abs(stall.heart.mechanical().phase - normal.heart.mechanical().phase) < 1e-8,
    'bulk catch-up preserves the beat schedule');
  const before = stall.heart._t;
  for (const dt of [0, -1, NaN, Infinity]) stall.clock.advance(dt);
  assert.equal(stall.heart._t, before, 'paused/invalid deltas do not advance playback');
}

console.log('PASS: ECG and echo phase retain real time across 50–160 bpm, sinus/AF/PVC, slow frames and long stalls.');
