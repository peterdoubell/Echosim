// One elapsed-time clock for the ECG samples and mechanical animation. Slow
// frames may omit intermediate images, but must never slow the heart rate.
export class PlaybackClock {
  constructor(heart, onSample, sampleRate = 250, historySeconds = 6) {
    this.heart = heart;
    this.onSample = onSample;
    this.step = 1 / sampleRate;
    this.maxSamples = Math.ceil(sampleRate * historySeconds);
    this.remainder = 0;
    this.time = 0;
  }

  advance(seconds) {
    if (!Number.isFinite(seconds) || seconds <= 0) return;
    this.time += seconds;
    this.remainder += seconds;
    const count = Math.floor((this.remainder + 1e-10) / this.step);
    // Following a long foreground stall, only the visible ECG history needs
    // sampling. Advance older beats directly, keeping their rhythm schedule.
    const skipped = Math.max(0, count - this.maxSamples);
    if (skipped) this.heart.advance(skipped * this.step);
    for (let i = skipped; i < count; i++) {
      this.heart.advance(this.step);
      this.heart.leads();
      this.onSample(this.heart.vcg);
    }
    this.remainder = Math.max(0, this.remainder - count * this.step);
  }
}
