// tools/verify-feedback.mjs — checks for the bug report / feature request dialog.
//
// Part 1 (Node): calls api/feedback.js directly with a stubbed mail API and
// checks validation, the honeypot, the rate limit, the mailto fallback when no
// mail key is set, and the email that is sent (recipient, subject, reply-to,
// diagnostics).
// Part 2 (browser): opens the dialog from the header and, with /api/feedback
// stubbed, checks required-field validation, the bug/feature tabs, the
// diagnostics preview, success, the mailto fallback, the clipboard fallback
// when there is no backend, Escape and focus return, and no console errors.
//
// Usage: node tools/verify-feedback.mjs      (exit 0 = pass, 1 = fail)
import http from 'http';
import fs from 'fs';
import path from 'path';
import { createRequire } from 'module';
import { fileURLToPath } from 'url';

const results = [];
const check = (name, ok, note = '') => results.push({ name, ok: !!ok, note });
const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');

// ---- part 1: the serverless function ----------------------------------------------
const handler = createRequire(import.meta.url)('../api/feedback.js');
const call = async (body, { method = 'POST', ip = '1.1.1.1' } = {}) => {
  const out = { status: 0, json: null, headers: {} };
  const res = {
    status(c) { out.status = c; return this; },
    json(j) { out.json = j; return this; },
    setHeader(k, v) { out.headers[k] = v; },
  };
  await handler({ method, body, headers: { 'x-forwarded-for': ip } }, res);
  return out;
};
const sent = [];
const realFetch = globalThis.fetch;
globalThis.fetch = async (url, opts) => { sent.push({ url, opts, body: JSON.parse(opts.body) }); return { ok: true, status: 200, text: async () => '' }; };
const realErr = console.error; console.error = () => {};

const OK = { type: 'bug', summary: 'Doppler freezes', details: 'Steps…', email: 'learner@example.org', diagnostics: 'View: PLAX' };
delete process.env.FEEDBACK_TO; delete process.env.RESEND_API_KEY;
let r = await call(OK, { method: 'GET' });
check('GET is rejected with 405', r.status === 405 && r.headers.Allow === 'POST');
r = await call(OK, { ip: '2.0.0.1' });
check('no FEEDBACK_TO: 503, nothing sent', r.status === 503 && !r.json.to && sent.length === 0);
process.env.FEEDBACK_TO = 'owner@example.com';
r = await call(OK, { ip: '2.0.0.2' });
check('no mail key: 503 with mailto fallback to FEEDBACK_TO', r.status === 503 && r.json.fallback === 'mailto' && r.json.to === 'owner@example.com' && sent.length === 0);
process.env.RESEND_API_KEY = 'test-key';
r = await call({ ...OK, summary: '   ' }, { ip: '2.0.0.3' });
check('empty summary: 400', r.status === 400);
r = await call({ ...OK, email: 'not an email' }, { ip: '2.0.0.4' });
check('bad reply-to: 400', r.status === 400);
r = await call('{bad json', { ip: '2.0.0.5' });
check('malformed body: 400', r.status === 400);
r = await call({ ...OK, website: 'http://spam' }, { ip: '2.0.0.6' });
check('honeypot: fake success, nothing sent', r.status === 200 && r.json.ok && sent.length === 0);
r = await call({ ...OK, summary: 'Line one\nBcc: x@y.z' }, { ip: '2.0.0.7' });
const m = sent[0];
check('valid bug: 200 and one email sent', r.status === 200 && r.json.ok && sent.length === 1);
check('email goes to FEEDBACK_TO via Resend', m && m.url === 'https://api.resend.com/emails' && m.body.to[0] === 'owner@example.com' && /Bearer test-key/.test(m.opts.headers.Authorization));
check('subject is tagged and single-line', m && m.body.subject === '[EchoSim bug] Line one Bcc: x@y.z');
check('reply-to is the learner', m && m.body.reply_to === 'learner@example.org');
check('body carries details and diagnostics', m && /Steps…/.test(m.body.text) && /View: PLAX/.test(m.body.text));
await call({ ...OK, type: 'feature', email: '' }, { ip: '2.0.0.8' });
check('feature request: tagged, no reply-to', sent[1] && /^\[EchoSim feature\]/.test(sent[1].body.subject) && !('reply_to' in sent[1].body));
await call({ ...OK, details: 'x'.repeat(20000) }, { ip: '2.0.0.9' });
check('long details are clipped', sent[2] && sent[2].body.text.length < 9000);
let codes = [];
for (let i = 0; i < 7; i++) codes.push((await call(OK, { ip: '9.9.9.9' })).status);
check('rate limit: 6th message in 10 min gets 429', codes.slice(0, 5).every((c) => c === 200) && codes[5] === 429, codes.join(','));
globalThis.fetch = async () => ({ ok: false, status: 422, text: async () => 'bad' });
r = await call(OK, { ip: '2.0.1.0' });
check('mail API failure: 502', r.status === 502 && !r.json.ok);
globalThis.fetch = realFetch; console.error = realErr;

// ---- part 2: the dialog in the browser ---------------------------------------------
let chromium;
try { ({ chromium } = await import('playwright')); } catch (e) { chromium = null; }
if (!chromium) {
  check('browser part ran (playwright importable)', false, 'playwright not found: npm install playwright');
} else {
  const MIME = { '.html': 'text/html', '.js': 'text/javascript', '.css': 'text/css', '.svg': 'image/svg+xml', '.png': 'image/png', '.webmanifest': 'application/manifest+json' };
  const PORT = +(process.env.PORT || 8097);
  const server = http.createServer((req, res) => {
    let p = decodeURIComponent(req.url.split('?')[0]);
    if (p === '/') p = '/index.html';
    fs.readFile(path.join(ROOT, p), (e, d) => {
      if (e) { res.writeHead(404); res.end('not found'); return; }
      res.writeHead(200, { 'Content-Type': MIME[path.extname(p)] || 'application/octet-stream' });
      res.end(d);
    });
  });
  await new Promise((ok) => server.listen(PORT, ok));
  const browser = await chromium.launch({
    executablePath: process.env.CHROME_PATH || [
      '/opt/pw-browsers/chromium', '/opt/pw-browsers/chromium-1194/chrome-linux/chrome',
    ].find((p) => fs.existsSync(p)),
    args: ['--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader', '--no-sandbox'],
  });
  const ctx = await browser.newContext({ viewport: { width: 1440, height: 880 }, permissions: ['clipboard-read', 'clipboard-write'] });
  const page = await ctx.newPage();
  const errors = [];
  page.on('console', (msg) => { if (msg.type() === 'error' && !/favicon|api\/feedback|status of 404|status of 503/.test(msg.text())) errors.push(msg.text()); });
  page.on('pageerror', (e) => errors.push(e.message));
  let mode = 'ok', posted = null;
  await page.route('**/api/feedback', async (route) => {
    posted = JSON.parse(route.request().postData() || '{}');
    if (mode === 'ok') return route.fulfill({ status: 200, contentType: 'application/json', body: '{"ok":true}' });
    if (mode === 'mailto') return route.fulfill({ status: 503, contentType: 'application/json', body: '{"ok":false,"fallback":"mailto","to":"owner@example.com"}' });
    return route.fulfill({ status: 404, body: 'not found' });
  });
  await page.goto(`http://localhost:${PORT}/`, { waitUntil: 'networkidle' });
  await page.waitForTimeout(1200);
  const got = await page.$('#helpGot'); if (got) await got.click().catch(() => {});
  await page.waitForTimeout(300);

  const status = () => page.$eval('#fbStatus', (el) => el.textContent);
  await page.selectOption('#pathology', 'as');
  await page.click('#feedbackBtn');
  check('dialog opens with focus in the summary', await page.evaluate(() => !document.getElementById('feedbackOverlay').classList.contains('hidden') && document.activeElement.id === 'fbSummary'));
  const diag = await page.$eval('#fbDiagText', (el) => el.textContent);
  check('diagnostics preview shows the case and browser', /Pathology: as/.test(diag) && /View: /.test(diag) && /Browser: /.test(diag) && /WebGL: /.test(diag), diag.split('\n')[0]);
  const EMAILISH = /[\w.+-]+@[\w-]+\.[\w.]+/;
  check('no recipient address in the shipped client code', ['index.html', 'js/feedback.js', 'js/main.js'].every((f) => !EMAILISH.test(fs.readFileSync(path.join(ROOT, f), 'utf8'))));
  await page.keyboard.press('Space');
  check('typing in the dialog does not trigger shortcuts', await page.$eval('#fbSummary', (el) => el.value === ' '));
  await page.fill('#fbSummary', '');
  await page.click('#fbSend');
  check('empty summary is refused before sending', posted === null && (await page.$eval('#fbSummary', (el) => el.getAttribute('aria-invalid'))) === 'true');
  await page.click('#ftab-feature');
  check('feature tab relabels the form and unticks diagnostics', /would you like/.test(await page.$eval('#fbDetailsLabel', (el) => el.textContent)) && !(await page.$eval('#fbDiag', (el) => el.checked)));
  await page.click('#ftab-bug');
  await page.fill('#fbSummary', 'Colour Doppler flickers in A4C');
  await page.fill('#fbDetails', '1. Open A4C\n2. Turn on colour');
  await page.fill('#fbEmail', 'learner@example.org');
  await page.click('#fbSend');
  await page.waitForFunction(() => /Sent/.test(document.getElementById('fbStatus').textContent), null, { timeout: 5000 });
  check('successful send posts the form', posted && posted.type === 'bug' && posted.summary === 'Colour Doppler flickers in A4C' && posted.email === 'learner@example.org' && /Pathology: as/.test(posted.diagnostics) && posted.website === '');
  await page.waitForTimeout(1900);
  check('dialog closes after sending and focus returns', await page.evaluate(() => document.getElementById('feedbackOverlay').classList.contains('hidden') && document.activeElement.id === 'feedbackBtn'));

  mode = 'mailto'; posted = null;
  await page.click('#feedbackBtn');
  await page.fill('#fbSummary', 'Add a PLAX zoom view');
  await page.click('#ftab-feature');
  await page.click('#fbSend');
  await page.waitForTimeout(600);
  check('no mail key: falls back to the user\'s mail app', /email app/.test(await status()) && posted && posted.type === 'feature' && posted.diagnostics === '');

  mode = 'none';
  await page.fill('#fbSummary', 'Offline message');
  await page.click('#fbSend');
  await page.waitForTimeout(600);
  const clip = await page.evaluate(() => navigator.clipboard.readText().catch(() => ''));
  check('no backend: message kept on the clipboard', /clipboard/.test(await status()) && /Offline message/.test(clip));
  await page.keyboard.press('Escape');
  check('Escape closes the dialog', await page.evaluate(() => document.getElementById('feedbackOverlay').classList.contains('hidden')));
  if (process.env.SHOT) {
    await page.click('#feedbackBtn'); await page.click('#ftab-bug');
    await page.fill('#fbSummary', 'Colour Doppler flickers in A4C');
    await page.fill('#fbDetails', 'Steps to reproduce:\n1. Open A4C\n2. Turn on colour');
    await page.click('.fb-diag summary');
    await page.screenshot({ path: process.env.SHOT });
  }
  check('no console errors', errors.length === 0, errors.slice(0, 3).join(' | '));
  await browser.close(); server.close();
}

for (const x of results) console.log(`  ${x.ok ? 'PASS' : 'FAIL'} ${x.name}${x.note ? '  [' + x.note + ']' : ''}`);
const failed = results.filter((x) => !x.ok).length;
console.log(`\n  ${results.length - failed} passed, ${failed} failed.`);
process.exit(failed ? 1 : 0);
