// Vercel serverless function: emails a bug report / feature request from the
// in-app Feedback dialog (js/feedback.js) to the maintainer.
//
// Environment (Vercel project settings → Environment Variables):
//   FEEDBACK_TO      recipient address (required)
//   RESEND_API_KEY   Resend API key (https://resend.com); without it the
//                    function answers 503 with a mailto fallback so the dialog
//                    opens the user's own mail app instead
//   FEEDBACK_FROM    optional sender, e.g. "EchoSim <feedback@your-domain>";
//                    defaults to Resend's shared onboarding sender, which can
//                    only deliver to the Resend account owner's own address

const LIMITS = { summary: 200, details: 5000, email: 200, diagnostics: 4000 };
const WINDOW_MS = 10 * 60 * 1000, MAX_PER_WINDOW = 5;
const hits = new Map(); // best-effort per-instance rate limit: ip -> [timestamps]

const clip = (v, n) => (typeof v === 'string' ? v.trim().slice(0, n) : '');
const EMAIL_RE = /^[^\s@<>()",;]+@[^\s@<>()",;]+\.[^\s@<>()",;]+$/;

function rateLimited(ip) {
  const now = Date.now();
  const list = (hits.get(ip) || []).filter((t) => now - t < WINDOW_MS);
  list.push(now);
  hits.set(ip, list);
  if (hits.size > 5000) hits.clear();
  return list.length > MAX_PER_WINDOW;
}

module.exports = async function handler(req, res) {
  if (req.method !== 'POST') {
    res.setHeader('Allow', 'POST');
    return res.status(405).json({ ok: false, error: 'Method not allowed' });
  }
  let body = req.body;
  if (typeof body === 'string') { try { body = JSON.parse(body); } catch (e) { body = null; } }
  if (!body || typeof body !== 'object') return res.status(400).json({ ok: false, error: 'Invalid request' });

  // honeypot filled: pretend success so bots learn nothing
  if (body.website) return res.status(200).json({ ok: true });

  const type = body.type === 'feature' ? 'feature' : 'bug';
  const summary = clip(body.summary, LIMITS.summary).replace(/[\r\n]+/g, ' ');
  const details = clip(body.details, LIMITS.details);
  const email = clip(body.email, LIMITS.email);
  const diagnostics = clip(body.diagnostics, LIMITS.diagnostics);
  if (!summary) return res.status(400).json({ ok: false, error: 'Please add a short summary.' });
  if (email && !EMAIL_RE.test(email)) return res.status(400).json({ ok: false, error: 'That email address does not look right.' });

  const ip = String(req.headers['x-forwarded-for'] || '').split(',')[0].trim() || 'unknown';
  if (rateLimited(ip)) return res.status(429).json({ ok: false, error: 'Too many messages — please try again later.' });

  const to = process.env.FEEDBACK_TO;
  const key = process.env.RESEND_API_KEY;
  if (!to) return res.status(503).json({ ok: false, error: 'Feedback is not configured on this deployment.' });
  if (!key) return res.status(503).json({ ok: false, fallback: 'mailto', to });

  const label = type === 'bug' ? 'Bug report' : 'Feature request';
  const commit = (process.env.VERCEL_GIT_COMMIT_SHA || '').slice(0, 7);
  const text = [
    `${label} from EchoSim${commit ? ` (build ${commit})` : ''}`,
    '',
    `Summary: ${summary}`,
    `Reply to: ${email || '(not given)'}`,
    '',
    details || '(no details)',
    '',
    diagnostics ? `--- Diagnostics ---\n${diagnostics}` : '(diagnostics not included)',
  ].join('\n');

  try {
    const r = await fetch('https://api.resend.com/emails', {
      method: 'POST',
      headers: { Authorization: `Bearer ${key}`, 'Content-Type': 'application/json' },
      body: JSON.stringify({
        from: process.env.FEEDBACK_FROM || 'EchoSim Feedback <onboarding@resend.dev>',
        to: [to],
        subject: `[EchoSim ${type}] ${summary}`,
        text,
        ...(email ? { reply_to: email } : {}),
      }),
    });
    if (!r.ok) {
      console.error('Resend error', r.status, await r.text().catch(() => ''));
      return res.status(502).json({ ok: false, error: 'The mail service rejected the message. Please try again later.' });
    }
    return res.status(200).json({ ok: true });
  } catch (e) {
    console.error('Feedback send failed', e);
    return res.status(502).json({ ok: false, error: 'Could not reach the mail service. Please try again later.' });
  }
};
