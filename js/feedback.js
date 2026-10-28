// Bug report / feature request dialog.
//
// POSTs to api/feedback.js (a Vercel function), which emails the maintainer.
// The recipient lives only in the deployment's FEEDBACK_TO variable, never in
// this file. If the function has no mail key configured it answers 503 with a
// mailto fallback, and the message opens in the user's own mail app instead;
// if there is no function at all (a local static server) the message is copied
// to the clipboard so nothing typed is lost.

const MAX_ERRORS = 10;
const recentErrors = [];
const note = (msg) => {
  recentErrors.push(`${new Date().toISOString().slice(11, 19)} ${String(msg).slice(0, 300)}`);
  if (recentErrors.length > MAX_ERRORS) recentErrors.shift();
};
// start collecting as soon as the module loads
window.addEventListener('error', (e) => note(`${e.message} (${(e.filename || '').split('/').pop()}:${e.lineno})`));
window.addEventListener('unhandledrejection', (e) => note(`Unhandled rejection: ${e.reason && e.reason.message || e.reason}`));

const DETAILS = {
  bug: { label: 'What happened?', placeholder: 'Steps to reproduce:\n1. \n2. \n\nWhat you expected:\n\nWhat happened instead:' },
  feature: { label: 'What would you like, and why?', placeholder: 'Describe the feature and how it would help your training or teaching.' },
};

function webglRenderer() {
  try {
    const gl = document.createElement('canvas').getContext('webgl');
    if (!gl) return 'none';
    const ext = gl.getExtension('WEBGL_debug_renderer_info');
    return ext ? gl.getParameter(ext.UNMASKED_RENDERER_WEBGL) : gl.getParameter(gl.RENDERER);
  } catch (e) { return 'unavailable'; }
}

// getContext() returns a flat { label: value } object describing the app state
export function initFeedback({ getContext, isBlocked = () => false }) {
  const $ = (id) => document.getElementById(id);
  const overlay = $('feedbackOverlay'), card = overlay.querySelector('.overlay-card');
  const btn = $('feedbackBtn'), form = $('feedbackForm');
  const summary = $('fbSummary'), details = $('fbDetails'), email = $('fbEmail');
  const diagOn = $('fbDiag'), diagText = $('fbDiagText'), send = $('fbSend'), status = $('fbStatus');
  let type = 'bug', opener = null, busy = false;
  const renderer = webglRenderer();

  const diagnostics = () => {
    const ctx = { ...getContext(),
      URL: location.origin + location.pathname,
      Browser: navigator.userAgent,
      Screen: `${screen.width}×${screen.height} @${devicePixelRatio}x, window ${innerWidth}×${innerHeight}`,
      WebGL: renderer,
      Language: navigator.language,
    };
    const lines = Object.entries(ctx).map(([k, v]) => `${k}: ${v}`);
    lines.push(`Recent errors: ${recentErrors.length ? '\n  ' + recentErrors.join('\n  ') : 'none'}`);
    return lines.join('\n');
  };
  const setStatus = (t, err = false) => { status.textContent = t; status.classList.toggle('error', err); };
  const setType = (t) => {
    type = t;
    overlay.querySelectorAll('[data-ftab]').forEach((b) => { const on = b.dataset.ftab === t; b.classList.toggle('active', on); b.setAttribute('aria-selected', on ? 'true' : 'false'); });
    $('fbDetailsLabel').textContent = DETAILS[t].label;
    details.placeholder = DETAILS[t].placeholder;
    diagOn.checked = t === 'bug';
  };
  const open = () => {
    opener = document.activeElement;
    diagText.textContent = diagnostics();
    setStatus('');
    overlay.classList.remove('hidden');
    summary.focus();
  };
  const close = () => {
    overlay.classList.add('hidden');
    const back = opener || btn;
    if (back && typeof back.focus === 'function') back.focus();
    opener = null;
  };
  const plainText = (p) => [
    `Type: ${p.type === 'bug' ? 'Bug report' : 'Feature request'}`,
    `Summary: ${p.summary}`,
    p.email ? `Reply to: ${p.email}` : null,
    '', p.details || '(no details)',
    p.diagnostics ? `\n--- Diagnostics ---\n${p.diagnostics}` : null,
  ].filter((l) => l != null).join('\n');

  btn.addEventListener('click', open);
  $('feedbackClose').addEventListener('click', close);
  overlay.addEventListener('click', (e) => { if (e.target === overlay) close(); });
  // on window, so Escape works wherever focus is while the dialog is open
  window.addEventListener('keydown', (e) => {
    if (overlay.classList.contains('hidden')) return;
    if (e.key === 'Escape') { e.preventDefault(); close(); return; }
    if (e.key === 'Tab') {
      const f = [...card.querySelectorAll('button, [href], input:not(.fb-hp), select, textarea, summary')].filter((el) => el.offsetParent !== null && !el.disabled);
      if (!f.length) return;
      const first = f[0], last = f[f.length - 1];
      if (e.shiftKey && document.activeElement === first) { e.preventDefault(); last.focus(); }
      else if (!e.shiftKey && document.activeElement === last) { e.preventDefault(); first.focus(); }
    }
  });
  overlay.querySelectorAll('[data-ftab]').forEach((b) => b.addEventListener('click', () => setType(b.dataset.ftab)));
  diagOn.addEventListener('change', () => { diagText.textContent = diagOn.checked ? diagnostics() : '(not included)'; });
  summary.addEventListener('input', () => summary.removeAttribute('aria-invalid'));
  setType('bug');

  form.addEventListener('submit', async (e) => {
    e.preventDefault();
    if (busy || isBlocked()) return;
    const p = {
      type,
      summary: summary.value.trim(),
      details: details.value.trim(),
      email: email.value.trim(),
      diagnostics: diagOn.checked ? diagnostics() : '',
      website: $('fbWebsite').value, // honeypot
    };
    if (!p.summary) { summary.setAttribute('aria-invalid', 'true'); summary.focus(); setStatus('Please add a short summary.', true); return; }
    if (p.email && !email.checkValidity()) { email.focus(); setStatus('That email address does not look right.', true); return; }
    busy = true; send.setAttribute('aria-disabled', 'true'); setStatus('Sending…');
    try {
      const res = await fetch('/api/feedback', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(p) });
      const data = await res.json().catch(() => ({}));
      if (res.ok && data.ok) {
        setStatus('Sent — thank you!');
        summary.value = ''; details.value = '';
        setTimeout(() => { if (!overlay.classList.contains('hidden')) close(); }, 1600);
      } else if (data.fallback === 'mailto' && data.to) {
        const subject = `[EchoSim ${p.type === 'bug' ? 'bug' : 'feature'}] ${p.summary}`;
        location.href = `mailto:${data.to}?subject=${encodeURIComponent(subject)}&body=${encodeURIComponent(plainText(p).slice(0, 1800))}`;
        setStatus('Opening your email app to send it…');
      } else if (res.status === 429) {
        setStatus('Too many messages — please try again in a few minutes.', true);
      } else if (res.status === 404 || res.status === 405 || (res.status === 503 && !data.to)) {
        // no backend here, or it has no recipient configured
        throw new Error('no endpoint');
      } else {
        setStatus(data.error || 'Could not send. Please try again.', true);
      }
    } catch (err) {
      // no backend (local static server) or offline: keep the text
      try { await navigator.clipboard.writeText(plainText(p)); setStatus('Could not reach the server — your message is copied to the clipboard.', true); }
      catch (e2) { setStatus('Could not reach the server. Please try again later.', true); }
    } finally {
      busy = false; send.removeAttribute('aria-disabled');
    }
  });

  return { isOpen: () => !overlay.classList.contains('hidden') };
}
