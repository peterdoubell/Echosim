// Learning pathway dialog: modules (learn / practise / test), case questions
// played in the simulator, a timed mock exam, and progress saved per browser.
// Content and scoring live in js/pathway.js; main.js passes in a small bridge
// to the simulator ({ load, end, openReport, caseActive }).

import { MODULES, MODULE_BY_ID, PASS_MARK, EXAM, DISCLAIMER, moduleQuiz, mockExam, score, weakModules } from './pathway.js';

const KEY = 'echosim_pathway_v1';
const esc = (t) => String(t).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' })[c]);
const pct = (x) => `${Math.round(x * 100)}%`;

function loadStore() {
  try {
    const s = JSON.parse(localStorage.getItem(KEY) || 'null');
    if (s && s.modules) return { exams: [], ...s };
  } catch (e) {}
  return { modules: {}, exams: [], last: null };
}

export function initPathway(sim) {
  const $ = (id) => document.getElementById(id);
  const overlay = $('pathwayOverlay'), card = overlay.querySelector('.overlay-card');
  const nav = $('pwNav'), main = $('pwMain'), bar = $('pwCaseBar');
  let store = loadStore();
  let view = { kind: 'overview' };
  let quiz = null;      // { mode, moduleId, items, answers, revealed, idx, deadline, timer }
  let opener = null;

  const save = () => { try { localStorage.setItem(KEY, JSON.stringify(store)); } catch (e) {} };
  const mod = (id) => store.modules[id] || (store.modules[id] = { best: 0, attempts: 0, passed: false, tasks: {} });
  const isOpen = () => !overlay.classList.contains('hidden');

  // ---------- navigation column + overall progress ----------
  function renderNav() {
    const passed = MODULES.filter((m) => store.modules[m.id] && store.modules[m.id].passed).length;
    const bestExam = store.exams.reduce((b, e) => Math.max(b, e.pct), 0);
    $('pwBar').style.width = `${(passed / MODULES.length) * 100}%`;
    $('pwProgressText').textContent = `${passed} of ${MODULES.length} modules passed` + (store.exams.length ? ` · best mock exam ${pct(bestExam)}` : '');
    const item = (key, label, icon, status, active) =>
      `<button class="pw-nav-item${active ? ' active' : ''}" data-nav="${key}" ${active ? 'aria-current="page"' : ''}>
        <span class="pw-nav-icon" aria-hidden="true">${icon}</span><span class="pw-nav-label">${esc(label)}</span>${status}</button>`;
    const cur = view.kind === 'module' ? view.id : view.kind === 'exam-intro' ? 'exam' : view.kind === 'quiz' || view.kind === 'result' ? (quiz && quiz.mode === 'module' ? quiz.moduleId : 'exam') : view.kind;
    const status = (id) => {
      const s = store.modules[id];
      if (!s || !s.attempts) return '<span class="pw-st" aria-label="not started"></span>';
      return s.passed ? `<span class="pw-st ok" aria-label="passed, best ${pct(s.best)}">✓</span>` : `<span class="pw-st part" aria-label="best ${pct(s.best)}">${pct(s.best)}</span>`;
    };
    nav.innerHTML = item('overview', 'Overview', '🎓', '', cur === 'overview') +
      MODULES.map((m, i) => item(m.id, `${i + 1}. ${m.title}`, m.icon, status(m.id), cur === m.id)).join('') +
      item('exam', 'Mock exam', '⏱️', store.exams.length ? `<span class="pw-st${store.exams.some((e) => e.passed) ? ' ok' : ' part'}">${pct(bestExam)}</span>` : '', cur === 'exam');
    nav.querySelectorAll('[data-nav]').forEach((b) => b.addEventListener('click', () => go(b.dataset.nav)));
  }

  function go(key) {
    if (quiz && !leaveQuiz()) return;
    if (key === 'overview') view = { kind: 'overview' };
    else if (key === 'exam') view = { kind: 'exam-intro' };
    else { view = { kind: 'module', id: key, tab: 'learn' }; store.last = key; save(); }
    render();
  }
  // abandoning an attempt part-way loses it; ask first
  function leaveQuiz() {
    if (!quiz) return true;
    const done = quiz.finished;
    if (!done && !window.confirm('Leave this test? Your answers in this attempt will not be saved.')) return false;
    stopQuiz();
    return true;
  }
  function stopQuiz() {
    if (quiz && quiz.timer) clearInterval(quiz.timer);
    quiz = null;
    sim.end(true);
    hideBar();
  }

  // ---------- pages ----------
  function render() {
    renderNav();
    if (view.kind === 'overview') renderOverview();
    else if (view.kind === 'module') renderModule();
    else if (view.kind === 'exam-intro') renderExamIntro();
    else if (view.kind === 'quiz') renderQuestion();
    else if (view.kind === 'result') renderResult();
    main.scrollTop = 0;
  }

  function renderOverview() {
    const next = MODULES.find((m) => !(store.modules[m.id] && store.modules[m.id].passed));
    const resume = store.last && MODULE_BY_ID[store.last] ? MODULE_BY_ID[store.last] : next;
    const history = store.exams.slice(-5).reverse().map((e) =>
      `<li>${new Date(e.at).toLocaleDateString()} — <b>${e.correct}/${e.total}</b> (${pct(e.pct)}) ${e.passed ? '<span class="sev-badge sev-normal">pass</span>' : '<span class="sev-badge sev-moderate">below pass mark</span>'}</li>`).join('');
    main.innerHTML = `
      <h3>Prepare for the BSE TTE accreditation exam</h3>
      <p>Work through ten modules that follow the main topics of the transthoracic written exam. Each module has three parts:</p>
      <ul class="pw-list">
        <li><b>Learn</b> — objectives and the key facts and numbers, using the BSE 2020 reference intervals.</li>
        <li><b>Practise</b> — tasks that load the matching case into the simulator.</li>
        <li><b>Test</b> — questions with explanations. Pass at ${pct(PASS_MARK)}. Case questions play a hidden case in the simulator, like the exam's video questions.</li>
      </ul>
      <p>Then sit the <b>mock exam</b>: ${EXAM.questions} questions across all modules in ${EXAM.minutes} minutes, with a breakdown by topic.</p>
      <p class="pw-note">Accreditation also needs a supervised case logbook and a practical assessment, which a simulator cannot replace. ${esc(DISCLAIMER)}</p>
      <div class="pw-actions">
        ${resume ? `<button class="pw-primary" data-go="${resume.id}">${store.last ? 'Continue' : 'Start'}: ${esc(resume.title)} →</button>` : '<span class="pw-done">All modules passed — try the mock exam.</span>'}
        <button data-go="exam">⏱️ Mock exam</button>
      </div>
      ${history ? `<h3>Recent mock exams</h3><ul class="pw-history">${history}</ul>` : ''}
      <p class="pw-small">Progress is saved in this browser only. <button class="pw-link" id="pwReset">Reset progress</button></p>`;
    main.querySelectorAll('[data-go]').forEach((b) => b.addEventListener('click', () => go(b.dataset.go)));
    $('pwReset').addEventListener('click', () => {
      if (!window.confirm('Reset all pathway progress and exam history in this browser?')) return;
      store = { modules: {}, exams: [], last: null }; save(); render();
    });
  }

  function renderModule() {
    const m = MODULE_BY_ID[view.id], s = mod(view.id);
    const tabs = [['learn', 'Learn'], ['practise', 'Practise'], ['test', 'Test']].map(([k, l]) =>
      `<button class="seg-btn${view.tab === k ? ' active' : ''}" role="tab" aria-selected="${view.tab === k}" data-mtab="${k}">${l}</button>`).join('');
    let body = '';
    if (view.tab === 'learn') {
      body = `<h4>Objectives</h4><ul class="pw-list">${m.objectives.map((o) => `<li>${esc(o)}</li>`).join('')}</ul>
        <h4>Key points</h4><ul class="pw-list pw-key">${m.learn.map((l) => `<li>${l}</li>`).join('')}</ul>
        <div class="pw-actions"><button class="pw-primary" data-mtab="practise">Practise →</button></div>`;
    } else if (view.tab === 'practise') {
      body = `<p>Each task opens the case in the simulator. Close the pathway to work on it, and tick it off when you are done.</p>
        <ul class="pw-tasks">${m.tasks.map((t) => `<li><label class="pw-task"><input type="checkbox" data-task="${t.id}" ${s.tasks[t.id] ? 'checked' : ''}/> <span>${esc(t.text)}</span></label>
          <button data-open="${t.id}">${t.action === 'report' ? '📝 Open case + report' : '▶ Open in simulator'}</button></li>`).join('')}</ul>
        <div class="pw-actions"><button class="pw-primary" data-mtab="test">Take the test →</button></div>`;
    } else {
      const n = m.questions.length, cases = m.questions.filter((q) => q.case).length;
      body = `<p>${n} questions${cases ? `, including ${cases} case question${cases > 1 ? 's' : ''} played in the simulator` : ''}. You see the answer and explanation after each one. Pass mark ${pct(PASS_MARK)}.</p>
        <p>${s.attempts ? `Attempts: ${s.attempts} · best ${pct(s.best)} ${s.passed ? '<span class="sev-badge sev-normal">passed</span>' : ''}` : 'Not attempted yet.'}</p>
        <div class="pw-actions"><button class="pw-primary" id="pwStartModule">Start test</button></div>`;
    }
    main.innerHTML = `<h3>${m.icon} ${esc(m.title)}</h3><div class="report-tabs seg" role="tablist" aria-label="Module sections">${tabs}</div><div class="pw-body">${body}</div>`;
    main.querySelectorAll('[data-mtab]').forEach((b) => b.addEventListener('click', () => { view.tab = b.dataset.mtab; render(); }));
    main.querySelectorAll('[data-task]').forEach((c) => c.addEventListener('change', () => { s.tasks[c.dataset.task] = c.checked; save(); }));
    main.querySelectorAll('[data-open]').forEach((b) => b.addEventListener('click', () => {
      const t = m.tasks.find((x) => x.id === b.dataset.open);
      sim.load({ ...t.sim, hidden: false });
      hide();
      if (t.action === 'report') setTimeout(() => sim.openReport(), 50);
    }));
    const start = $('pwStartModule');
    if (start) start.addEventListener('click', () => startQuiz('module', view.id));
  }

  function renderExamIntro() {
    const last = store.exams[store.exams.length - 1];
    main.innerHTML = `<h3>⏱️ Mock exam</h3>
      <p>${EXAM.questions} questions drawn evenly from all ten modules, including case questions played in the simulator. You have <b>${EXAM.minutes} minutes</b>; the exam submits itself when time runs out.</p>
      <p>Answers are not shown until you finish. You can move back and forth and change answers. Pass mark ${pct(PASS_MARK)}.</p>
      <p class="pw-note">The real exam's format, length and pass mark are set by the BSE and may differ; this is practice.</p>
      ${last ? `<p>Last attempt: ${last.correct}/${last.total} (${pct(last.pct)}).</p>` : ''}
      <div class="pw-actions"><button class="pw-primary" id="pwStartExam">Start mock exam</button></div>`;
    $('pwStartExam').addEventListener('click', () => startQuiz('exam'));
  }

  // ---------- quiz runner (module tests and the mock exam) ----------
  function startQuiz(mode, moduleId) {
    const items = mode === 'exam' ? mockExam() : moduleQuiz(moduleId);
    quiz = { mode, moduleId, items, answers: items.map(() => null), revealed: items.map(() => false), idx: 0, finished: false };
    if (mode === 'exam') {
      quiz.deadline = Date.now() + EXAM.minutes * 60000;
      quiz.timer = setInterval(tick, 1000);
    }
    view = { kind: 'quiz' };
    render();
  }
  function tick() {
    if (!quiz || !quiz.deadline) return;
    const left = Math.max(0, quiz.deadline - Date.now());
    const t = `${Math.floor(left / 60000)}:${String(Math.floor(left / 1000) % 60).padStart(2, '0')}`;
    document.querySelectorAll('.pw-timer').forEach((el) => { el.textContent = `⏱ ${t}`; el.classList.toggle('low', left < 5 * 60000); });
    if (left <= 0) finishQuiz();
  }
  // put the current question's case (if any) into the simulator
  function syncCase() {
    const q = quiz.items[quiz.idx];
    const want = `${quiz.idx}:${quiz.revealed[quiz.idx]}`;
    if (quiz.loaded === want) return;     // re-renders keep the learner's probe adjustments
    quiz.loaded = want;
    if (q.case && !(quiz.mode === 'module' && quiz.revealed[quiz.idx])) sim.load({ ...q.case, hidden: true });
    else if (q.case) sim.end(false);            // answered: show the case unlocked
    else sim.end(true);
  }
  function renderQuestion() {
    const q = quiz.items[quiz.idx], n = quiz.items.length, i = quiz.idx;
    const exam = quiz.mode === 'exam', revealed = quiz.revealed[i], chosen = quiz.answers[i];
    syncCase();
    const title = exam ? 'Mock exam' : `Test: ${MODULE_BY_ID[quiz.moduleId].title}`;
    const dots = quiz.items.map((_, k) => {
      const cls = k === i ? 'cur' : exam ? (quiz.answers[k] != null ? 'done' : '') : quiz.revealed[k] ? (quiz.answers[k] === quiz.items[k].correct ? 'ok' : 'no') : '';
      return `<button class="pw-dot ${cls}" data-jump="${k}" aria-label="Question ${k + 1}" ${exam ? '' : 'tabindex="-1"'}></button>`;
    }).join('');
    const opts = q.shown.map((o, k) => {
      let cls = '';
      if (revealed) cls = k === q.correct ? 'correct' : k === chosen ? 'wrong' : '';
      else if (k === chosen) cls = 'sel';
      return `<button class="pw-opt ${cls}" data-opt="${k}" ${revealed ? 'disabled' : ''} aria-pressed="${k === chosen}"><span class="pw-opt-l">${'ABCDE'[k]}</span>${esc(o)}</button>`;
    }).join('');
    const caseBox = q.case ? `<div class="pw-case">📺 <b>Case question.</b> The simulator is playing a case${q.case.hideView ? ' with the view name hidden' : ' with the diagnosis hidden'}.
      <button id="pwWatch">Watch the case</button></div>` : '';
    const feedback = revealed ? `<div class="pw-fb ${chosen === q.correct ? 'ok' : 'no'}" role="status"><b>${chosen === q.correct ? '✓ Correct' : `✗ The answer is ${'ABCDE'[q.correct]}: ${esc(q.shown[q.correct])}`}</b><br>${esc(q.explain)}</div>` : '';
    const last = i === n - 1;
    let actions;
    if (exam) {
      actions = `<button id="pwPrev" ${i === 0 ? 'disabled' : ''}>← Previous</button>
        ${last ? '' : '<button id="pwNext" class="pw-primary">Next →</button>'}
        <button id="pwFinish" class="${last ? 'pw-primary' : ''}">Finish exam</button>`;
    } else {
      actions = revealed ? `<button id="pwNext" class="pw-primary">${last ? 'See results' : 'Next question →'}</button>` : '<span class="pw-small">Choose an answer.</span>';
    }
    main.innerHTML = `<div class="pw-qhead"><b>${esc(title)}</b><span>Question ${i + 1} of ${n}</span>${exam ? '<span class="pw-timer"></span>' : ''}</div>
      <div class="pw-dots">${dots}</div>
      <p class="pw-q">${esc(q.q)}</p>${caseBox}
      <div class="pw-opts" role="group" aria-label="Answers">${opts}</div>${feedback}
      <div class="pw-actions">${actions}</div>`;
    if (exam) tick();
    main.querySelectorAll('[data-opt]').forEach((b) => b.addEventListener('click', () => choose(+b.dataset.opt)));
    main.querySelectorAll('[data-jump]').forEach((b) => b.addEventListener('click', () => { if (exam) { quiz.idx = +b.dataset.jump; render(); } }));
    const on = (id, f) => { const el = $(id); if (el) el.addEventListener('click', f); };
    on('pwWatch', showBar);
    on('pwPrev', () => { quiz.idx--; render(); });
    on('pwNext', () => { if (!exam && last) finishQuiz(); else { quiz.idx++; render(); } });
    on('pwFinish', () => {
      const blank = quiz.answers.filter((a) => a == null).length;
      if (blank && !window.confirm(`${blank} question${blank > 1 ? 's are' : ' is'} unanswered. Finish anyway?`)) return;
      finishQuiz();
    });
    // keep keyboard focus inside the question area
    const focus = main.querySelector(revealed ? '#pwNext' : '.pw-opt.sel, .pw-opt');
    if (focus && isOpen()) focus.focus({ preventScroll: true });
  }
  function choose(k) {
    const i = quiz.idx;
    if (quiz.revealed[i]) return;
    quiz.answers[i] = k;
    if (quiz.mode === 'module') quiz.revealed[i] = true;
    if (!bar.classList.contains('hidden')) { hideBar(); show(); }
    render();
  }
  function finishQuiz() {
    if (!quiz || quiz.finished) return;
    if (quiz.timer) clearInterval(quiz.timer);
    quiz.finished = true;
    const r = score(quiz.items, quiz.answers);
    quiz.result = r;
    if (quiz.mode === 'module') {
      const s = mod(quiz.moduleId);
      s.attempts++; s.best = Math.max(s.best, r.pct); s.passed = s.passed || r.passed;
    } else {
      store.exams.push({ at: Date.now(), correct: r.correct, total: r.total, pct: r.pct, passed: r.passed, byModule: r.byModule });
      store.exams = store.exams.slice(-20);
    }
    save();
    sim.end(true);
    hideBar();
    if (!isOpen()) show();
    view = { kind: 'result' };
    render();
  }
  function renderResult() {
    const r = quiz.result, exam = quiz.mode === 'exam';
    const rows = Object.entries(r.byModule).map(([id, b]) => {
      const p = b.correct / b.total;
      return `<tr><th scope="row">${MODULE_BY_ID[id].icon} ${esc(MODULE_BY_ID[id].title)}</th><td>${b.correct}/${b.total}</td>
        <td><div class="pw-mini"><span class="${p >= PASS_MARK ? 'ok' : 'no'}" style="width:${p * 100}%"></span></div></td></tr>`;
    }).join('');
    const weak = weakModules(r);
    const review = quiz.items.map((q, i) => {
      const a = quiz.answers[i], ok = a === q.correct;
      return `<details class="pw-rev ${ok ? 'ok' : 'no'}"><summary>${ok ? '✓' : '✗'} ${i + 1}. ${esc(q.q)}${q.case ? ' <em>(case)</em>' : ''}</summary>
        <p>Your answer: ${a == null ? '<i>not answered</i>' : esc(q.shown[a])}<br>Correct: <b>${esc(q.shown[q.correct])}</b></p><p>${esc(q.explain)}</p></details>`;
    }).join('');
    main.innerHTML = `<h3>${exam ? '⏱️ Mock exam result' : `Test result: ${esc(MODULE_BY_ID[quiz.moduleId].title)}`}</h3>
      <div class="pw-score ${r.passed ? 'ok' : 'no'}"><b>${r.correct}/${r.total}</b> (${pct(r.pct)}) — ${r.passed ? 'pass' : `below the ${pct(PASS_MARK)} pass mark`}</div>
      ${exam ? `<table class="report-table pw-bymod"><caption>By module</caption><tbody>${rows}</tbody></table>` : ''}
      ${weak.length ? `<p>Revise: ${weak.map((id) => `<button class="pw-link" data-go="${id}">${esc(MODULE_BY_ID[id].title)}</button>`).join(' · ')}</p>` : ''}
      <div class="pw-actions">
        <button class="pw-primary" id="pwRetry">${exam ? 'New mock exam' : 'Retake test'}</button>
        <button data-go="${exam ? 'overview' : quiz.moduleId}">${exam ? 'Back to overview' : 'Back to module'}</button>
      </div>
      <h4>Review your answers</h4>${review}`;
    $('pwRetry').addEventListener('click', () => { const m = quiz.mode, id = quiz.moduleId; quiz = null; startQuiz(m, id); });
    main.querySelectorAll('[data-go]').forEach((b) => b.addEventListener('click', () => { quiz = null; go(b.dataset.go); }));
  }

  // ---------- the case bar: answer while watching the simulator ----------
  function showBar() {
    const q = quiz.items[quiz.idx];
    $('pwCbTitle').textContent = `${quiz.mode === 'exam' ? 'Mock exam' : 'Test'} · question ${quiz.idx + 1} of ${quiz.items.length}`;
    $('pwCbQ').textContent = q.q;
    $('pwCbOpts').innerHTML = q.shown.map((o, k) =>
      `<button class="pw-opt${quiz.answers[quiz.idx] === k ? ' sel' : ''}" data-cbopt="${k}"><span class="pw-opt-l">${'ABCDE'[k]}</span>${esc(o)}</button>`).join('');
    $('pwCbOpts').querySelectorAll('[data-cbopt]').forEach((b) => b.addEventListener('click', () => choose(+b.dataset.cbopt)));
    $('pwCbTimer').className = quiz.mode === 'exam' ? 'pw-timer' : '';
    overlay.classList.add('hidden');
    bar.classList.remove('hidden');
    if (quiz.mode === 'exam') tick();
    $('pwCbBack').focus();
  }
  function hideBar() { bar.classList.add('hidden'); }
  $('pwCbBack').addEventListener('click', () => { hideBar(); show(); render(); });

  // ---------- open / close ----------
  function show() {
    overlay.classList.remove('hidden');
  }
  function open() {
    opener = document.activeElement;
    if (!bar.classList.contains('hidden')) hideBar();
    show();
    render();
    const f = main.querySelector('button');
    (f || $('pathwayClose')).focus();
  }
  function hide() {
    overlay.classList.add('hidden');
    const back = opener || $('pathwayBtn');
    if (back && typeof back.focus === 'function') back.focus();
    opener = null;
  }
  // closing mid-attempt: a case question minimises to the bar; otherwise ask
  function close() {
    if (quiz && !quiz.finished) {
      if (quiz.items[quiz.idx].case && !quiz.revealed[quiz.idx]) { showBar(); return; }
      if (!leaveQuiz()) return;
    } else if (quiz && quiz.finished) { quiz = null; view = { kind: 'overview' }; }
    hide();
  }

  document.querySelectorAll('[data-open-pathway]').forEach((b) => b.addEventListener('click', open));
  $('pathwayClose').addEventListener('click', close);
  overlay.addEventListener('click', (e) => { if (e.target === overlay) close(); });
  window.addEventListener('keydown', (e) => {
    if (!isOpen()) return;
    if (e.key === 'Escape') { e.preventDefault(); close(); return; }
    if (e.key === 'Tab') {
      const f = [...card.querySelectorAll('button, [href], input, select, textarea, summary')].filter((el) => el.offsetParent !== null && !el.disabled && el.tabIndex !== -1);
      if (!f.length) return;
      const first = f[0], last = f[f.length - 1];
      if (e.shiftKey && document.activeElement === first) { e.preventDefault(); last.focus(); }
      else if (!e.shiftKey && document.activeElement === last) { e.preventDefault(); first.focus(); }
    }
  });

  return { isOpen, open };
}
