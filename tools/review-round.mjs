export const meta = {
  name: 'echosim-board-round',
  description: 'One board-review + technical-staff implementation round for EchoSim',
  phases: [
    { title: 'Board Review' },
    { title: 'Synthesis' },
    { title: 'Implementation' },
    { title: 'Verify' },
  ],
};

// args = { iter: number, shotsDir: string, targets: number }
const iter = args?.iter ?? 1;
const shotsDir = args?.shotsDir ?? '.review/iter0';
const TARGET = args?.targets ?? 9.5;

const CODE = [
  'index.html', 'css/style.css',
  'js/cardiac-model.js', 'js/echo.js', 'js/heart3d.js', 'js/main.js',
];
const FULLSHOTS = [
  'normal_plax_dia', 'normal_psax_sys', 'normal_a4c_dia', 'mr_a4c_sys',
  'as_plax_sys', 'dcm_a4c_sys', 'effusion_plax', 'vsd_a4c_sys', 'asd_subcostal',
].map((n) => `${shotsDir}/${n}.png`);
const ECHOSHOTS = ['normal_plax_dia', 'mr_a4c_sys', 'as_plax_sys', 'normal_psax_sys']
  .map((n) => `${shotsDir}/${n}_echo.png`);

const BOARD = [
  { metric: 'clinical', persona: 'a Consultant Echocardiographer (25 yrs)',
    title: 'Clinical & Anatomical Accuracy',
    focus: 'Are chambers, valves, great vessels, wall thicknesses and their motion anatomically credible? Do the pathologies show the correct, recognisable findings (e.g. MR = systolic jet into LA; AS = restricted valve + LVH; effusion = echo-free rind)? Do standard views resemble their real counterparts?' },
  { metric: 'physics', persona: 'a Diagnostic Ultrasound Physicist',
    title: 'Ultrasound & Doppler Physics Fidelity',
    focus: 'B-mode sector geometry, speckle, depth attenuation/TGC, log compression. Colour Doppler: BART convention, angle dependence (cos theta), aliasing past Nyquist, turbulence/variance mosaic. Spectral trace correctness. Physical plausibility of velocities.' },
  { metric: 'pedagogy', persona: 'a Medical Educator who teaches echo',
    title: 'Educational & Pedagogical Value',
    focus: 'Does linking the 3D plane to the 2D image build spatial understanding? Are teaching notes accurate and useful? Is there enough guidance/labelling to learn from? What would make it a better trainer (annotations, quiz, measurements, comparisons)?' },
  { metric: 'visual', persona: 'a Medical Visualization Artist',
    title: 'Visual Realism & Design',
    focus: 'Does the echo image read as a genuine ultrasound? Is the 3D heart attractive and legible? Colour palette, contrast, layout polish, light/dark balance, absence of visual glitches or floating/overlapping artefacts.' },
  { metric: 'ux', persona: 'a Senior UX Designer',
    title: 'Interaction & UX',
    focus: 'Discoverability and clarity of controls, feedback, affordances, responsiveness of layout, information hierarchy, onboarding. Are the probe/plane controls intuitive? Any confusing or dead UI?' },
  { metric: 'robustness', persona: 'a QA & Performance Engineer',
    title: 'Technical Robustness & Performance',
    focus: 'Runtime errors (see errors.json in the shots dir), frame stability, resize handling, edge cases (extreme sliders), memory/alloc in hot loops, offline self-containment. Does anything break or degrade?' },
  { metric: 'code', persona: 'a Software Architect',
    title: 'Code Quality & Maintainability',
    focus: 'Module boundaries, readability, naming, duplication, magic numbers, comments that explain intent, correctness of the maths, and how easily a new pathology/view could be added.' },
];

const SCORE_SCHEMA = {
  type: 'object',
  required: ['metric', 'score', 'justification', 'findings'],
  properties: {
    metric: { type: 'string' },
    score: { type: 'number', description: 'one decimal, 0-10' },
    justification: { type: 'string' },
    findings: {
      type: 'array',
      items: {
        type: 'object',
        required: ['severity', 'owner', 'problem', 'fix'],
        properties: {
          severity: { type: 'string', enum: ['high', 'medium', 'low'] },
          owner: { type: 'string', enum: ['model', 'echo', 'three', 'ui'], description: 'which subsystem file-owner should fix it' },
          problem: { type: 'string' },
          fix: { type: 'string', description: 'concrete, implementable change' },
        },
      },
    },
  },
};

phase('Board Review');
log(`Iteration ${iter}: convening ${BOARD.length}-member board (target ${TARGET}/10 all metrics)`);

const rubric = `Score to ONE decimal place, 0-10, calibrated for a HIGH bar:
10.0 = flawless, nothing to improve; ${TARGET} = production-grade for a serious teaching tool, only trivial cosmetic nitpicks remain; 9.0 = strong but with at least one clear, worth-fixing gap; 8.0 = good with several gaps; <=7 = notable deficiencies.
Be a demanding but FAIR domain expert. Reward genuine quality; do not inflate. If you score below ${TARGET}, you MUST list specific, implementable findings (each assigned to an owner: model=js/cardiac-model.js, echo=js/echo.js, three=js/heart3d.js, ui=index.html/css/style.css/js/main.js) that, if fixed, would raise the score to >=${TARGET}. Keep findings concrete and buildable in one round.`;

const boardResults = await parallel(BOARD.map((b) => () =>
  agent(
    `You are ${b.persona} sitting on the review board for "EchoSim", a browser-based 3D echocardiography training simulator. This is iteration ${iter}.

YOUR METRIC: "${b.title}".
What to judge: ${b.focus}

Read the source files: ${CODE.join(', ')}.
Look at the rendered screenshots (full app): ${FULLSHOTS.join(', ')}.
Close-up echo crops: ${ECHOSHOTS.join(', ')}.
Also read ${shotsDir}/errors.json (runtime errors captured this round).

${rubric}

Return ONLY the structured score object for your metric.`,
    { label: `board:${b.metric}`, phase: 'Board Review', schema: SCORE_SCHEMA, effort: 'medium' }
  ).then((r) => ({ ...r, metric: b.metric }))
));

const scored = boardResults.filter(Boolean);
const scoreLine = scored.map((s) => `${s.metric}=${s.score}`).join('  ');
log(`Scores: ${scoreLine}`);

phase('Synthesis');
const SYNTH_SCHEMA = {
  type: 'object',
  required: ['overallReady', 'perMetric', 'tasks'],
  properties: {
    overallReady: { type: 'boolean', description: `true only if every metric >= ${TARGET}` },
    summary: { type: 'string' },
    perMetric: {
      type: 'array',
      items: { type: 'object', required: ['metric', 'score'], properties: { metric: { type: 'string' }, score: { type: 'number' } } },
    },
    tasks: {
      type: 'array',
      description: 'deduplicated, prioritised, implementable tasks for this round',
      items: {
        type: 'object',
        required: ['owner', 'title', 'detail', 'priority'],
        properties: {
          owner: { type: 'string', enum: ['model', 'echo', 'three', 'ui'] },
          title: { type: 'string' },
          detail: { type: 'string', description: 'precise change, incl. any cross-file contract (e.g. new model field echo must read)' },
          priority: { type: 'number', description: '1=highest' },
        },
      },
    },
  },
};

const chair = await agent(
  `You are the Chair of the EchoSim review board. Here are this round's expert scores and findings as JSON:

${JSON.stringify(scored, null, 2)}

Synthesise into a single implementation plan for the technical staff. Deduplicate overlapping findings, resolve conflicts, and PRIORITISE the changes most likely to lift every metric to >= ${TARGET}. Assign each task to exactly one owner (model/echo/three/ui) that owns the file(s) it touches; if a change spans subsystems, put the primary edit under one owner and describe the contract other owners must honour in 'detail'. Keep the list focused (roughly 4-9 tasks) and each task concretely buildable in one round. Set overallReady = true ONLY if every metric is already >= ${TARGET}.`,
  { label: 'chair:synthesis', phase: 'Synthesis', schema: SYNTH_SCHEMA, effort: 'high' }
);

if (chair.overallReady) {
  log(`Chair: all metrics >= ${TARGET}. No implementation needed.`);
  return { iter, scored, chair, implemented: [] };
}

phase('Implementation');
const owners = ['model', 'echo', 'three', 'ui'];
const ownerFiles = {
  model: 'js/cardiac-model.js',
  echo: 'js/echo.js',
  three: 'js/heart3d.js',
  ui: 'index.html, css/style.css, js/main.js',
};
const tasksByOwner = Object.fromEntries(owners.map((o) => [o, (chair.tasks || []).filter((t) => t.owner === o)]));

function implPrompt(owner) {
  const mine = tasksByOwner[owner];
  return `You are the ${owner.toUpperCase()} engineer on the EchoSim technical staff (iteration ${iter}).
You may ONLY edit these file(s): ${ownerFiles[owner]}. Do NOT touch any other file (other engineers own them in parallel).

Apply these prioritised tasks from the board Chair:
${JSON.stringify(mine, null, 2)}

Full board plan for context (honour cross-file contracts noted in 'detail', but only edit YOUR files):
${JSON.stringify(chair.tasks, null, 2)}

Rules:
- Read your file(s) first, then make surgical, correct edits. Preserve existing working behaviour and public function/field names other modules rely on (classify, velocityAt, geometryAt, TISSUE, FLOW, EchoView, Heart3D and their fields) unless the plan explicitly changes a contract.
- Keep the app a self-contained ES-module static site. No new network deps.
- Match the surrounding code style and comment density.
- Do NOT run the app, npm, or playwright. Just edit files.
- If a task is infeasible or risky, skip it and say why.

Return a short bullet list of exactly what you changed (file + change), and note any contract you exposed for other owners.`;
}

// model first (others depend on its data contract), then echo/three/ui in parallel
let modelReport = null;
if (tasksByOwner.model.length) {
  modelReport = await agent(implPrompt('model'), { label: 'staff:model', phase: 'Implementation', effort: 'high' });
}
const parallelOwners = ['echo', 'three', 'ui'].filter((o) => tasksByOwner[o].length);
const reports = await parallel(parallelOwners.map((o) => () =>
  agent(implPrompt(o), { label: `staff:${o}`, phase: 'Implementation', effort: 'high' }).then((r) => ({ owner: o, report: r }))
));

phase('Verify');
const verify = await agent(
  `Verify the EchoSim JS still parses after this round's edits. Run:
  node --check js/cardiac-model.js && node --check js/echo.js && node --check js/heart3d.js && node --check js/main.js
Report PASS/FAIL and paste any syntax error. If a file FAILS to parse, open it and fix ONLY the syntax error you introduced (do not change behaviour). Then re-run the checks and confirm PASS.`,
  { label: 'verify:syntax', phase: 'Verify', effort: 'medium' }
);

return {
  iter,
  scored,
  perMetric: chair.perMetric,
  tasks: chair.tasks,
  implemented: [modelReport && { owner: 'model', report: modelReport }, ...reports].filter(Boolean),
  verify,
};
