/* ── Overlays section ────────────────────── */
async function loadOverlays() {
  try {
    const r = await fetch('/api/overlays');
    allOverlays = (await r.json()) || [];
  } catch { allOverlays = []; }
  renderOverlays();
}

// renderOverlays renders the .sqf overlay table, keeping only entries whose
// name matches every whitespace-separated term of filter.
function renderOverlays(filter) {
  const terms = searchTerms(filter);
  const q     = terms.length;
  const sqfs  = allOverlays.filter(o => o.type === 'sqf');
  const rows  = q ? sqfs.filter(o => matchesAllTerms(o.name, terms)) : sqfs;
  const tbody = gid('ov-tbody');
  if (!rows.length) {
    tbody.innerHTML = '<tr><td colspan="4" class="td-empty">' +
      (q ? 'No matches.' : 'No module overlays found.') + '</td></tr>';
    return;
  }
  tbody.innerHTML = rows.map(o =>
    '<tr>' +
      '<td class="td-name">' + escHtml(o.name) + '</td>' +
      '<td class="mono td-muted td-fit">' + fmtSize(o.size) + '</td>' +
      '<td class="mono path-cell td-muted">' + pathTailHtml(o.path) + '</td>' +
      '<td class="td-nowrap td-fit">' +
        '<button class="btn btn-sm btn-flat btn-icon" onclick="copyStr(\'' + escHtml(o.path) + '\',this)" title="Copy path">' + iconSvg('content_copy') + '</button>' +
        ' <button class="btn btn-sm btn-flat btn-icon btn-icon-danger" onclick="deleteOverlay(\'' + escHtml(o.path) + '\',\'' + escHtml(o.name) + '\')" title="Delete overlay">' + iconSvg('delete') + '</button>' +
      '</td>' +
    '</tr>'
  ).join('');
}

// deleteOverlay confirms, then removes an installed overlay via POST /api/remove.
async function deleteOverlay(path, name) {
  const ok = await askConfirm('Delete overlay "' + name + '"? This cannot be undone.',
    { title: 'Delete overlay', okLabel: 'Delete' });
  if (!ok) return;
  let r;
  try {
    r = await fetch('/api/remove', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ path }),
    });
  } catch (e) { showProgressError('Delete ' + name, String(e)); return; }
  if (!r.ok) { showProgressError('Delete ' + name, await r.text()); return; }
  loadOverlays();
}

function _addModuleOverlay(path) {
  const entry = allOverlays.find(o => o.path === path);
  if (entry && !selectedModules.find(m => m.path === path)) {
    selectedModules.push(entry);
    renderModuleChips();
  }
}

/* ── Overlay Picker Modal ────────────────── */
let _opFileController = null;

function openOverlayPicker(type) {
  opTargetType = type;
  gid('op-app-search').value = '';
  filterAppOverlays('');
  opSwitchTab('module');
  gid('op-modal').classList.add('open');
}

function opSwitchTab(tab) {
  const isModule = tab === 'module';
  gid('op-tab-module-btn').classList.toggle('active', isModule);
  gid('op-tab-file-btn').classList.toggle('active', !isModule);
  gid('op-tab-app').classList.toggle('active', isModule);
  gid('op-tab-file').classList.toggle('active', !isModule);
  if (!isModule) opFileNavigate(gid('cfg-cwd')?.value || srvScratch || srvHome || '/');
}

function filterAppOverlays(q) {
  const sqfOvs = allOverlays.filter(o => o.type === 'sqf');
  const terms  = searchTerms(q);
  const hits   = terms.length ? sqfOvs.filter(o => matchesAllTerms(o.name, terms)) : sqfOvs;
  gid('op-app-list').innerHTML = hits.length
    ? hits.map(o =>
        '<div class="modal-row" onclick="selectOverlay(\'' + escHtml(o.path) + '\')">' +
          '<span class="modal-row-icon">' + iconSvg('hexagon') + '</span>' +
          '<span class="modal-row-name fw-bold">' + escHtml(o.name) + '</span>' +
          '<span class="modal-row-size mono">' + fmtSize(o.size) + '</span>' +
          '<span class="modal-row-sel"><button class="btn btn-sm btn-primary">Select</button></span>' +
        '</div>'
      ).join('')
    : '<div class="modal-empty">' +
        (q ? 'No matches.' : 'No .sqf overlays found.') + '</div>';
}

async function opFileNavigate(path) {
  if (_opFileController) _opFileController.abort();
  _opFileController = new AbortController();
  const url = path ? '/api/fs?path=' + encodeURIComponent(path) : '/api/fs';
  try {
    const r       = await fetch(url, { signal: _opFileController.signal });
    const entries = (await r.json()) || [];
    const curPath = entries.length > 0 ? entries[0].path.replace(/\/[^/]+$/, '') || '/' : (path || '/');
    _renderBreadcrumb('op-file-bc', curPath, 'opFileNavigate');
    const dirs  = entries.filter(e =>  e.is_dir && !e.name.startsWith('.'));
    const files = entries.filter(e => !e.is_dir && e.name.endsWith('.sqf'));
    gid('op-file-list').innerHTML =
      dirs.map(d =>
        '<div class="modal-row" onclick="opFileNavigate(\'' + escHtml(d.path) + '\')">' +
          '<span class="modal-row-icon">' + iconSvg('folder', 'icon-folder', true) + '</span>' +
          '<span class="modal-row-name">' + escHtml(d.name) + '</span>' +
          '<span class="modal-row-sel">' + iconSvg('chevron_right') + '</span>' +
        '</div>'
      ).join('') +
      (files.length
        ? files.map(f =>
            '<div class="modal-row" onclick="selectExternalOverlay(\'' + escHtml(f.path) + '\')">' +
              '<span class="modal-row-icon">' + iconSvg('draft') + '</span>' +
              '<span class="modal-row-name">' + escHtml(f.name) + '</span>' +
              '<span class="modal-row-size mono">' + fmtSize(f.size) + '</span>' +
              '<span class="modal-row-sel"><button class="btn btn-sm btn-primary" onclick="event.stopPropagation();selectExternalOverlay(\'' + escHtml(f.path) + '\')">Select</button></span>' +
            '</div>'
          ).join('')
        : (!dirs.length ? '<div class="modal-empty">No .sqf files here.</div>' : ''));
  } catch(err) {
    if (err.name === 'AbortError') return;
    gid('op-file-list').innerHTML = '<div class="modal-error">Error loading directory.</div>';
  }
}

function selectOverlay(path) {
  _addModuleOverlay(path);
  closeModal('op-modal');
}

function selectExternalOverlay(path) {
  _addExternalOverlay(path);
  closeModal('op-modal');
}

/* ── Overlay Create dialog ──────────────── */
// Browse /api/avail entries (build scripts, templates collapsed); when no
// script matches, look up an exact conda package, then install via POST /api/create.
let _ocAvail = null;       // /api/avail entries (null until loaded)
let _ocShown = [];         // entries currently rendered in the browse list
let _ocConda = [];         // exact conda match(es) for _ocCondaQuery
let _ocCondaQuery = '';
let _ocCondaTimer = null;
let _ocSel = null;         // {kind:'script', entry, vars} | {kind:'conda', conda, version}
let _ocInspectRes = null;  // last /api/avail/inspect response
let _ocInspectTimer = null;
let _ocInspectPending = false; // true while an inspect request is scheduled/in flight

async function openOverlayCreate() {
  gid('oc-search').value = '';
  _ocConda = []; _ocCondaQuery = '';
  ocShowBrowse();
  gid('oc-modal').classList.add('open');
  gid('oc-list').innerHTML = '<div class="modal-empty">Loading…</div>';
  try {
    const r = await fetch('/api/avail');
    _ocAvail = (await r.json()) || [];
  } catch { _ocAvail = []; }
  ocFilter(gid('oc-search').value);
}

function ocShowBrowse() {
  gid('oc-browse').classList.add('active');
  gid('oc-detail-wrap').classList.remove('active');
  gid('oc-back-btn').hidden = true;
  gid('oc-create-btn').hidden = true;
  _ocSel = null; _ocInspectRes = null;
  clearTimeout(_ocInspectTimer);
}

// ocFilter filters the avail list: entries must match every whitespace-
// separated term of q (across name/alias/whatis/target).
function ocFilter(q) {
  if (_ocAvail === null) return;
  const lower = (q || '').toLowerCase().trim();
  const terms = searchTerms(lower);
  let hits = terms.length
    ? _ocAvail.filter(e =>
        matchesAllTerms([e.name, e.alias, e.whatis, e.target_template].join(' '), terms))
    : _ocAvail;
  if (terms.length) {
    // Exact-first ordering (mirrors the CLI's exact-first search): full-name or
    // alias match, then path-segment match, then name prefix/substring. Ranked
    // by the first term when the query has several.
    const key = terms[0];
    const rank = e => {
      const n = e.name.toLowerCase();
      if (n === key || (e.alias || '').toLowerCase() === key) return 0;
      if (n.split('/').includes(key)) return 1;
      if (n.startsWith(key)) return 2;
      if (n.includes(key)) return 3;
      return 4; // matched via whatis/target only
    };
    hits = hits.slice().sort((a, b) => rank(a) - rank(b));
  }
  // Always check conda by exact name (debounced) so its result shows alongside
  // build scripts — unless a script already owns that exact name (shadowed
  // below). A multi-term query is never a package name, so it skips conda.
  clearTimeout(_ocCondaTimer);
  if (terms.length === 1 && lower.length >= 2 && lower !== _ocCondaQuery && !_ocNameOwnedByScript(lower)) {
    _ocCondaTimer = setTimeout(() => _ocCondaSearch(lower), 400);
  }
  _ocUpdateAnacondaLink(lower);
  _ocRenderList(hits, lower);
}

// _ocNameOwnedByScript reports whether an avail entry already covers name —
// as an exact name, a default-distro alias, or a "name/<version>" entry.
// Such a conda result is redundant and gets shadowed (e.g. cytoscape, igv).
function _ocNameOwnedByScript(name) {
  return (_ocAvail || []).some(e => {
    const n = e.name.toLowerCase();
    return n === name || n.startsWith(name + '/') || (e.alias || '').toLowerCase() === name;
  });
}

async function _ocCondaSearch(q) {
  try {
    const r = await fetch('/api/search?q=' + encodeURIComponent(q));
    const data = r.ok ? await r.json() : null;
    _ocConda = (data && data.results) || [];
  } catch { _ocConda = []; }
  _ocCondaQuery = q;
  const cur = gid('oc-search').value.toLowerCase().trim();
  if (cur === q) ocFilter(cur);
}

// _ocUpdateAnacondaLink points the footer link at anaconda.org's own search
// for the current query (fuzzy discovery, which this tool intentionally omits).
function _ocUpdateAnacondaLink(q) {
  const a = gid('oc-anaconda-link');
  if (a) a.href = 'https://anaconda.org/search' + (q ? '?q=' + encodeURIComponent(q) : '');
}

// _ocCondaRowsHtml renders the current conda results as clickable rows.
function _ocCondaRowsHtml() {
  return _ocConda.map((r, i) =>
    '<div class="modal-row" onclick="ocSelectConda(' + i + ')">' +
      '<span class="modal-row-icon">' + iconSvg('draft') + '</span>' +
      '<span class="modal-row-name"><span class="oc-row-name">' + escHtml(r.name) +
        ' <span class="td-muted">(' + escHtml(r.owner) + ')</span></span>' +
        (r.summary ? '<span class="oc-whatis">' + escHtml(r.summary) + '</span>' : '') + '</span>' +
    '</div>').join('');
}

function _ocBadges(e) {
  const badges = [];
  if (e.installed)   badges.push('<span class="badge badge-green">installed</span>');
  if (e.is_template) badges.push('<span class="badge badge-amber">template</span>');
  if (e.container)   badges.push('<span class="badge badge-grey">container</span>');
  if (e.remote)      badges.push('<span class="badge badge-grey">remote</span>');
  return badges.join(' ');
}

function _ocRenderList(hits, q) {
  _ocShown = hits;
  let html = hits.map((e, i) =>
    '<div class="modal-row" onclick="ocSelect(' + i + ')">' +
      '<span class="modal-row-icon">' + iconSvg('hexagon') + '</span>' +
      '<span class="modal-row-name"><span class="oc-row-name">' + escHtml(e.name) + '</span>' +
        (e.alias ? ' <span class="oc-alias">[' + escHtml(e.alias) + ']</span>' : '') +
        (e.whatis ? '<span class="oc-whatis">' + escHtml(e.whatis) + '</span>' : '') + '</span>' +
      '<span class="modal-row-size">' + _ocBadges(e) + '</span>' +
    '</div>'
  ).join('');

  // A script owning the exact name shadows the conda result (e.g. cytoscape).
  // Multi-term queries never hit conda (see ocFilter), so they get plain
  // "No matches." instead of the conda messaging.
  const condaable  = q.length >= 2 && !/\s/.test(q);
  const shadowed   = condaable && _ocNameOwnedByScript(q);
  const condaReady = condaable && !shadowed && q === _ocCondaQuery;
  const condaHdr   = '<div class="oc-group-hdr">Conda package (anaconda.org)</div>';
  const condaRows  = condaReady && _ocConda.length ? condaHdr + _ocCondaRowsHtml() : '';

  if (!hits.length) {
    if (condaRows) {
      html = condaRows;
    } else if (condaReady) {
      html = '<div class="modal-empty">No build script or conda package named “' + escHtml(q) + '”.</div>';
    } else if (condaable) {
      html = '<div class="modal-empty">No build scripts match. Checking conda…</div>';
    } else {
      html = '<div class="modal-empty">' + (q ? 'No matches.' : 'No build scripts found.') + '</div>';
    }
  } else if (condaRows) {
    html = '<div class="oc-group-hdr">Build scripts</div>' + html + condaRows;
  }
  gid('oc-list').innerHTML = html;
}

function _ocShowDetail() {
  gid('oc-browse').classList.remove('active');
  gid('oc-detail-wrap').classList.add('active');
  gid('oc-back-btn').hidden = false;
  gid('oc-create-btn').hidden = false;
}

function ocSelect(i) {
  const e = _ocShown[i];
  if (!e) return;
  // vars hold user-entered values; empty means "use the default (latest)".
  _ocSel = { kind: 'script', entry: e, vars: {}, defaults: {} };
  if (e.is_template) {
    for (const k of (e.pl_order || [])) {
      const vals = (e.pl || {})[k] || [];
      _ocSel.defaults[k] = vals.find(v => v !== '*') || '';
      _ocSel.vars[k] = '';
    }
  }
  _ocRenderDetail();
  _ocShowDetail();
  _ocScheduleInspect();
}

function ocSelectConda(i) {
  const c = _ocConda[i];
  if (!c) return;
  // version empty = use latest.
  _ocSel = { kind: 'conda', conda: c, version: '', latest: (c.versions && c.versions[0]) || '' };
  _ocRenderDetail();
  _ocShowDetail();
  ocValidate();
}

// _ocConcreteName returns the name/version spec that will be passed to `create`.
// Empty placeholder/version inputs fall back to their default (latest) value.
function _ocConcreteName() {
  if (!_ocSel) return '';
  if (_ocSel.kind === 'conda') {
    const c = _ocSel.conda;
    const v = _ocSel.version || _ocSel.latest;
    // CondaSearchResult serializes the channel as "owner" (anaconda.org naming).
    return c.owner + '::' + c.name + (v ? '/' + v : '');
  }
  const e = _ocSel.entry;
  if (!e.is_template) return e.name;
  return (e.target_template || '').replace(/\{([^}]+)\}/g,
    (m, k) => (_ocSel.vars[k] || _ocSel.defaults[k] || m));
}

function _ocRenderDetail() {
  let html = '';
  if (_ocSel.kind === 'conda') {
    const c = _ocSel.conda;
    html += '<div class="oc-hdr"><div class="oc-title">' + escHtml(c.name) +
      ' <span class="td-muted">(' + escHtml(c.owner) + ')</span>' +
      ' <span class="badge badge-grey">conda</span></div>' +
      (c.summary ? '<div class="oc-whatis">' + escHtml(c.summary) + '</div>' : '') + '</div>';
    html += '<div class="field"><label class="field-label">Version</label>' +
      makeComboboxHtml('oc-conda-ver', '', _ocSel.latest ? _ocSel.latest + ' (default)' : 'version') +
      '</div>';
  } else {
    const e = _ocSel.entry;
    html += '<div class="oc-hdr"><div class="oc-title">' + escHtml(e.name) +
      (e.alias ? ' <span class="oc-alias">[' + escHtml(e.alias) + ']</span>' : '') +
      ' ' + _ocBadges(e) + '</div>' +
      (e.whatis ? '<div class="oc-whatis">' + escHtml(e.whatis) + '</div>' : '') + '</div>';
    if (e.is_template) {
      html += '<div class="f-divider">Placeholders</div>';
      for (const k of (e.pl_order || [])) {
        const dflt = _ocSel.defaults[k];
        html += '<div class="field"><label class="field-label">' + escHtml(k) + '</label>' +
          makeComboboxHtml('oc-pl-' + k, '', dflt ? dflt + ' (default)' : 'value') +
          '</div>';
      }
    }
  }
  html += '<div id="oc-inspect"></div>';
  html += '<div class="f-divider">Will install</div>' +
    '<div class="oc-concrete" id="oc-concrete">' + escHtml(_ocConcreteName()) + '</div>';
  gid('oc-detail').innerHTML = html;

  // Wire comboboxes after the HTML is in the DOM.
  if (_ocSel.kind === 'conda') {
    initCombobox('oc-conda-ver', _ocSel.conda.versions || [], false, () => {
      _ocSel.version = gid('oc-conda-ver').value.trim();
      gid('oc-concrete').textContent = _ocConcreteName();
    });
  } else if (_ocSel.entry.is_template) {
    const e = _ocSel.entry;
    for (const k of (e.pl_order || [])) {
      const all  = (e.pl || {})[k] || [];
      const vals = all.filter(v => v !== '*');
      const open = all.includes('*'); // '*' = free-form values allowed
      initCombobox('oc-pl-' + k, vals, open,
        () => ocSetVar(k, gid('oc-pl-' + k).value.trim()));
    }
  }
}

function ocSetVar(k, v) {
  _ocSel.vars[k] = v;
  gid('oc-concrete').textContent = _ocConcreteName();
  _ocScheduleInspect();
}

// _ocScheduleInspect re-checks the concrete name against /api/avail/inspect
// (interactive prompts, scheduler submission) after placeholder changes.
// The Create button stays disabled while the check is in flight.
function _ocScheduleInspect() {
  if (!_ocSel || _ocSel.kind === 'conda') return;
  clearTimeout(_ocInspectTimer);
  _ocInspectRes = null;
  _ocInspectPending = true;
  ocValidate();
  const name = _ocConcreteName();
  const box = gid('oc-inspect');
  if (box) box.innerHTML = '<div class="td-muted">Checking build script…</div>';
  _ocInspectTimer = setTimeout(async () => {
    let info = null;
    try {
      const r = await fetch('/api/avail/inspect?name=' + encodeURIComponent(name));
      if (r.ok) info = await r.json();
    } catch {}
    if (_ocConcreteName() !== name) return; // selection changed while fetching
    _ocInspectRes = info;
    _ocInspectPending = false;
    _ocRenderInspect(info);
    ocValidate();
  }, 300);
}

function _ocRenderInspect(info) {
  const box = gid('oc-inspect');
  if (!box) return;
  if (!info) { box.innerHTML = ''; return; }
  let html = '';
  if ((info.prompts || []).length) html += '<div class="f-divider">Required input</div>';
  (info.prompts || []).forEach((p, i) => {
    // Runs on escaped text: URLs become clickable, casing stays intact
    // (oc-prompt disables the field-label uppercase transform).
    const label = escHtml(p).replace(/\\\\n|\\n/g, '<br>')
      .replace(/https?:\/\/[^\s<]+/g, u =>
        '<a href="' + u + '" target="_blank" rel="noopener">' + u + '</a>');
    html += '<div class="field"><label class="field-label oc-prompt">' + label + '</label>' +
      '<input class="field-input" id="oc-ans-' + i + '" placeholder="required" ' +
        'oninput="this.dataset.touched=1;ocValidate()"></div>';
  });
  if (info.conda_fallback) {
    html += '<div class="oc-note">' + iconSvg('info') +
      '<span>No build script found for this name — it will be installed from conda via micromamba.</span></div>';
  } else if (info.will_submit) {
    html += '<div class="oc-note">' + iconSvg('info') +
      '<span>This build will be submitted as a ' + escHtml(info.scheduler) + ' job. ' +
      'The overlay appears in the list once the job finishes.</span></div>';
  }
  box.innerHTML = html;
}

// ocValidate enables the Create button only when the selection is ready:
// no inspect in flight and every interactive prompt answered.
function ocValidate() {
  const btn = gid('oc-create-btn');
  if (!btn || !_ocSel) return;
  let ok = true;
  if (_ocSel.kind === 'script') {
    if (_ocInspectPending) ok = false;
    const n = (_ocInspectRes && _ocInspectRes.prompts) ? _ocInspectRes.prompts.length : 0;
    for (let i = 0; i < n; i++) {
      const el = gid('oc-ans-' + i);
      const empty = !el || !el.value.trim();
      if (el) el.classList.toggle('combo-invalid', empty && el.dataset.touched === '1');
      if (empty) ok = false;
    }
  }
  btn.disabled = !ok;
}

async function ocCreate() {
  const name = _ocConcreteName();
  if (!name) return;
  const answers = [];
  const nPrompts = (_ocInspectRes && _ocInspectRes.prompts) ? _ocInspectRes.prompts.length : 0;
  for (let i = 0; i < nPrompts; i++) {
    const el = gid('oc-ans-' + i);
    const v = el ? el.value.trim() : '';
    if (!v) { // required — mark and abort
      if (el) { el.dataset.touched = '1'; el.classList.add('combo-invalid'); el.focus(); }
      return;
    }
    answers.push(v);
  }
  const btn = gid('oc-create-btn');
  btn.disabled = true;
  let r;
  try {
    r = await fetch('/api/create', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ name, answers }),
    });
  } catch (e) {
    btn.disabled = false;
    showProgressError('Create ' + name, String(e));
    return;
  }
  btn.disabled = false;
  if (!r.ok) { showProgressError('Create ' + name, await r.text()); return; }
  const { id } = await r.json();
  closeModal('oc-modal');
  openLogProgress('Create ' + name, id, msg => {
    if (msg.ok && !msg.submitted) loadOverlays();
  });
}
