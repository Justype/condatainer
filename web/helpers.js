/* ── Progress stack ──────────────────────── */
// A stack of independent floating progress "pills" in #prog-stack, one per
// operation. Two kinds:
//   'log' — a collapsible scrolling TermView, for operations with text
//           output (helper starts, overlay create/edit, script updates).
//   'bar' — a compact percentage bar, for "N of M" operations (delete,
//           move, upload).
let _progPillSeq = 0;

// _createPill builds a pill. Its close button is a Cancel button while the
// task runs (calls onCancel once, then disables) and becomes a plain
// dismiss action once the pill finishes; pass onCancel=null for a pill
// with nothing to cancel.
function _createPill(title, kind, onCancel) {
  const el = document.createElement('div');
  el.className = 'prog-pill ' + kind + '-type';
  el.id = 'prog-pill-' + (++_progPillSeq);
  el.innerHTML =
    '<div class="prog-pill-hdr">' +
      '<span class="prog-pill-status"></span>' +
      '<span class="prog-pill-title"></span>' +
      (kind === 'log' ? '<span class="prog-pill-chevron">' + iconSvg('expand_circle_down') + '</span>' : '') +
      '<button class="btn btn-ghost btn-sm btn-icon prog-pill-close" title="' + (onCancel ? 'Cancel' : 'Close') + '"><svg class="icon" aria-hidden="true"><use href="#i-close"></use></svg></button>' +
    '</div>' +
    (kind === 'bar'
      ? '<div class="prog-pill-bar-wrap"><div class="prog-pill-track"><div class="prog-pill-fill"></div></div><div class="prog-pill-label"></div></div>'
      : '<div class="prog-pill-log"></div>');
  gid('prog-stack').appendChild(el);

  el.querySelector('.prog-pill-title').textContent = title;
  el.querySelector('.prog-pill-status').innerHTML = iconSvg('progress_activity', 'spin');

  const closeBtn = el.querySelector('.prog-pill-close');
  closeBtn.disabled = !onCancel;
  closeBtn.onclick = () => {
    if (el.classList.contains('finished')) { el.remove(); return; }
    if (onCancel) { closeBtn.disabled = true; onCancel(); }
  };

  if (kind === 'log') {
    el.classList.add('expanded'); // show the log by default while actively running
    el.querySelector('.prog-pill-hdr').addEventListener('click', e => {
      if (e.target.closest('.prog-pill-close')) return;
      const expanded = el.classList.toggle('expanded');
      el.querySelector('.prog-pill-chevron').innerHTML = iconSvg(expanded ? 'expand_circle_down' : 'expand_circle_up');
    });
  }
  return el;
}

// _pillFinish sets a pill's terminal visual state. Successful pills
// auto-remove after a short delay; failures and cancellations stay until
// dismissed.
function _pillFinish(el, ok) {
  el.classList.add('finished');
  el.querySelector('.prog-pill-status').innerHTML = ok ? iconSvg('check') : iconSvg('close', 'icon-danger');
  const closeBtn = el.querySelector('.prog-pill-close');
  closeBtn.disabled = false;
  closeBtn.title = 'Close';
  if (ok) setTimeout(() => el.remove(), 2500);
}

// cancelTask POSTs /api/tasks/{id}/cancel.
function cancelTask(taskId) {
  fetch('/api/tasks/' + encodeURIComponent(taskId) + '/cancel', { method: 'POST' }); //nolint
}

// openLogProgress opens a 'log' pill driven by a server task's SSE stream
// (GET /api/tasks/{id}/stream, {"t":"out"|"done"} events).
function openLogProgress(title, taskId, onDone) {
  const el = _createPill(title, 'log', () => cancelTask(taskId));
  const tv = new TermView(el.querySelector('.prog-pill-log'));

  const es = new EventSource('/api/tasks/' + encodeURIComponent(taskId) + '/stream');
  es.onmessage = e => {
    try {
      const msg = JSON.parse(e.data);
      if (msg.t === 'out') {
        tv.write(msg.d);
      } else if (msg.t === 'done') {
        tv.done(msg.cancelled ? true : msg.ok, msg.cancelled ? 'Cancelled' : msg.err);
        _pillFinish(el, msg.ok);
        es.close();
        if (onDone) onDone(msg);
      }
    } catch {}
  };
  es.onerror = () => {
    tv.done(false, 'Connection lost');
    _pillFinish(el, false);
    es.close();
    if (onDone) onDone({ ok: false, err: 'Connection lost' });
  };
  return el;
}

// openStreamLogProgress opens a 'log' pill the caller feeds directly.
// Returns {write(chunk), done(ok, errMsg)}.
function openStreamLogProgress(title) {
  const el = _createPill(title, 'log', null);
  const tv = new TermView(el.querySelector('.prog-pill-log'));
  let finished = false;
  return {
    write(chunk) {
      if (!finished) tv.write(chunk);
    },
    done(ok, errMsg) {
      if (finished) return;
      finished = true;
      tv.done(ok, errMsg);
      _pillFinish(el, ok);
    },
  };
}

// openBarProgress opens a 'bar' pill. With taskId it follows the task's
// SSE stream ({"t":"progress"|"done"} events) and its close button cancels
// the task; without taskId the caller drives it via the returned
// .update(current,total) / .done(ok,errMsg), and may pass opts.onCancel.
// opts.formatLabel(current,total) customizes the label text (default
// "current / total"). Returns {update, done, el}.
function openBarProgress(title, taskId, onDone, opts) {
  opts = opts || {};
  const formatLabel = opts.formatLabel || ((cur, tot) => cur + ' / ' + tot);
  const el = _createPill(title, 'bar', taskId ? () => cancelTask(taskId) : (opts.onCancel || null));
  const fill  = el.querySelector('.prog-pill-fill');
  const label = el.querySelector('.prog-pill-label');
  label.textContent = 'Starting…';

  function update(current, total) {
    const pct = total > 0 ? Math.min(100, Math.round((current / total) * 100)) : 0;
    fill.style.width = pct + '%';
    label.textContent = formatLabel(current, total) + ' (' + pct + '%)';
  }
  function done(ok, errMsg, cancelled) {
    if (ok) fill.style.width = '100%';
    label.textContent = ok ? 'Done' : (cancelled ? 'Cancelled' : (errMsg || 'Failed'));
    _pillFinish(el, ok);
  }

  if (taskId) {
    const es = new EventSource('/api/tasks/' + encodeURIComponent(taskId) + '/stream');
    es.onmessage = e => {
      try {
        const msg = JSON.parse(e.data);
        if (msg.t === 'progress') update(msg.current, msg.total);
        else if (msg.t === 'done') { done(msg.ok, msg.err, msg.cancelled); es.close(); if (onDone) onDone(msg); }
      } catch {}
    };
    es.onerror = () => { done(false, 'Connection lost'); es.close(); if (onDone) onDone({ ok: false, err: 'Connection lost' }); };
  }

  return { update, done, el };
}

// showProgressError opens an already-failed 'log' pill, for errors that
// happen before a server-side task exists.
function showProgressError(title, msg) {
  const el = _createPill(title, 'log');
  const tv = new TermView(el.querySelector('.prog-pill-log'));
  tv.done(false, msg);
  _pillFinish(el, false);
}

/* ── TermView — terminal renderer ────────── */
class TermView {
  constructor(el) {
    const pre = document.createElement('pre');
    pre.className = 'tv-pre';
    el.appendChild(pre);
    this._pre = pre;
    this._lines = [''];
  }
  write(chunk) {
    // strip ANSI escape codes
    chunk = chunk.replace(/\x1b\[[0-9;]*[a-zA-Z]/g, '');
    for (const ch of chunk) {
      if      (ch === '\r') this._lines[this._lines.length - 1] = '';
      else if (ch === '\n') this._lines.push('');
      else                  this._lines[this._lines.length - 1] += ch;
    }
    this._pre.textContent = this._lines.join('\n');
    this._pre.parentElement.scrollTop = this._pre.parentElement.scrollHeight;
  }
  done(ok, errMsg) {
    if (!ok && errMsg) {
      const div = document.createElement('div');
      div.className = 'tv-err';
      div.textContent = 'ERROR: ' + errMsg;
      this._pre.parentElement.appendChild(div);
    }
  }
}

/* ── Helpers / Start section ─────────────── */
let _helperSavedConfig = {}; // saved config for the currently selected helper
let _gpuOptions = [];        // available GPU types from scheduler
let _gpuOptionsLoaded = false;
let _gpuOptionsLoading = null;
let _helpersLoaded = false;
let _helpersLoading = null;

async function loadGpuOptions(force) {
  if (!force && (_gpuOptionsLoaded || _gpuOptionsLoading)) return _gpuOptionsLoading;
  _gpuOptionsLoading = (async () => {
    try {
      const opts = await fetch('/api/gpu-options').then(r => r.json());
      if (Array.isArray(opts) && opts.length) {
        _gpuOptions = opts;
        initCombobox('cfg-gpu', opts, true, null);
      }
      _gpuOptionsLoaded = true;
    } catch {}
  })();
  try {
    await _gpuOptionsLoading;
  } finally {
    _gpuOptionsLoading = null;
  }
}

async function loadHelpers(force) {
  loadHelperBookmarks();
  if (!force && _helpersLoading) {
    await _helpersLoading;
    return;
  }
  if (!force && _helpersLoaded) {
    renderHelperList(gid('h-search')?.value);
    _applyStartLock();
    if (_formOpen || _opRunning) return;
    if (_startBusy && _startRunHelperName) {
      _setStartBusy(true, iconSvg('play_arrow') + ' Starting…');
      return;
    }
    if (_applyPendingStartSelection(true)) return;
    _restoreSelection();
    loadGpuOptions();
    return;
  }

  _helpersLoading = (async () => {
    try {
      const r = await fetch('/api/helpers/available' + (force ? '?refresh=1' : ''));
      allHelpers = (await r.json()) || [];
      _helpersLoaded = true;
    } catch { allHelpers = []; }
  })();
  try {
    await _helpersLoading;
  } finally {
    _helpersLoading = null;
  }
  if (_startBusy && _startRunHelperName) {
    const h = allHelpers.find(h => h.name === _startRunHelperName);
    if (h) selectedHelper = h;
  }
  renderHelperList(gid('h-search')?.value);
  _applyStartLock();
  if (_formOpen || _opRunning) return;
  if (_startBusy && _startRunHelperName) {
    _setStartBusy(true, iconSvg('play_arrow') + ' Starting…');
    return;
  }
  if (_applyPendingStartSelection(true)) return;
  _restoreSelection();
  loadGpuOptions(force);
}

async function refreshStartHelpers() {
  const btn = gid('start-refresh-btn');
  const old = btn && btn.innerHTML;
  if (btn) {
    btn.disabled = true;
    btn.innerHTML = iconSvg('refresh', 'spin') + ' Refreshing';
  }
  try {
    _helpersLoaded = false;
    await loadHelpers(true);
  } finally {
    if (btn) {
      btn.disabled = false;
      btn.innerHTML = old;
    }
  }
}

async function updateHelperScripts() {
  const btn = gid('start-update-btn');
  if (btn) btn.disabled = true;
  _opRunning = true;
  try {
    const r = await fetch('/api/helpers/update', { method: 'POST' });
    if (!r.ok) {
      const msg = await r.text();
      showProgressError('Updating helpers', msg || String(r.status));
      if (btn) btn.disabled = false;
      _opRunning = false;
      return;
    }
    const { id } = await r.json();
    openLogProgress('Updating helpers', id, async msg => {
      if (btn) btn.disabled = false;
      _opRunning = false;
      if (msg.ok) await refreshStartHelpers();
    });
  } catch (e) {
    if (btn) btn.disabled = false;
    _opRunning = false;
    showProgressError('Updating helpers', String(e));
  }
}

/* ── Selection persistence ───────────────── */
const _selKey = 'cnt_start_selection';
let _pendingStartSelection = null;
function _saveSelection() {
  if (!selectedHelper) return;
  const params = {};
  _helperParamKeys.forEach(key => {
    const el = gid('hparam-' + key);
    if (el && el.value) params[key] = el.value;
  });
  localStorage.setItem(_selKey, JSON.stringify({ name: selectedHelper.name, params }));
}
function _restoreSelection() {
  try {
    const saved = JSON.parse(localStorage.getItem(_selKey) || 'null');
    if (!saved || !saved.name) return;
    if (!allHelpers.find(h => h.name === saved.name)) return;
    selectHelper(saved.name, Object.keys(saved.params || {}).length ? { params: saved.params } : undefined);
  } catch { /* ignore */ }
}

function queueStartSelection(name, overrides) {
  _pendingStartSelection = { name, overrides };
  localStorage.setItem(_selKey, JSON.stringify({ name, params: (overrides && overrides.params) || {} }));
  _applyPendingStartSelection(false);
}

function _applyPendingStartSelection(clear) {
  const pending = _pendingStartSelection;
  if (!pending || !allHelpers.find(h => h.name === pending.name)) return false;
  selectHelper(pending.name, pending.overrides);
  if (clear) _pendingStartSelection = null;
  return true;
}

/* ── Apply saved-config as placeholders ───── */
function _parseParamOpts(pv, key) {
  const values = (pv || {})[key] || [];
  const isOpen = values.length > 0 && values[values.length - 1] === '*';
  const opts   = isOpen ? values.slice(0, -1) : values;
  return { opts, isOpen };
}

function _computeParams(params, pv) {
  return (params || []).map(p => {
    const { opts, isOpen } = _parseParamOpts(pv, p.key);
    return { p, opts, isOpen };
  });
}

function _rdPh(rd) {
  return {
    cpus: rd.cpus ? String(rd.cpus) : '',
    mem:  rd.mem      || '',
    time: rd.walltime || '',
    gpu:  rd.gpu      || '',
  };
}

function _paramPlaceholder(param, opts) {
  if (!param) return '';
  if (param.default) return param.default;
  if (opts && opts.length) return opts[0];
  return param.optional ? '(optional)' : '';
}

function _applyConfigPlaceholders() {
  const cfg = _helperSavedConfig;
  const ph  = _rdPh((selectedHelper && selectedHelper.resource_defaults) || {});
  gid('cfg-cpus').placeholder = cfg.cpus || ph.cpus;
  gid('cfg-mem').placeholder  = cfg.mem  || ph.mem;
  gid('cfg-wall').placeholder = cfg.time || ph.time;
  gid('cfg-gpu').placeholder  = cfg.gpu  || ph.gpu;
  _helperParamKeys.forEach(key => {
    const el = gid('hparam-' + key);
    if (!el) return;
    const param = (selectedHelper?.params || []).find(p => p.key === key);
    const { opts } = _parseParamOpts(selectedHelper?.param_values, key);
    el.placeholder = cfg[key.toLowerCase()] || _paramPlaceholder(param, opts);
  });
  renderModuleChips();
}

/* ── Saved-config gear button ────────────── */
function _renderHelperConfigWrap() {
  const wrap = gid('helper-config-wrap');
  if (!wrap) return;
  if (!selectedHelper) { wrap.innerHTML = ''; return; }
  wrap.innerHTML =
    '<button class="btn btn-ghost btn-sm" onclick="openHelperConfigModal()">' +
    iconSvg('settings') + ' Defaults</button>';
}

function openHelperConfigModal() {
  if (!selectedHelper) return;
  gid('hcfg-title').textContent = 'Saved Defaults — ' + selectedHelper.name;
  const body   = gid('hcfg-body');
  const rd     = selectedHelper.resource_defaults || {};
  const params = selectedHelper.params || [];
  const pv     = selectedHelper.param_values || {};

  const rph = _rdPh(rd);
  const resRows = [
    { key: 'cpus', label: 'CPUs',          ph: rph.cpus },
    { key: 'mem',  label: 'Memory',        ph: rph.mem  },
    { key: 'time', label: 'Walltime',      ph: rph.time },
    { key: 'gpu',  label: 'GPU (optional)', ph: rph.gpu },
  ];

  const paramComputed = _computeParams(params, pv);

  let html = '<div class="param-grid">';

  resRows.forEach(({ key, label, ph }) => {
    const saved = _helperSavedConfig[key] || '';
    const useCombo = key === 'gpu' && _gpuOptions.length > 0;
    html += '<span class="param-name">' + escHtml(label) + '</span>' +
      '<div class="param-cell">' +
        (useCombo
          ? '<div class="combo-wrap" id="combo-hcfg-gpu">' +
              '<input class="field-input combo-input" id="hcfg-gpu" autocomplete="off"' +
              ' value="' + escHtml(saved) + '" data-value="' + escHtml(saved) + '"' +
              ' placeholder="' + escHtml(ph) + '">' +
              '<div class="combo-list" id="combo-list-hcfg-gpu"></div>' +
            '</div>'
          : '<input class="field-input" id="hcfg-' + key + '" value="' + escHtml(saved) +
              '" placeholder="' + escHtml(ph) + '">') +
      '</div>';
  });

  if (params.length) {
    html += '<div class="f-divider" style="grid-column:1/-1;margin:10px 0 2px">Helper params</div>';
  }

  paramComputed.forEach(({ p, opts }) => {
    const saved = _helperSavedConfig[p.key.toLowerCase()] || '';
    const ph    = _paramPlaceholder(p, opts);
    html += '<span class="param-name">' + escHtml(p.key) + '</span>' +
      '<div class="param-cell">' +
        (opts.length
          ? '<div class="combo-wrap" id="combo-hcfg-' + escHtml(p.key) + '">' +
              '<input class="field-input combo-input" id="hcfg-' + escHtml(p.key) + '" autocomplete="off"' +
              ' value="' + escHtml(saved) + '" data-value="' + escHtml(saved) + '"' +
              ' placeholder="' + escHtml(ph) + '">' +
              '<div class="combo-list" id="combo-list-hcfg-' + escHtml(p.key) + '"></div>' +
            '</div>'
          : '<input class="field-input" id="hcfg-' + escHtml(p.key) + '" value="' + escHtml(saved) +
              '" placeholder="' + escHtml(ph) + '">') +
        (p.desc ? '<span class="param-desc">' + escHtml(p.desc) + '</span>' : '') +
      '</div>';
  });

  html += '</div>';
  body.innerHTML = html;

  if (_gpuOptions.length) initCombobox('hcfg-gpu', _gpuOptions, true, null);
  paramComputed.forEach(({ p, opts, isOpen }) => {
    if (opts.length) initCombobox('hcfg-' + p.key, opts, isOpen, null);
  });

  gid('hcfg-modal').classList.add('open');
}

async function saveHelperConfigModal() {
  if (!selectedHelper) return;
  const data = {};
  ['cpus', 'mem', 'time', 'gpu'].forEach(k => {
    const el = gid('hcfg-' + k);
    if (el) data[k] = el.value.trim();
  });
  (selectedHelper.params || []).forEach(p => {
    const el = gid('hcfg-' + p.key);
    if (!el) return;
    data[p.key.toLowerCase()] = (el.dataset.value !== undefined ? el.dataset.value : el.value).trim();
  });
  try {
    const r = await fetch('/api/helpers/' + encodeURIComponent(selectedHelper.name) + '/config', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(data),
    });
    if (!r.ok) return;
    _helperSavedConfig = await r.json();
    _applyConfigPlaceholders();
    closeModal('hcfg-modal');
  } catch {}
}

function resetStartForm(helper) {
  const h = helper || selectedHelper;

  selectedModules          = [];
  selectedExternalOverlays = [];
  renderModuleChips();

  // Clear typed values; placeholders already reflect the defaults.
  gid('cfg-cpus').value = '';
  gid('cfg-mem').value  = '';
  gid('cfg-wall').value = '';
  gid('cfg-gpu').value  = '';
  _setCwd('');
  gid('cfg-cwd').dispatchEvent(new Event('input'));
  _setVal('cfg-overlay', '');
  gid('cfg-overlay').dispatchEvent(new Event('input'));

  _helperParamKeys.forEach(key => {
    const input = gid('hparam-' + key);
    if (!input) return;
    input.value = '';
    if (input.dataset.value !== undefined) input.dataset.value = '';
    input.classList.remove('combo-invalid');
  });

  gid('ov-create-form').classList.remove('visible');
  gid('ov-name').value = '';
  gid('ov-size').value = '';
  gid('ov-pkgs').value = '';
  gid('ov-create-btn').disabled = true;

  const _pv = (h && h.param_values) || {};
  _ovTokenKeys.forEach(key => {
    const input = gid('ovtok-' + key);
    if (!input) return;
    const { opts } = _parseParamOpts(_pv, key);
    const param = (h?.params || []).find(p => p.key === key);
    const value = (param && param.default) || (opts[0] || '');
    input.value = value;
    input.dataset.value = value;
    input.classList.remove('combo-invalid');
  });
  updateOvPreview();
}

/* ── Helper bookmark state (server-persisted, not localStorage) ── */
let helperBookmarks       = new Set();
let _helperBookmarksLoaded = false;

async function loadHelperBookmarks() {
  if (_helperBookmarksLoaded) return;
  _helperBookmarksLoaded = true;
  try {
    const list = (await (await fetch('/api/helpers/bookmarks')).json()) || [];
    helperBookmarks = new Set(list);
  } catch { helperBookmarks = new Set(); }
  renderHelperList(gid('h-search')?.value);
}

async function toggleHelperBookmark(name, e) {
  e.stopPropagation();
  const isBookmarked = helperBookmarks.has(name);
  try {
    const r = isBookmarked
      ? await fetch('/api/helpers/bookmarks?name=' + encodeURIComponent(name), { method: 'DELETE' })
      : await fetch('/api/helpers/bookmarks', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ name })
        });
    if (!r.ok) { alert('Could not update bookmark: ' + (await r.text())); return; }
    helperBookmarks = new Set((await r.json()) || []);
  } catch (err) {
    alert('Could not update bookmark: ' + err);
    return;
  }
  renderHelperList(gid('h-search')?.value);
}

function renderHelperList(filter) {
  const list      = gid('helper-list');
  const q         = (filter || '').toLowerCase().trim();
  const bookmarked = helperBookmarks;
  let helpers  = q
    ? allHelpers.filter(h =>
        h.name.toLowerCase().includes(q) ||
        (h.whatis || '').toLowerCase().includes(q))
    : allHelpers;
  // Bookmarked items float to the top
  helpers = [
    ...helpers.filter(h => bookmarked.has(h.name)),
    ...helpers.filter(h => !bookmarked.has(h.name)),
  ];
  if (!helpers.length) {
    list.innerHTML = '<div class="modal-empty">' +
      (q ? 'No matches.' : 'No helpers found.') + '</div>';
    return;
  }
  list.innerHTML = helpers.map(h => {
    const isBookmarked = bookmarked.has(h.name);
    return '<div class="h-item" data-name="' + escHtml(h.name) + '" onclick="selectHelper(\'' + escHtml(h.name) + '\')">' +
      '<div class="h-item-content">' +
        '<div class="h-item-name">' + escHtml(h.name) + '</div>' +
        '<div class="h-item-desc">'  + escHtml(h.whatis || '') + '</div>' +
      '</div>' +
      '<button class="btn btn-sm btn-ghost h-item-star' + (isBookmarked ? ' starred' : '') + '" onclick="toggleHelperBookmark(\'' + escHtml(h.name) + '\',event)" title="' + (isBookmarked ? 'Remove bookmark' : 'Bookmark') + '">' +
        iconSvg('star', null, isBookmarked) +
      '</button>' +
    '</div>';
  }).join('');
  // Restore active state for currently selected helper
  if (selectedHelper) {
    document.querySelectorAll('.h-item').forEach(e =>
      e.classList.toggle('active', e.dataset.name === selectedHelper.name));
  }
}

function selectHelper(name, overrides) {
  if ((_formOpen || _opRunning) && selectedHelper && selectedHelper.name !== name) return;
  selectedHelper = allHelpers.find(h => h.name === name);
  if (!selectedHelper) return;
  setHash('#start/' + encodeURIComponent(name), false);
  if (!overrides) localStorage.setItem(_selKey, JSON.stringify({ name, params: {} }));
  selectedModules          = overrides?.modules   ?? [];
  selectedExternalOverlays = overrides?.externals ?? [];
  _helperParamKeys = [];
  renderModuleChips();

  document.querySelectorAll('.h-item').forEach(e =>
    e.classList.toggle('active', e.dataset.name === name));
  gid('config-placeholder').style.display = 'none';
  gid('config-panel').classList.add('visible');
  gid('start-btn').innerHTML = iconSvg('play_arrow') + ' Start ' + escHtml(name);
  gid('start-term').classList.remove('visible');
  gid('start-term').innerHTML = '';
  if (_startBusy) _setStartBusy(true, iconSvg('play_arrow') + ' Starting…');
  if (!overrides) {
    ['cfg-cpus', 'cfg-mem', 'cfg-wall', 'cfg-gpu'].forEach(id => {
      const el = gid(id);
      if (!el) return;
      el.value = '';
      if (el.dataset.value !== undefined) el.dataset.value = '';
    });
  }

  // Helper-specific params
  const params    = selectedHelper.params || [];
  _helperParamKeys = params.map(p => p.key);
  const paramWrap  = gid('helper-params-wrap');
  if (params.length) {
    const pv = selectedHelper.param_values || {};
    const computed = _computeParams(params, pv);
    paramWrap.innerHTML =
      '<div class="f-divider">Helper params</div>' +
      '<div class="param-grid">' +
      computed.map(({ p, opts }) => {
        const ph      = _paramPlaceholder(p, opts);
        const isRequired = !p.optional && !p.default;
        const initVal = (p.optional && opts.length) ? opts[0] : '';
        const input  = opts.length
          ? makeComboboxHtml('hparam-' + p.key, initVal, ph)
          : '<input class="field-input" id="hparam-' + escHtml(p.key) + '"' +
              (isRequired ? ' data-required="1"' : '') +
              ' placeholder="' + escHtml(ph) + '">';
        return '<span class="param-name">' + escHtml(p.key) + '</span>' +
          '<div class="param-cell">' +
            input +
            (p.desc ? '<span class="param-desc">' + escHtml(p.desc) + '</span>' : '') +
          '</div>';
      }).join('') +
      '</div>';
    // Wire comboboxes after HTML is in the DOM
    computed.forEach(({ p, opts, isOpen }) => {
      if (opts.length) initCombobox('hparam-' + p.key, opts, isOpen, () => renderModuleChips());
      const isRequired = !p.optional && !p.default;
      const el = gid('hparam-' + p.key);
      if (!el) return;
      el.addEventListener('input', () => { renderModuleChips(); _saveSelection(); });
      if (!isRequired) return;
      el.dataset.required = '1';
      // Plain inputs need their own blur listener; combobox inputs use initCombobox's
      if (!opts.length) {
        el.addEventListener('blur', () =>
          el.classList.toggle('combo-invalid', !el.value.trim()));
      }
    });
    // Restore param values when rerunning
    if (overrides && overrides.params) {
      Object.entries(overrides.params).forEach(([k, v]) => {
        const el = gid('hparam-' + k);
        if (el) { el.value = v; el.dataset.value = v; }
      });
    }
  } else {
    paramWrap.innerHTML = '';
  }

  // Load resource defaults from script headers (#NCPUS:/#MEM:/#TIME:)
  fetch('/api/helpers/' + encodeURIComponent(name) + '/resources')
    .then(r => r.json())
    .then(d => {
      if (!selectedHelper || selectedHelper.name !== name) return;
      selectedHelper.resource_defaults = d || {};
      if (overrides) {
        _setVal('cfg-cpus',  overrides.cpus  != null ? overrides.cpus  : (d.cpus  || ''));
        _setVal('cfg-mem',   overrides.mem   || d.mem      || '');
        _setVal('cfg-wall',  overrides.wall  || d.walltime || '');
        _setVal('cfg-gpu',   overrides.gpu   != null ? overrides.gpu   : (d.gpu      || ''));
      } else {
        // Script defaults loaded — re-apply placeholders (saved config may already be set).
        _applyConfigPlaceholders();
      }
    }).catch(() => {});

  // Load saved config; apply resource + param values and render the gear button.
  _helperSavedConfig = {};
  _renderHelperConfigWrap();
  fetch('/api/helpers/' + encodeURIComponent(name) + '/config')
    .then(r => r.json())
    .then(cfg => {
      if (!selectedHelper || selectedHelper.name !== name) return;
      _helperSavedConfig = cfg || {};
      if (!overrides) {
        _applyConfigPlaceholders();
      }
      _renderHelperConfigWrap();
    }).catch(() => {});

  // Reset create form
  gid('ov-create-form').classList.remove('visible');
  gid('ov-name').value = '';
  gid('ov-create-btn').disabled = true;
  renderOverlayCreateForm();

  // When rerunning, apply saved cwd/overlay directly and skip auto-resolution.
  if (overrides) {
    _setCwd(overrides.cwd || '');
    let overlayPath = overrides.overlay || '';
    if (overlayPath && !overlayPath.startsWith('/')) {
      overlayPath = (overrides.cwd || '').replace(/\/+$/, '') + '/' + overlayPath;
    }
    _setVal('cfg-overlay', overlayPath);
    gid('cfg-overlay').dispatchEvent(new Event('input'));
    renderModuleChips();
    return;
  }

  // Keep a user-selected env overlay when switching helpers; auto-resolve only
  // when the field is blank.
  if (selectedHelper.img_required) {
    gid('ov-pkgs').placeholder = 'additional packages';
  }
  if (gid('cfg-overlay').value) {
    gid('cfg-overlay').dispatchEvent(new Event('input'));
  } else {
    refreshEnvOverlayFromCWD();
  }
}

async function refreshEnvOverlayFromCWD() {
  if (!selectedHelper) return;
  if (gid('cfg-overlay').value) return; // don't overwrite a manually-set or rerun overlay
  const helperName = selectedHelper.name;
  const cwd = gid('cfg-cwd').value || '';
  const url = '/api/helpers/' + encodeURIComponent(helperName) +
    '/find-env?cwd=' + encodeURIComponent(cwd);
  try {
    const d = await fetch(url).then(r => r.json());
    if (!selectedHelper || selectedHelper.name !== helperName) return;
    _setVal('cfg-overlay', d.path || '');
  } catch {
    _setVal('cfg-overlay', '');
  }
  gid('cfg-overlay').dispatchEvent(new Event('input'));
}

function clearEnvOverlay() {
  _setVal('cfg-overlay', '');
  gid('cfg-overlay').dispatchEvent(new Event('input'));
}

gid('cfg-overlay').addEventListener('input', function () {
  const hasOverlay = !!this.value;
  gid('ov-clear-btn').disabled = !hasOverlay;
  _scheduleOvCheck(this.value.trim());

  const notice  = gid('ov-notice');
  const startBtn = gid('start-btn');
  const required = !!(selectedHelper && selectedHelper.img_required);

  if (!hasOverlay && selectedHelper) {
    notice.classList.add('visible');
    if (required) {
      notice.classList.add('error');
      notice.classList.remove('info');
      gid('ov-notice-icon').setAttribute('href', '#i-error');
      gid('ov-notice-text').textContent = 'No env overlay selected — required by ' + selectedHelper.name + ' helper';
      if (startBtn) startBtn.disabled = true;
    } else {
      notice.classList.add('info');
      notice.classList.remove('error');
      gid('ov-notice-icon').setAttribute('href', '#i-info');
      gid('ov-notice-text').textContent = 'No env overlay selected — helpful for storing packages and conda env';
      if (startBtn) startBtn.disabled = _startBusy || _formOpen;
    }
  } else {
    notice.classList.remove('visible', 'error', 'info');
    if (startBtn) startBtn.disabled = _startBusy || _formOpen;
  }
});

gid('cfg-cwd').addEventListener('input', refreshEnvOverlayFromCWD);

function _addExternalOverlay(path) {
  if (path && !selectedExternalOverlays.includes(path)) {
    selectedExternalOverlays.push(path);
    renderModuleChips();
  }
}

function _collectParams() {
  const out = {};
  _helperParamKeys.forEach(k => {
    const e = gid('hparam-' + k);
    if (!e) return;
    const v = (e.dataset.value !== undefined ? e.dataset.value : e.value).trim();
    if (v) out[k] = v;
  });
  return out;
}

/* ── Overlay create form ─────────────────── */
let _ovTokenKeys = [];
let _ovExists    = false; // true when cfg-overlay points to an existing .img
let _ovWritable  = false; // true when that .img is writable

/* ── Overlay state check (New vs Edit button) ── */
let _ovCheckTimer = null;
function _scheduleOvCheck(path) {
  clearTimeout(_ovCheckTimer);
  _ovCheckTimer = setTimeout(() => _checkOverlayState(path), 300);
}

async function _checkOverlayState(path) {
  const btn = gid('ov-new-btn');
  if (!path) {
    _ovExists = false; _ovWritable = false;
    btn.textContent = 'New';
    btn.onclick = toggleOverlayCreate;
    gid('ov-edit-form').classList.remove('visible');
    // Don't unfreeze if the create form is still open (user opened it intentionally).
    if (!gid('ov-create-form').classList.contains('visible')) {
      _setFormOpen(false);
    }
    return;
  }
  try {
    const r = await fetch('/api/env/check?path=' + encodeURIComponent(path));
    if (!r.ok) return;
    const st = await r.json();
    _ovExists   = !!st.exists;
    _ovWritable = !!st.writable;
    if (_ovExists && _ovWritable) {
      btn.textContent = 'Edit';
      btn.onclick = toggleOverlayEdit;
      if (gid('ov-edit-form').classList.contains('visible')) {
        _loadOverlayInfo();
      }
      if (gid('ov-create-form').classList.contains('visible')) {
        gid('ov-create-form').classList.remove('visible');
        _setFormOpen(false);
      }
    } else {
      btn.textContent = 'New';
      btn.onclick = toggleOverlayCreate;
      if (gid('ov-edit-form').classList.contains('visible')) {
        gid('ov-edit-form').classList.remove('visible');
        _setFormOpen(false);
      }
    }
  } catch (_) {}
}

/* ── Overlay edit form ───────────────────── */

// Shared busy state for when create/edit forms are open.
// Locks the helper list and start button; on close, re-evaluates start button
// state via the cfg-overlay input event (handles required-overlay validation).
function _setFormOpen(open) {
  _formOpen = open;
  _applyStartLock();
  if (!open) gid('cfg-overlay').dispatchEvent(new Event('input'));
}

function toggleOverlayEdit() {
  if (_startBusy || _opRunning) return;
  const form = gid('ov-edit-form');
  const wasHidden = !form.classList.contains('visible');
  form.classList.toggle('visible');
  _setFormOpen(wasHidden);
  if (wasHidden) _loadOverlayInfo();
}

async function _loadOverlayInfo() {
  const path = gid('cfg-overlay').value.trim();
  if (!path) return;
  try {
    const r = await fetch('/api/env/info?path=' + encodeURIComponent(path));
    if (!r.ok) return;
    const info = await r.json();
    const curSize = gid('ov-edit-cur-size');
    if (curSize && info.size_mb) curSize.textContent = '(current: ' + info.size_mb + ' MiB)';
    const pkgWrap = gid('ov-edit-pkgs-wrap');
    const pkgList = gid('ov-edit-pkg-list');
    if (info.specs && info.specs.length) {
      pkgWrap.style.display = '';
      pkgList.innerHTML = info.specs.map(s => '<span>' + escHtml(s) + '</span>').join(' ');
    } else {
      pkgWrap.style.display = 'none';
    }
  } catch (_) {}
}

async function editOverlay() {
  if (_startBusy || _opRunning) return;
  const path       = gid('cfg-overlay').value.trim();
  const size       = gid('ov-edit-size').value.trim();
  const addPkgs    = gid('ov-edit-add').value.trim();
  const removePkgs = gid('ov-edit-remove').value.trim();
  if (!path) return;
  const title = 'Editing ' + path.split('/').pop();
  if (!size && !addPkgs && !removePkgs) {
    showProgressError(title, 'Nothing to change.');
    return;
  }

  const editBtn = gid('ov-edit-btn');
  editBtn.disabled = true;
  _opRunning = true;
  _applyStartLock();
  try {
    const r = await fetch('/api/overlay/edit', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ path, size, add_packages: addPkgs, remove_packages: removePkgs }),
    });
    if (!r.ok) {
      const e = await r.json().catch(() => ({}));
      showProgressError(title, e.detail || e.error || String(r.status));
      _opRunning = false;
      _applyStartLock();
      editBtn.disabled = false;
      return;
    }
    const { id } = await r.json();
    openLogProgress(title, id, async msg => {
      _opRunning = false;
      _applyStartLock();
      editBtn.disabled = false;
      if (msg.ok) {
        gid('ov-edit-size').value   = '';
        gid('ov-edit-add').value    = '';
        gid('ov-edit-remove').value = '';
        gid('ov-edit-form').classList.remove('visible');
        _setFormOpen(false);
      }
    });
  } catch (e) {
    _opRunning = false;
    _applyStartLock();
    editBtn.disabled = false;
    showProgressError(title, String(e));
  }
}

function toggleOverlayCreate() {
  if (_startBusy || _opRunning) return;
  const form = gid('ov-create-form');
  const wasHidden = !form.classList.contains('visible');
  form.classList.toggle('visible');
  _setFormOpen(wasHidden);
  if (wasHidden) {
    const nameEl = gid('ov-name');
    if (!nameEl.value.trim()) {
      nameEl.value = 'env.img';
      gid('ov-create-btn').disabled = false;
    }
  }
}

// Renders token inputs (comboboxes/text) driven by img_packages {KEY} tokens.
function renderOverlayCreateForm() {
  _ovTokenKeys = [];
  const wrap = gid('ov-token-wrap');
  const prevWrap = gid('ov-pkg-preview-wrap');
  if (!selectedHelper) {
    wrap.innerHTML = '';
    prevWrap.style.display = 'none';
    return;
  }
  if (!selectedHelper.img_packages) {
    wrap.innerHTML = '';
    prevWrap.style.display = 'none';
    return;
  }

  const tokens = extractTokens(selectedHelper.img_packages);
  _ovTokenKeys = tokens;

  if (!tokens.length) {
    wrap.innerHTML = '';
    prevWrap.style.display = 'block';
    gid('ov-pkgs').addEventListener('input', updateOvPreview);
    updateOvPreview();
    return;
  }

  const paramValues = selectedHelper.param_values || {};
  const params      = selectedHelper.params || [];

  const computed = tokens.map(key => {
    const { opts, isOpen } = _parseParamOpts(paramValues, key);
    const param  = params.find(p => p.key === key);
    return { key, opts, isOpen, param };
  });

  wrap.innerHTML = '<div class="form-grid">' +
    computed.map(({ key, opts, param }) => {
      const defVal = (param && param.default) || (opts[0] || '');
      const desc   = param ? param.desc : '';
      return '<div class="field">' +
        '<label class="field-label">' + escHtml(key) + '</label>' +
        (opts.length
          ? makeComboboxHtml('ovtok-' + key, '', defVal)
          : '<input class="field-input" id="ovtok-' + escHtml(key) + '" placeholder="' + escHtml(defVal) + '">') +
        (desc ? '<span class="param-desc">' + escHtml(desc) + '</span>' : '') +
      '</div>';
    }).join('') +
  '</div>';

  // Wire comboboxes and inputs to preview updates
  computed.forEach(({ key, opts, isOpen }) => {
    if (opts.length) {
      initCombobox('ovtok-' + key, opts, isOpen, updateOvPreview);
    } else {
      const el = gid('ovtok-' + key);
      if (el) el.addEventListener('input', updateOvPreview);
    }
  });

  prevWrap.style.display = 'block';
  gid('ov-pkgs').addEventListener('input', updateOvPreview);
  updateOvPreview();
}

function extractTokens(tmpl) {
  const seen = new Set(), result = [];
  let m; const re = /\{([A-Z][A-Z0-9_]*)\}/g;
  while ((m = re.exec(tmpl)) !== null) {
    if (!seen.has(m[1])) { seen.add(m[1]); result.push(m[1]); }
  }
  return result;
}

function updateOvPreview() {
  if (!selectedHelper) return;
  let s = selectedHelper.img_packages || '';
  _ovTokenKeys.forEach(k => {
    const el       = gid('ovtok-' + k);
    const v        = el ? (el.dataset.value !== undefined ? el.dataset.value : el.value) : '';
    const fallback = el ? (el.placeholder || ('{' + k + '}')) : ('{' + k + '}');
    s = s.replaceAll('{' + k + '}', v || fallback);
  });
  const extra = (gid('ov-pkgs')?.value || '').trim();
  if (extra) s += ' ' + extra;
  gid('ov-pkg-preview').textContent = s;
}

/* ── Searchable combobox ─────────────────── */
// If placeholder is provided, the input starts empty and shows placeholder as hint.
// Otherwise defVal is used as both value and placeholder (existing behaviour).
function makeComboboxHtml(id, defVal, placeholder) {
  const val = placeholder !== undefined ? '' : defVal;
  const ph  = placeholder !== undefined ? placeholder : defVal;
  return '<div class="combo-wrap" id="combo-' + escHtml(id) + '">' +
    '<input class="field-input combo-input" id="' + escHtml(id) + '" autocomplete="off"' +
      ' value="' + escHtml(val) + '" data-value="' + escHtml(val) + '" placeholder="' + escHtml(ph) + '">' +
    '<div class="combo-list" id="combo-list-' + escHtml(id) + '"></div>' +
  '</div>';
}

function initCombobox(id, options, isOpen, onChange) {
  const input = gid(id);
  const list  = gid('combo-list-' + id);
  if (!input || !list) return;

  function renderList(filter) {
    const f = filter.toLowerCase();
    const matches = options.filter(o => o.toLowerCase().includes(f));
    if (!matches.length && !isOpen) {
      list.innerHTML = '<div class="combo-item combo-empty">No matches</div>';
    } else {
      list.innerHTML = matches.map(o =>
        '<div class="combo-item" data-val="' + escHtml(o) + '">' + escHtml(o) + '</div>'
      ).join('');
    }
    list.querySelectorAll('.combo-item[data-val]').forEach(el =>
      el.addEventListener('mousedown', e => {
        e.preventDefault();
        setValue(el.dataset.val);
        list.classList.remove('open');
      })
    );
  }

  function setValue(v) {
    input.value = v;
    input.dataset.value = v;
    if (onChange) onChange();
  }

  input.addEventListener('focus', () => {
    renderList(input.value);
    list.classList.add('open');
  });
  input.addEventListener('input', () => {
    input.dataset.value = input.value;
    renderList(input.value);
    list.classList.add('open');
    if (onChange) onChange();
  });
  input.addEventListener('blur', () => {
    setTimeout(() => list.classList.remove('open'), 150);
    // Only mark invalid when explicitly required and empty
    if (input.dataset.required && !input.value.trim()) {
      input.classList.add('combo-invalid');
    } else {
      input.classList.remove('combo-invalid');
    }
  });
  input.addEventListener('keydown', e => {
    if (e.key === 'Escape') { list.classList.remove('open'); input.blur(); }
    if (e.key === 'Enter') {
      const first = list.querySelector('.combo-item[data-val]');
      if (first) { setValue(first.dataset.val); list.classList.remove('open'); }
    }
  });

  // Initialize data-value
  input.dataset.value = input.value;
}

/* ── Shared busy state ───────────────────── */
let _startBusy  = false;
let _startBusyLabel = '';
let _startRunHelperName = '';
let _formOpen = false;
let _opRunning  = false; // true only during an active create/edit/start operation
window.addEventListener('beforeunload', e => {
  if (_opRunning) e.preventDefault();
});

function _setStartBusy(busy, label) {
  _startBusy = busy;
  _startBusyLabel = busy ? (label || _startBusyLabel) : '';
  _applyStartLock();
}

function _applyStartLock() {
  const busy = _startBusy || _formOpen || _opRunning;
  const btn  = gid('start-btn');
  const wrap = gid('helper-list-wrap');
  if (btn) {
    btn.disabled = busy;
    btn.innerHTML = _startBusy && _startBusyLabel
      ? _startBusyLabel
      : iconSvg('play_arrow') + ' Start ' + escHtml(selectedHelper?.name || '');
  }
  if (wrap) wrap.classList.toggle('locked', busy);

  const overlayBtn = gid('ov-new-btn');
  if (overlayBtn) overlayBtn.disabled = _startBusy || _opRunning;
  const cwdBrowseBtn = gid('cfg-cwd-browse-btn');
  if (cwdBrowseBtn) cwdBrowseBtn.disabled = busy;
  const overlayBrowseBtn = gid('cfg-overlay-browse-btn');
  if (overlayBrowseBtn) overlayBrowseBtn.disabled = busy;
  const clearBtn = gid('ov-clear-btn');
  if (clearBtn) clearBtn.disabled = busy || !gid('cfg-overlay')?.value;
  const createBtn = gid('ov-create-btn');
  if (createBtn && (_startBusy || _opRunning)) createBtn.disabled = true;
  const editBtn = gid('ov-edit-btn');
  if (editBtn && (_startBusy || _opRunning)) editBtn.disabled = true;
  if (editBtn && !_startBusy && !_opRunning) editBtn.disabled = false;
}

/* ── Create overlay ──────────────────────── */
async function createOverlay() {
  if (_startBusy || _opRunning) return;
  const name = gid('ov-name').value || 'my-env.img';
  const size  = gid('ov-size').value || '20G';

  // Validate closed comboboxes — non-empty values must be from the list.
  for (const key of _ovTokenKeys) {
    const paramValues = (selectedHelper && selectedHelper.param_values) || {};
    const values = paramValues[key] || [];
    const isOpen = values.length > 0 && values[values.length - 1] === '*';
    const opts   = isOpen ? values.slice(0, -1) : values;
    if (!isOpen && opts.length) {
      const el = gid('ovtok-' + key);
      const v  = el ? el.value : '';
      if (v && !opts.includes(v)) {
        el && el.classList.add('combo-invalid');
        showProgressError('Creating ' + name, key + ': select a value from the list.');
        return;
      }
    }
  }

  // Resolve {KEY} tokens in img_packages; empty → fall back to placeholder (default).
  const params = _collectParams();
  let basePkgs = selectedHelper?.img_packages || '';
  _ovTokenKeys.forEach(k => {
    const el       = gid('ovtok-' + k);
    const v        = el ? (el.dataset.value !== undefined ? el.dataset.value : el.value) : '';
    const fallback = el ? (el.placeholder || ('{' + k + '}')) : ('{' + k + '}');
    basePkgs = basePkgs.replaceAll('{' + k + '}', v || fallback);
  });
  // Also resolve helper params ({KEY} in img_packages that come from #PARAM:)
  for (const [k, v] of Object.entries(params)) {
    basePkgs = basePkgs.replaceAll('{' + k + '}', v);
  }

  const extraPkgs = gid('ov-pkgs').value;
  const body = {
    name,
    size,
    cwd:              gid('cfg-cwd').value.trim(),
    base_packages:    basePkgs,
    extra_packages:   extraPkgs,
    post_install_cmd: selectedHelper?.post_install_cmd || '',
  };

  const createUseBtn = gid('ov-create-btn');
  const _setBusy = busy => {
    _opRunning = busy;
    _applyStartLock();
    if (createUseBtn) createUseBtn.disabled = busy;
  };

  const createTitle = 'Creating ' + name;
  _setBusy(true);
  try {
    const r = await fetch('/api/overlay/create', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(body),
    });
    if (!r.ok) {
      const e = await r.json().catch(() => ({}));
      showProgressError(createTitle, e.detail || e.error || String(r.status));
      _setBusy(false);
      return;
    }
    const { id } = await r.json();
    openLogProgress(createTitle, id, async msg => {
      _setBusy(false);
      if (msg.ok && msg.path) {
        _setVal('cfg-overlay', msg.path);
        gid('cfg-overlay').dispatchEvent(new Event('input'));
        await loadOverlays();
        // collapse create form
        gid('ov-create-form').classList.remove('visible');
        _setFormOpen(false);
      }
    });
  } catch (e) {
    _setBusy(false);
    showProgressError(createTitle, String(e));
  }
}

/* ── Start helper job ────────────────────── */
async function startHelper() {
  if (!selectedHelper) return;
  if (selectedHelper.img_required && !gid('cfg-overlay').value) return;
  // Validate required params
  let hasEmpty = false;
  _helperParamKeys.forEach(key => {
    const el = gid('hparam-' + key);
    if (!el || !el.dataset.required) return;
    const empty = !el.value.trim();
    el.classList.toggle('combo-invalid', empty);
    if (empty) hasEmpty = true;
  });
  if (hasEmpty) return;
  const _rv  = id => { const e = gid(id); return e.value.trim() || e.placeholder; };
  const _mem  = v => /^\d+$/.test(v) ? v + 'G' : v;
  const _wall = v => /^\d+$/.test(v) ? v + ':00:00' : /^\d+:\d+$/.test(v) ? v + ':00' : v;
  const body = {
    cpus:    parseInt(_rv('cfg-cpus')) || 0,
    mem:     _mem(_rv('cfg-mem'))  || '',
    time:    _wall(_rv('cfg-wall')) || '',
    gpu:     gid('cfg-gpu').value.trim() || gid('cfg-gpu').placeholder || '',
    cwd:     gid('cfg-cwd').value || '',
    overlay: gid('cfg-overlay').value || '',
    overlays: [...selectedModules.map(m => m.path), ...selectedExternalOverlays],
    params:  _collectParams(),
  };

  const term = openStreamLogProgress('Starting ' + selectedHelper.name);
  _opRunning = true;
  _startRunHelperName = selectedHelper.name;
  _setStartBusy(true, iconSvg('play_arrow') + ' Starting…');

  try {
    const r = await fetch('/api/helpers/' + encodeURIComponent(selectedHelper.name) + '/start', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(body),
    });

    if (!r.ok) {
      const e = await r.json().catch(() => ({}));
      let msg = 'Server error (' + r.status + ')';
      if (e.error === 'env_in_use') {
        msg = 'Env overlay is currently in use by another running instance. Stop it first, or choose a different overlay.';
      } else if (e.error === 'singleton_running') {
        msg = 'A ' + selectedHelper.name + ' instance is already running.';
      } else if (e.detail) {
        msg = e.detail;
      } else if (e.error) {
        msg = e.error;
      }
      term.done(false, msg);
      return;
    }

    const reader = r.body.getReader();
    const dec    = new TextDecoder();
    let hasError = false;
    let output = '';
    while (true) {
      const { done, value } = await reader.read();
      if (done) break;
      const chunk = dec.decode(value, { stream: true });
      output += chunk;
      if (output.includes('ERROR:')) hasError = true;
      term.write(chunk);
    }
    const rest = dec.decode();
    if (rest) {
      output += rest;
      if (output.includes('ERROR:')) hasError = true;
      term.write(rest);
    }
    term.done(!hasError);
    await loadJobs();
    if (!hasError) {
      resetStartForm(selectedHelper);
      if (gid('sec-start')?.classList.contains('active'))
        setTimeout(() => navigate('jobs'), 500);
    }
  } catch (e) {
    term.done(false, String(e));
  } finally {
    _opRunning = false;
    _startRunHelperName = '';
    _setStartBusy(false);
  }
}

/* ── Module chips ────────────────────────── */
function renderModuleChips() {
  const list = gid('cfg-modules-list');
  if (!list) return;

  // Required overlays from the helper script — greyed out, not removable.
  // Resolve {KEY} tokens using current param values (or defaults).
  const params = _collectParams();
  const pv     = selectedHelper?.param_values || {};
  (selectedHelper?.params || []).forEach(p => {
    if (!(p.key in params)) {
      const vals = pv[p.key] || [];
      params[p.key] = _helperSavedConfig[p.key.toLowerCase()] || p.default || (vals[0] || '');
    }
  });
  const requiredNames = (selectedHelper?.required_overlays || '').trim().split(/\s+/).filter(Boolean);
  const requiredChips = requiredNames.map(name => {
    let resolved = name;
    Object.entries(params).forEach(([k, v]) => { resolved = resolved.replaceAll('{' + k + '}', v); });
    return '<div class="module-chip required">' +
      '<span class="module-chip-name">' + escHtml(resolved) + '</span>' +
      '<span class="module-chip-label">default</span>' +
    '</div>';
  }).join('');

  const userChips = selectedModules.map((m, i) =>
    '<div class="module-chip">' +
      '<span class="module-chip-name">' + escHtml(m.name) + '</span>' +
      '<button class="btn btn-ghost btn-sm btn-icon" onclick="removeModule(' + i + ')">' + iconSvg('close') + '</button>' +
    '</div>'
  ).join('');

  const extChips = selectedExternalOverlays.map((p, i) => {
    const label = p.split('/').pop();
    return '<div class="module-chip" title="' + escHtml(p) + '">' +
      iconSvg('draft') +
      '<span class="module-chip-name">' + escHtml(label) + '</span>' +
      '<button class="btn btn-ghost btn-sm btn-icon" onclick="removeExternalOverlay(' + i + ')">' + iconSvg('close') + '</button>' +
    '</div>';
  }).join('');

  list.innerHTML = requiredChips + userChips + extChips;
}

function removeModule(i) {
  selectedModules.splice(i, 1);
  renderModuleChips();
}

function removeExternalOverlay(i) {
  selectedExternalOverlays.splice(i, 1);
  renderModuleChips();
}
