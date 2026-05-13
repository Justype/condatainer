/* ── Progress window ─────────────────────── */
let _progES = null;

function openProgressWindow(title, taskId, onDone) {
  gid('prog-win-title').textContent = title;
  gid('prog-win-status').textContent = '⏳';
  gid('prog-win-close').disabled = true;
  gid('prog-window').classList.add('visible');

  const body = gid('prog-win-body');
  body.innerHTML = '';
  const tv = new TermView(body);

  if (_progES) { _progES.close(); _progES = null; }
  const es = new EventSource('/api/tasks/' + encodeURIComponent(taskId) + '/stream');
  _progES = es;

  es.onmessage = e => {
    try {
      const msg = JSON.parse(e.data);
      if (msg.t === 'out') {
        tv.write(msg.d);
      } else if (msg.t === 'done') {
        tv.done(msg.ok, msg.err);
        gid('prog-win-status').textContent = msg.ok ? '✓' : '✗';
        gid('prog-win-close').disabled = false;
        es.close(); _progES = null;
        if (onDone) onDone(msg);
      }
    } catch {}
  };
  es.onerror = () => {
    tv.done(false, 'Connection lost');
    gid('prog-win-status').textContent = '✗';
    gid('prog-win-close').disabled = false;
    es.close(); _progES = null;
  };
}

function closeProgressWindow() {
  if (_progES) { _progES.close(); _progES = null; }
  gid('prog-window').classList.remove('visible');
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
    if (!ok && errMsg) this.write('\nERROR: ' + errMsg + '\n');
  }
}

/* ── Helpers / Start section ─────────────── */
async function loadHelpers() {
  try {
    const r = await fetch('/api/helpers/available');
    allHelpers = (await r.json()) || [];
  } catch { allHelpers = []; }
  renderHelperList();
}

function _defaultStartResources(helper) {
  const defaults = helper && helper.resource_defaults ? helper.resource_defaults : {};
  return {
    cpus: defaults.cpus || 4,
    mem: defaults.mem || '16G',
    walltime: defaults.walltime || '08:00:00',
  };
}

function resetStartForm(helper) {
  const h = helper || selectedHelper;
  const defaults = _defaultStartResources(h);

  selectedModules = [];
  renderModuleChips();

  _setVal('cfg-cpus', defaults.cpus);
  _setVal('cfg-mem', defaults.mem);
  _setVal('cfg-wall', defaults.walltime);
  _setVal('cfg-gpu', '');
  _setVal('cfg-cwd', '');
  _setVal('cfg-overlay', '');
  gid('cfg-overlay').dispatchEvent(new Event('input'));

  _helperParamKeys.forEach(key => {
    const input = gid('hparam-' + key);
    if (!input) return;
    const param = (h?.params || []).find(p => p.key === key);
    const value = param?.default || '';
    input.value = value;
    input.dataset.value = value;
    input.classList.remove('combo-invalid');
  });

  gid('ov-create-form').classList.remove('visible');
  gid('ov-name').value = '';
  gid('ov-size').value = '8G';
  gid('ov-pkgs').value = '';
  gid('ov-create-btn').disabled = true;

  _ovTokenKeys.forEach(key => {
    const input = gid('ovtok-' + key);
    if (!input) return;
    const paramValues = (h && h.param_values) || {};
    const values = paramValues[key] || [];
    const isOpen = values.length > 0 && values[values.length - 1] === '*';
    const opts = isOpen ? values.slice(0, -1) : values;
    const param = (h?.params || []).find(p => p.key === key);
    const value = (param && param.default) || (opts[0] || '');
    input.value = value;
    input.dataset.value = value;
    input.classList.remove('combo-invalid');
  });
  updateOvPreview();
}

function renderHelperList(filter) {
  const list = gid('helper-list');
  const q = (filter || '').toLowerCase().trim();
  const helpers = q
    ? allHelpers.filter(h =>
        h.name.toLowerCase().includes(q) ||
        (h.whatis || '').toLowerCase().includes(q))
    : allHelpers;
  if (!helpers.length) {
    list.innerHTML = '<div style="padding:20px;color:var(--muted);font-size:12px;text-align:center;">' +
      (q ? 'No matches.' : 'No helpers found.') + '</div>';
    return;
  }
  list.innerHTML = helpers.map(h =>
    '<div class="h-item" data-name="' + escHtml(h.name) + '" onclick="selectHelper(\'' + escHtml(h.name) + '\')">' +
      '<div class="h-item-name">' + escHtml(h.name) + '</div>' +
      '<div class="h-item-desc">'  + escHtml(h.whatis || '') + '</div>' +
    '</div>'
  ).join('');
}

function selectHelper(name, overrides) {
  selectedHelper = allHelpers.find(h => h.name === name);
  if (!selectedHelper) return;
  selectedModules  = overrides && overrides.modules ? overrides.modules : [];
  _helperParamKeys = [];
  renderModuleChips();

  document.querySelectorAll('.h-item').forEach(e =>
    e.classList.toggle('active', e.dataset.name === name));
  gid('config-placeholder').style.display = 'none';
  gid('config-panel').classList.add('visible');
  gid('start-btn').innerHTML = iconSvg('play_arrow') + ' Start ' + escHtml(name);
  gid('start-term').classList.remove('visible');
  gid('start-term').innerHTML = '';

  // Helper-specific params
  const params    = selectedHelper.params || [];
  _helperParamKeys = params.map(p => p.key);
  const paramWrap  = gid('helper-params-wrap');
  if (params.length) {
    const pv = selectedHelper.param_values || {};
    paramWrap.innerHTML =
      '<div class="f-divider">Helper params</div>' +
      params.map(p => {
        const values = pv[p.key] || [];
        const isOpen = values.length > 0 && values[values.length - 1] === '*';
        const opts   = isOpen ? values.slice(0, -1) : values;
        const defVal = p.default || (opts[0] || '');
        const input  = opts.length
          ? makeComboboxHtml('hparam-' + p.key, defVal)
          : '<input class="field-input" id="hparam-' + escHtml(p.key) + '"' +
              ' placeholder="' + escHtml(p.default || '') + '">';
        return '<div class="param-row">' +
          '<span class="param-name">' + escHtml(p.key) + '</span>' +
          '<div class="param-cell">' +
            input +
            (p.desc ? '<span class="param-desc">' + escHtml(p.desc) + '</span>' : '') +
          '</div>' +
        '</div>';
      }).join('');
    // Wire comboboxes after HTML is in the DOM
    params.forEach(p => {
      const values = pv[p.key] || [];
      const isOpen = values.length > 0 && values[values.length - 1] === '*';
      const opts   = isOpen ? values.slice(0, -1) : values;
      if (opts.length) initCombobox('hparam-' + p.key, opts, isOpen, null);
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
      selectedHelper.resource_defaults = d || {};
      if (overrides) {
        _setVal('cfg-cpus',  overrides.cpus  != null ? overrides.cpus  : (d.cpus     || 4));
        _setVal('cfg-mem',   overrides.mem   || d.mem   || '16G');
        _setVal('cfg-wall',  overrides.wall  || d.walltime || '08:00:00');
        _setVal('cfg-gpu',   overrides.gpu   != null ? overrides.gpu   : (d.gpu      || ''));
      } else {
        if (d.cpus)     _setVal('cfg-cpus', d.cpus);
        if (d.mem)      _setVal('cfg-mem',  d.mem);
        if (d.walltime) _setVal('cfg-wall', d.walltime);
        if (d.gpu)      _setVal('cfg-gpu',  d.gpu);
      }
    }).catch(() => {});

  // Reset create form
  gid('ov-create-form').classList.remove('visible');
  gid('ov-name').value = '';
  gid('ov-create-btn').disabled = true;
  renderOverlayCreateForm();

  // When rerunning, apply saved cwd/overlay directly and skip auto-resolution.
  if (overrides) {
    _setVal('cfg-cwd',     overrides.cwd     || '');
    _setVal('cfg-overlay', overrides.overlay || '');
    gid('cfg-overlay').dispatchEvent(new Event('input'));
    renderModuleChips();
    return;
  }

  // Keep a user-selected env overlay when switching helpers; auto-resolve only
  // when the field is blank.
  if (selectedHelper.img_packages) {
    gid('ov-pkgs').placeholder = selectedHelper.img_packages;
    if (gid('cfg-overlay').value) {
      gid('cfg-overlay').dispatchEvent(new Event('input'));
    } else {
      refreshEnvOverlayFromCWD();
    }
  } else {
    gid('cfg-overlay').dispatchEvent(new Event('input'));
  }
}

async function refreshEnvOverlayFromCWD() {
  if (!selectedHelper) return;
  const helperName = selectedHelper.name;
  const cwd = gid('cfg-cwd').value || '';
  const url = '/api/helpers/' + encodeURIComponent(helperName) +
    '/find-env?cwd=' + encodeURIComponent(cwd);
  try {
    const d = await fetch(url).then(r => r.json());
    if (!selectedHelper || selectedHelper.name !== helperName) return;
    if (d.path) {
      _setVal('cfg-overlay', d.path);
      gid('cfg-overlay').dispatchEvent(new Event('input'));
    } else if (selectedHelper.img_packages) {
      gid('ov-notice').classList.add('visible');
    }
  } catch {
    if (selectedHelper.img_packages) gid('ov-notice').classList.add('visible');
  }
}

function clearEnvOverlay() {
  _setVal('cfg-overlay', '');
  gid('cfg-overlay').dispatchEvent(new Event('input'));
}

gid('cfg-overlay').addEventListener('input', function () {
  const hasOverlay = !!this.value;
  gid('ov-clear-btn').disabled = !hasOverlay;
  gid('ov-notice').classList.toggle(
    'visible',
    !!(selectedHelper && selectedHelper.img_packages && !hasOverlay)
  );
});

gid('cfg-cwd').addEventListener('input', refreshEnvOverlayFromCWD);

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

function toggleOverlayCreate() {
  const form = gid('ov-create-form');
  form.classList.toggle('visible');
}

// Renders token inputs (comboboxes/text) driven by img_packages {KEY} tokens.
function renderOverlayCreateForm() {
  _ovTokenKeys = [];
  const wrap = gid('ov-token-wrap');
  const prevWrap = gid('ov-pkg-preview-wrap');
  if (!selectedHelper || !selectedHelper.img_packages) {
    wrap.innerHTML = '';
    prevWrap.style.display = 'none';
    return;
  }

  const tokens = extractTokens(selectedHelper.img_packages);
  _ovTokenKeys = tokens;

  if (!tokens.length) {
    wrap.innerHTML = '';
    prevWrap.style.display = 'none';
    updateOvPreview();
    return;
  }

  const paramValues = selectedHelper.param_values || {};
  const params      = selectedHelper.params || [];

  wrap.innerHTML = '<div class="form-grid">' +
    tokens.map(key => {
      const values  = paramValues[key] || [];
      const isOpen  = values.length > 0 && values[values.length - 1] === '*';
      const opts    = isOpen ? values.slice(0, -1) : values;
      const param   = params.find(p => p.key === key);
      const defVal  = (param && param.default) || (opts[0] || '');
      const desc    = param ? param.desc : '';
      return '<div class="field">' +
        '<label class="field-label">' + escHtml(key) + '</label>' +
        (opts.length
          ? makeComboboxHtml('ovtok-' + key, defVal)
          : '<input class="field-input" id="ovtok-' + escHtml(key) + '" value="' + escHtml(defVal) + '">') +
        (desc ? '<span style="font-size:10px;color:var(--muted);">' + escHtml(desc) + '</span>' : '') +
      '</div>';
    }).join('') +
  '</div>';

  // Wire comboboxes and inputs to preview updates
  tokens.forEach(key => {
    const values  = paramValues[key] || [];
    const isOpen  = values.length > 0 && values[values.length - 1] === '*';
    const opts    = isOpen ? values.slice(0, -1) : values;
    if (opts.length) {
      initCombobox('ovtok-' + key, opts, isOpen, updateOvPreview);
    } else {
      const el = gid('ovtok-' + key);
      if (el) el.addEventListener('input', updateOvPreview);
    }
  });

  prevWrap.style.display = '';
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
  if (!selectedHelper || !selectedHelper.img_packages) return;
  let s = selectedHelper.img_packages;
  _ovTokenKeys.forEach(k => {
    const el = gid('ovtok-' + k);
    const v  = el ? (el.dataset.value !== undefined ? el.dataset.value : el.value) : '';
    s = s.replaceAll('{' + k + '}', v || ('{' + k + '}'));
  });
  gid('ov-pkg-preview').textContent = s;
}

/* ── Searchable combobox ─────────────────── */
function makeComboboxHtml(id, defVal) {
  return '<div class="combo-wrap" id="combo-' + escHtml(id) + '">' +
    '<input class="field-input combo-input" id="' + escHtml(id) + '" autocomplete="off"' +
      ' value="' + escHtml(defVal) + '" data-value="' + escHtml(defVal) + '">' +
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
    // Validate closed lists
    if (!isOpen && options.length && !options.includes(input.value)) {
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

/* ── Create overlay ──────────────────────── */
async function createOverlay() {
  const name = gid('ov-name').value || 'my-env.img';
  const size  = gid('ov-size').value || '8G';

  // Validate closed comboboxes
  for (const key of _ovTokenKeys) {
    const paramValues = (selectedHelper && selectedHelper.param_values) || {};
    const values = paramValues[key] || [];
    const isOpen = values.length > 0 && values[values.length - 1] === '*';
    const opts   = isOpen ? values.slice(0, -1) : values;
    if (!isOpen && opts.length) {
      const el = gid('ovtok-' + key);
      const v  = el ? el.value : '';
      if (!opts.includes(v)) {
        el && el.classList.add('combo-invalid');
        return;
      }
    }
  }

  // Resolve {KEY} tokens in img_packages
  const params = _collectParams();
  let basePkgs = selectedHelper?.img_packages || '';
  _ovTokenKeys.forEach(k => {
    const el = gid('ovtok-' + k);
    const v  = el ? (el.dataset.value !== undefined ? el.dataset.value : el.value) : '';
    basePkgs = basePkgs.replaceAll('{' + k + '}', v);
  });
  // Also resolve helper params ({KEY} in img_packages that come from #PARAM:)
  for (const [k, v] of Object.entries(params)) {
    basePkgs = basePkgs.replaceAll('{' + k + '}', v);
  }

  const extraPkgs = gid('ov-pkgs').value;
  const body = {
    name,
    size,
    base_packages:    basePkgs,
    extra_packages:   extraPkgs,
    post_install_cmd: selectedHelper?.post_install_cmd || '',
  };

  try {
    const r = await fetch('/api/overlay/create', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(body),
    });
    if (!r.ok) {
      const e = await r.json().catch(() => ({}));
      alert('Error: ' + (e.error || r.status));
      return;
    }
    const { id } = await r.json();
    openProgressWindow('Creating ' + name, id, async msg => {
      if (msg.ok && msg.path) {
        _setVal('cfg-overlay', msg.path);
        gid('cfg-overlay').dispatchEvent(new Event('input'));
        await loadOverlays();
        // collapse create form
        gid('ov-create-form').classList.remove('visible');
      }
    });
  } catch (e) {
    alert('Error: ' + e);
  }
}

/* ── Start helper job ────────────────────── */
async function startHelper() {
  if (!selectedHelper) return;
  const startBtn = gid('start-btn');
  const body = {
    cpus:    parseInt(gid('cfg-cpus').value)  || 4,
    mem:     gid('cfg-mem').value   || '16G',
    time:    gid('cfg-wall').value  || '08:00:00',
    gpu:     gid('cfg-gpu').value   || '',
    cwd:     gid('cfg-cwd').value   || '',
    overlay: gid('cfg-overlay').value || '',
    overlays: selectedModules.map(m => m.path),
    params:  _collectParams(),
  };

  const term = gid('start-term');
  term.classList.add('visible');
  term.innerHTML = '';
  if (startBtn) {
    startBtn.disabled = true;
    startBtn.innerHTML = iconSvg('play_arrow') + 'Starting…';
  }

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
      term.innerHTML = '<div class="tl-err">' + escHtml(msg) + '</div>';
      return;
    }

    const reader = r.body.getReader();
    const dec    = new TextDecoder();
    while (true) {
      const { done, value } = await reader.read();
      if (done) break;
      dec.decode(value).split('\n').forEach(line => {
        if (!line) return;
        const div = document.createElement('div');
        div.className = 'tl-info';
        div.textContent = line;
        term.appendChild(div);
        term.scrollTop = term.scrollHeight;
      });
    }
    await loadJobs();
    resetStartForm(selectedHelper);
    setTimeout(() => navigate('jobs'), 500);
  } catch (e) {
    term.innerHTML += '<div class="tl-err">Error: ' + escHtml(String(e)) + '</div>';
  } finally {
    if (startBtn) {
      startBtn.disabled = false;
      startBtn.innerHTML = iconSvg('play_arrow') + ' Start ' + escHtml(selectedHelper ? selectedHelper.name : '');
    }
  }
}

/* ── Module chips ────────────────────────── */
function renderModuleChips() {
  const list = gid('cfg-modules-list');
  if (!list) return;

  // Required overlays from the helper script — greyed out, not removable.
  const requiredNames = (selectedHelper?.required_overlays || '').trim().split(/\s+/).filter(Boolean);
  const requiredChips = requiredNames.map(name =>
    '<div style="display:flex;align-items:center;gap:6px;padding:4px 8px;' +
      'background:var(--surface2);border:1px solid var(--border);border-radius:var(--radius-sm);opacity:.5;">' +
      '<span style="font-family:var(--mono);font-size:11px;flex:1;">' + escHtml(name) + '</span>' +
      '<span style="font-size:10px;color:var(--muted);">default</span>' +
    '</div>'
  ).join('');

  const userChips = selectedModules.map((m, i) =>
    '<div style="display:flex;align-items:center;gap:6px;padding:4px 8px;' +
      'background:var(--surface2);border:1px solid var(--border);border-radius:var(--radius-sm);">' +
      '<span style="font-family:var(--mono);font-size:11px;flex:1;">' + escHtml(m.name) + '</span>' +
      '<button class="btn btn-ghost btn-sm" style="padding:0 4px;" onclick="removeModule(' + i + ')">' + iconSvg('close') + '</button>' +
    '</div>'
  ).join('');

  list.innerHTML = requiredChips + userChips;
}

function removeModule(i) {
  selectedModules.splice(i, 1);
  renderModuleChips();
}
