/* ── Overlays section ────────────────────── */
async function loadOverlays() {
  try {
    const r = await fetch('/api/overlays');
    allOverlays = (await r.json()) || [];
  } catch { allOverlays = []; }
  renderOverlays();
}

function renderOverlays(filter) {
  const q    = (filter || '').toLowerCase();
  const sqfs = allOverlays.filter(o => o.type === 'sqf');
  const rows = q ? sqfs.filter(o => o.name.toLowerCase().includes(q)) : sqfs;
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
      '<td><div class="td-actions">' +
        '<button class="btn btn-sm" onclick="pickOverlay(\'' + escHtml(o.path) + '\')">Use</button>' +
        '<button class="btn btn-sm btn-ghost btn-icon" onclick="copyStr(\'' + escHtml(o.path) + '\',this)" title="Copy path">' + iconSvg('content_copy') + '</button>' +
      '</div></td>' +
    '</tr>'
  ).join('');
}

function _addModuleOverlay(path) {
  const entry = allOverlays.find(o => o.path === path);
  if (entry && !selectedModules.find(m => m.path === path)) {
    selectedModules.push(entry);
    renderModuleChips();
  }
}

function pickOverlay(path) {
  _addModuleOverlay(path);
  navigate('start');
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
  const lower  = q.toLowerCase();
  const hits   = lower ? sqfOvs.filter(o => o.name.toLowerCase().includes(lower)) : sqfOvs;
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
