/* ── Files section ───────────────────────── */
let _navController  = null; // AbortController for the in-flight navigateFiles request
let showHiddenFiles = localStorage.getItem('showHiddenFiles') === 'true';
if (showHiddenFiles) { const cb = gid('show-hidden-files'); if (cb) cb.checked = true; }

async function navigateFiles(path) {
  if (_navController) _navController.abort();
  _navController = new AbortController();
  const url = path ? '/api/fs?path=' + encodeURIComponent(path) : '/api/fs';
  try {
    const r       = await fetch(url, { signal: _navController.signal });
    const entries = (await r.json()) || [];
    const truncated = r.headers.get('X-FS-Truncated') === 'true';
    // Resolve actual path from first entry's parent (API returns full paths)
    const resolved = entries.length > 0
      ? entries[0].path.replace(/\/[^/]+$/, '') || '/'
      : (path || '/');
    currentPath = resolved;
    renderFileBreadcrumb(resolved);
    _setVal('path-input', resolved);
    renderFileListing(entries, truncated);
  } catch(err) {
    if (err.name === 'AbortError') return;
    gid('file-tbody').innerHTML =
      '<tr><td colspan="3" class="td-error">Error loading path.</td></tr>';
  }
  renderFileTree();
}

function _renderBreadcrumb(containerId, path, navFn) {
  const parts = (path || '/').split('/').filter(Boolean);
  gid(containerId).innerHTML =
    '<span class="bc-seg" onclick="' + navFn + '(\'/\')">/</span>' +
    parts.map((p, i) => {
      const full = '/' + parts.slice(0, i + 1).join('/');
      const isLast = i === parts.length - 1;
      if (isLast) return '<span class="bc-cur">' + escHtml(p) + '</span>';
      return '<span class="bc-seg" onclick="' + navFn + '(\'' + escHtml(full) + '\')">' + escHtml(p) + '</span>' +
        '<span class="bc-sep">/</span>';
    }).join('');
}
function renderFileBreadcrumb(path) { _renderBreadcrumb('bc',    path, 'navigateFiles'); }

function renderFileListing(entries, truncated) {
  const tbody = gid('file-tbody');
  const vis = showHiddenFiles ? entries : entries.filter(e => !e.name.startsWith('.'));
  if (!vis.length) {
    tbody.innerHTML = '<tr><td colspan="3" class="td-empty">Empty directory</td></tr>';
    return;
  }
  const truncRow = truncated
    ? '<tr><td colspan="3" class="td-warn">' +
      iconSvg('warning') + ' Showing first 500 entries — use the path bar to navigate deeper.</td></tr>'
    : '';
  tbody.innerHTML = truncRow + vis.map(e => {
    if (e.is_dir) {
      return '<tr onclick="navigateFiles(\'' + escHtml(e.path) + '\')">' +
        '<td>' + iconSvg('folder') + ' ' + escHtml(e.name) + '</td>' +
        '<td class="mono td-muted">—</td>' +
        '<td class="td-muted">' + (e.modified_at ? fmtDate(e.modified_at) : '—') + '</td>' +
        '</tr>';
    }
    return '<tr>' +
      '<td>' + iconSvg('draft') + ' ' + escHtml(e.name) + '</td>' +
      '<td class="mono td-muted">' + fmtSize(e.size) + '</td>' +
      '<td class="td-muted">' + (e.modified_at ? fmtDate(e.modified_at) : '—') + '</td>' +
      '</tr>';
  }).join('');
}

function renderFileTree() {
  const bookmarks = [['/', '/ Root']];
  if (srvHome)    bookmarks.push([srvHome,    '$HOME']);
  if (srvScratch) bookmarks.push([srvScratch, '$SCRATCH']);
  gid('file-tree').innerHTML = bookmarks.map(([p, label]) =>
    '<div class="tree-item ' + (p === currentPath ? 'active' : '') + '"' +
      ' onclick="navigateFiles(\'' + escHtml(p) + '\')">' + escHtml(label) + '</div>'
  ).join('');
}

function setShowHidden(val, save) {
  showHiddenFiles = val;
  if (save) localStorage.setItem('showHiddenFiles', val);
  const cb1 = gid('show-hidden-files');
  const cb2 = gid('fp-show-hidden');
  if (cb1) cb1.checked = val;
  if (cb2) cb2.checked = val;
  navigateFiles(currentPath || '/');
  const fpModal = gid('fp-modal');
  if (fpModal && fpModal.classList.contains('open') && fpPath) fpNavigate(fpPath);
}

gid('path-input').addEventListener('keydown', e => {
  if (e.key === 'Enter') navigateFiles(e.target.value.trim() || '/');
});

/* ── File Picker Modal ───────────────────── */
function _parentDir(path) {
  const trimmed = (path || '').replace(/\/+$/, '');
  if (!trimmed || trimmed === '/') return '/';
  const parent = trimmed.replace(/\/[^/]+$/, '');
  return parent || '/';
}

function openFilePicker(targetId, mode, fallbackId, suffix) {
  if ((targetId === 'cfg-cwd' || targetId === 'cfg-overlay') &&
      ((typeof _startBusy !== 'undefined' && _startBusy) ||
       (typeof _formOpen !== 'undefined' && _formOpen) ||
       (typeof _opRunning !== 'undefined' && _opRunning))) {
    return;
  }
  fpTargetId = targetId;
  fpMode     = mode;
  fpSuffix   = suffix || '';
  const targetPath   = gid(targetId)?.value || '';
  const fallbackPath = (fallbackId && gid(fallbackId)?.value) || '';
  fpPath = mode === 'file' && targetPath
    ? _parentDir(targetPath)
    : (targetPath || fallbackPath || srvScratch || srvHome || '/');
  gid('fp-title').textContent = mode === 'dir' ? 'Select Directory' : 'Select File';
  const selDirBtn = gid('fp-select-dir-btn');
  if (selDirBtn) selDirBtn.style.display = mode === 'file' ? 'none' : '';
  fpNavigate(fpPath);
  gid('fp-modal').classList.add('open');
}

function _renderFpBreadcrumb(path) { _renderBreadcrumb('fp-bc', path, 'fpNavigate'); }

let _fpController = null; // AbortController for the in-flight fpNavigate request

async function fpNavigate(path) {
  if (_fpController) _fpController.abort();
  _fpController = new AbortController();
  fpPath = path;
  const url = path ? '/api/fs?path=' + encodeURIComponent(path) : '/api/fs';
  try {
    const r       = await fetch(url, { signal: _fpController.signal });
    const entries = (await r.json()) || [];
    // Resolve actual path from API response
    if (entries.length > 0) {
      fpPath = entries[0].path.replace(/\/[^/]+$/, '') || '/';
    }
    _renderFpBreadcrumb(fpPath);
    const dirs  = entries.filter(e =>  e.is_dir && (showHiddenFiles || !e.name.startsWith('.')));
    const files = fpMode !== 'dir' ? entries.filter(e => !e.is_dir && (showHiddenFiles || !e.name.startsWith('.')) && (!fpSuffix || e.name.endsWith(fpSuffix))) : [];
    gid('fp-body').innerHTML =
      dirs.map(d =>
        '<div class="modal-row" onclick="fpNavigate(\'' + escHtml(d.path) + '\')">' +
          '<span class="modal-row-icon">' + iconSvg('folder') + '</span>' +
          '<span class="modal-row-name">' + escHtml(d.name) + '</span>' +
          '<span class="modal-row-sel">' + iconSvg('chevron_right') + '</span>' +
        '</div>'
      ).join('') +
      files.map(f =>
        '<div class="modal-row" onclick="fpSelectDir(\'' + escHtml(f.path) + '\')">' +
          '<span class="modal-row-icon">' + iconSvg('draft') + '</span>' +
          '<span class="modal-row-name">' + escHtml(f.name) + '</span>' +
          '<span class="modal-row-size mono">' + fmtSize(f.size) + '</span>' +
          '<span class="modal-row-sel">' +
            '<button class="btn btn-sm btn-primary" onclick="event.stopPropagation();fpSelectDir(\'' + escHtml(f.path) + '\')">Select</button>' +
          '</span></div>'
      ).join('') ||
      '<div class="modal-empty">Empty directory</div>';
  } catch(err) {
    if (err.name === 'AbortError') return;
    gid('fp-body').innerHTML =
      '<div class="modal-error">Error loading directory</div>';
  }
}

function fpSelectDir(path) {
  _setVal(fpTargetId, path);
  const e = gid(fpTargetId);
  if (e) e.dispatchEvent(new Event('input'));
  closeModal('fp-modal');
}
function fpSelect() { fpSelectDir(fpPath); }
