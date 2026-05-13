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
      '<tr><td colspan="3" style="padding:20px;text-align:center;color:var(--danger);">Error loading path.</td></tr>';
  }
  renderFileTree();
}

function renderFileBreadcrumb(path) {
  const parts = (path || '/').split('/').filter(Boolean);
  gid('bc').innerHTML =
    '<span class="bc-seg" onclick="navigateFiles(\'/\')">/ </span>' +
    parts.map((p, i) => {
      const full = '/' + parts.slice(0, i + 1).join('/');
      return '<span class="bc-sep">/</span>' +
        '<span class="bc-seg" onclick="navigateFiles(\'' + escHtml(full) + '\')">' + escHtml(p) + '</span>';
    }).join('');
}

function renderFileListing(entries, truncated) {
  const tbody = gid('file-tbody');
  const vis = showHiddenFiles ? entries : entries.filter(e => !e.name.startsWith('.'));
  if (!vis.length) {
    tbody.innerHTML = '<tr><td colspan="3" style="padding:20px;text-align:center;color:var(--muted);">Empty directory</td></tr>';
    return;
  }
  const truncRow = truncated
    ? '<tr><td colspan="3" style="padding:6px 12px;text-align:center;color:var(--warn,#e6a817);font-size:12px;">' +
      '⚠ Showing first 500 entries — use the path bar to navigate deeper.</td></tr>'
    : '';
  tbody.innerHTML = truncRow + vis.map(e => {
    if (e.is_dir) {
      return '<tr style="cursor:pointer;" onclick="navigateFiles(\'' + escHtml(e.path) + '\')">' +
        '<td>' + iconSvg('folder') + ' ' + escHtml(e.name) + '</td>' +
        '<td class="mono" style="color:var(--muted);">—</td>' +
        '<td style="color:var(--muted);">' + (e.modified_at ? fmtDate(e.modified_at) : '—') + '</td>' +
        '</tr>';
    }
    return '<tr>' +
      '<td>' + iconSvg('draft') + ' ' + escHtml(e.name) + '</td>' +
      '<td class="mono" style="color:var(--muted);">' + fmtSize(e.size) + '</td>' +
      '<td style="color:var(--muted);">' + (e.modified_at ? fmtDate(e.modified_at) : '—') + '</td>' +
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

function openFilePicker(targetId, mode, fallbackId) {
  fpTargetId = targetId;
  fpMode     = mode;
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

function _renderFpBreadcrumb(path) {
  const parts = (path || '/').split('/').filter(Boolean);
  gid('fp-bc').innerHTML =
    '<span class="bc-seg" onclick="fpNavigate(\'/\')">/ </span>' +
    parts.map((p, i) => {
      const full = '/' + parts.slice(0, i + 1).join('/');
      return '<span class="bc-sep">/</span>' +
        '<span class="bc-seg" onclick="fpNavigate(\'' + escHtml(full) + '\')">' + escHtml(p) + '</span>';
    }).join('');
}

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
    const files = fpMode !== 'dir' ? entries.filter(e => !e.is_dir && (showHiddenFiles || !e.name.startsWith('.'))) : [];
    gid('fp-body').innerHTML =
      dirs.map(d =>
        '<div class="modal-row" onclick="fpNavigate(\'' + escHtml(d.path) + '\')">' +
          '<span class="modal-row-icon">' + iconSvg('folder') + '</span>' +
          '<span class="modal-row-name">' + escHtml(d.name) + '</span>' +
          '<span class="modal-row-sel">' + iconSvg('chevron_right') + '</span>' +
        '</div>'
      ).join('') +
      files.map(f =>
        '<div class="modal-row" onclick="fpSelectFile(\'' + escHtml(f.path) + '\')">' +
          '<span class="modal-row-icon">' + iconSvg('draft') + '</span>' +
          '<span class="modal-row-name">' + escHtml(f.name) + '</span>' +
          '<span class="modal-row-size mono">' + fmtSize(f.size) + '</span>' +
          '<span class="modal-row-sel">' +
            '<button class="btn btn-sm btn-primary" onclick="event.stopPropagation();fpSelectFile(\'' + escHtml(f.path) + '\')">Select</button>' +
          '</span></div>'
      ).join('') ||
      '<div style="padding:20px;text-align:center;color:var(--muted);">Empty directory</div>';
  } catch(err) {
    if (err.name === 'AbortError') return;
    gid('fp-body').innerHTML =
      '<div style="padding:20px;text-align:center;color:var(--danger);">Error loading directory</div>';
  }
}

function fpSelectDir(path)  {
  _setVal(fpTargetId, path);
  const e = gid(fpTargetId);
  if (e) e.dispatchEvent(new Event('input'));
  closeModal('fp-modal');
}
function fpSelectFile(path) {
  _setVal(fpTargetId, path);
  const e = gid(fpTargetId);
  if (e) e.dispatchEvent(new Event('input'));
  closeModal('fp-modal');
}
function fpSelect()         { fpSelectDir(fpPath); }
