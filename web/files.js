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
      '<tr><td colspan="4" class="td-error">Error loading path.</td></tr>';
  }
  renderFileTree();
  updateBookmarkStarBtn();
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

let _lastEntries   = [];
let _lastTruncated = false;

function renderFileListing(entries, truncated) {
  _lastEntries   = entries;
  _lastTruncated = truncated;
  const tbody = gid('file-tbody');
  const vis = showHiddenFiles ? entries : entries.filter(e => !e.name.startsWith('.'));
  if (!vis.length) {
    tbody.innerHTML = '<tr><td colspan="4" class="td-empty">Empty directory</td></tr>';
    return;
  }
  const truncRow = truncated
    ? '<tr><td colspan="4" class="td-warn">' +
      iconSvg('warning') + ' Showing first 500 entries — use the path bar to navigate deeper.</td></tr>'
    : '';
  tbody.innerHTML = truncRow + vis.map(e => {
    const dlTitle = e.is_dir ? 'Download as zip (symlinks inside the folder are skipped, not followed)' : 'Download';
    const dlBtn = '<button class="btn btn-sm btn-ghost" onclick="event.stopPropagation();downloadPath(\'' + escHtml(e.path) + '\')" title="' + escHtml(dlTitle) + '">' + iconSvg('download') + '</button>';
    const delBtn = '<button class="btn btn-sm btn-ghost" onclick="event.stopPropagation();deletePath(\'' + escHtml(e.path) + '\',\'' + escHtml(e.name) + '\')" title="Delete">' + iconSvg('delete') + '</button>';
    if (e.is_dir) {
      const isStarred = fileBookmarks.some(b => b.path === e.path);
      const bmBtn = '<button class="btn btn-sm btn-ghost' + (isStarred ? ' starred' : '') + '" onclick="event.stopPropagation();toggleBookmark(\'' + escHtml(e.path) + '\')" title="Bookmark">' +
        iconSvg('star', null, isStarred) + '</button>';
      return '<tr onclick="navigateFiles(\'' + escHtml(e.path) + '\')">' +
        '<td>' + iconSvg('folder', 'icon-folder', true) + ' ' + escHtml(e.name) + '</td>' +
        '<td class="mono td-muted">—</td>' +
        '<td class="td-muted">' + (e.modified_at ? fmtDate(e.modified_at) : '—') + '</td>' +
        '<td class="file-actions-cell">' + bmBtn + dlBtn + delBtn + '</td>' +
        '</tr>';
    }
    return '<tr onclick="viewPath(\'' + escHtml(e.path) + '\',' + (e.size || 0) + ')" title="Click to open">' +
      '<td>' + iconSvg('draft') + ' ' + escHtml(e.name) + '</td>' +
      '<td class="mono td-muted">' + fmtSize(e.size) + '</td>' +
      '<td class="td-muted">' + (e.modified_at ? fmtDate(e.modified_at) : '—') + '</td>' +
      '<td class="file-actions-cell">' + dlBtn + delBtn + '</td>' +
      '</tr>';
  }).join('');
}

// viewLargeFileWarnBytes is the size above which viewPath asks for
// confirmation before opening a file — clicking a listing row is easy to do
// by accident, and for a huge file that immediately kicks off a big
// transfer with no other feedback (unlike the explicit Download button).
const viewLargeFileWarnBytes = 1 << 30; // 1 GiB

// viewPath opens a file in a new tab with Content-Disposition: inline, so
// the browser renders it directly (plain text, a rendered HTML page, an
// image, ...) instead of downloading it. Distinct from downloadPath, which
// always forces a save-as download via the Download button. Types the
// browser can't render (e.g. binary formats) still fall back to its normal
// download prompt — that decision is made client-side by the browser based
// on the response's Content-Type, not by this app.
async function viewPath(path, size) {
  if (size > viewLargeFileWarnBytes) {
    const ok = await askConfirm('This file is ' + fmtSize(size) + '. Open it anyway?',
      { title: 'Open large file?', okLabel: 'Open' });
    if (!ok) return;
  }
  window.open('/api/fs/download?path=' + encodeURIComponent(path) + '&inline=1', '_blank');
}

function downloadPath(path) {
  window.location.href = '/api/fs/download?path=' + encodeURIComponent(path);
}

/* ── Delete ──────────────────────────────── */
let _deleteConfirmPath = null;

function deletePath(path, name) {
  _deleteConfirmPath = path;
  gid('delete-confirm-name').textContent = name;
  gid('delete-confirm-modal').classList.add('open');
}

function confirmDelete(mode) {
  gid('delete-confirm-modal').classList.remove('open');
  const path = _deleteConfirmPath;
  _deleteConfirmPath = null;
  if (mode !== 'yes' || !path) return;
  doDelete(path);
}
// Clicking the backdrop cancels, matching the "cancel" button (the generic
// .modal-backdrop click handler in app.js only hides the modal).
gid('delete-confirm-modal').addEventListener('click', e => {
  if (e.target.id === 'delete-confirm-modal') confirmDelete('cancel');
});

async function doDelete(path) {
  try {
    const r = await fetch('/api/fs?path=' + encodeURIComponent(path), { method: 'DELETE' });
    if (!r.ok) { alert('Could not delete: ' + (await r.text())); return; }
  } catch (e) {
    alert('Could not delete: ' + e);
    return;
  }
  // Drop a now-dangling bookmark for the deleted path, if any.
  if (fileBookmarks.some(b => b.path === path)) {
    try {
      fileBookmarks = (await (await fetch('/api/fs/bookmarks?path=' + encodeURIComponent(path), { method: 'DELETE' })).json()) || [];
    } catch { /* best-effort */ }
  }
  navigateFiles(currentPath);
}

function renderFileTree() {
  const permanent = [['/', '/ Root']];
  if (srvHome)    permanent.push([srvHome,    '$HOME']);
  if (srvScratch) permanent.push([srvScratch, '$SCRATCH']);
  let html = permanent.map(([p, label]) =>
    '<div class="tree-item ' + (p === currentPath ? 'active' : '') + '"' +
      ' onclick="navigateFiles(\'' + escHtml(p) + '\')">' + escHtml(label) + '</div>'
  ).join('');
  if (fileBookmarks.length) {
    html += '<div class="tree-sep"></div>' + fileBookmarks.map(b =>
      '<div class="tree-item tree-item-bm ' + (b.path === currentPath ? 'active' : '') + '"' +
        ' onclick="navigateFiles(\'' + escHtml(b.path) + '\')" title="' + escHtml(b.path) + '">' +
        '<span class="tree-item-label">' + escHtml(b.label || b.path) + '</span>' +
        '<button class="btn btn-ghost btn-sm tree-item-remove" onclick="event.stopPropagation();toggleBookmark(\'' + escHtml(b.path) + '\')" title="Remove bookmark">' + iconSvg('close') + '</button>' +
      '</div>'
    ).join('');
  }
  gid('file-tree').innerHTML = html;
}

/* ── Bookmarks (server-persisted, not localStorage) ─── */
let fileBookmarks   = [];
let _bookmarksLoaded = false;

async function loadFileBookmarks() {
  if (_bookmarksLoaded) return;
  _bookmarksLoaded = true;
  try {
    fileBookmarks = (await (await fetch('/api/fs/bookmarks')).json()) || [];
  } catch { fileBookmarks = []; }
  renderFileTree();
  updateBookmarkStarBtn();
}

async function toggleBookmark(path) {
  if (!path) return;
  const existing = fileBookmarks.some(b => b.path === path);
  try {
    const r = existing
      ? await fetch('/api/fs/bookmarks?path=' + encodeURIComponent(path), { method: 'DELETE' })
      : await fetch('/api/fs/bookmarks', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ path })
        });
    if (!r.ok) { alert('Could not update bookmark: ' + (await r.text())); return; }
    fileBookmarks = (await r.json()) || [];
  } catch (e) {
    alert('Could not update bookmark: ' + e);
    return;
  }
  renderFileTree();
  updateBookmarkStarBtn();
  renderFileListing(_lastEntries, _lastTruncated); // refresh per-row star icons too
}

function updateBookmarkStarBtn() {
  const btn = gid('bookmark-star-btn');
  if (!btn) return;
  const starred = fileBookmarks.some(b => b.path === currentPath);
  btn.innerHTML = iconSvg('star', null, starred);
  btn.classList.toggle('starred', starred);
  btn.title = starred ? 'Remove bookmark' : 'Bookmark this folder';
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

/* ── Upload ──────────────────────────────── */
let _pendingUpload      = null; // [{file, relPath}] awaiting conflict resolution / upload
let _uploadConflictResolve = null; // resolves the promise inside resolveConflictsSequentially
let _uploadBusy = false; // true from file selection through the end of the XHR (covers the conflict-resolution modal too)

// Locks/unlocks the upload buttons and shows/hides the progress bar. Locking
// starts as soon as files are picked (not just during the XHR) so a second
// upload can't be queued while the conflict-resolution modal is still being
// worked through — this dashboard only tracks one upload at a time, it
// doesn't support a concurrent multi-upload queue.
function _setUploadBusy(busy) {
  _uploadBusy = busy;
  const fileBtn = gid('upload-file-btn');
  const dirBtn  = gid('upload-dir-btn');
  if (fileBtn) fileBtn.disabled = busy;
  if (dirBtn)  dirBtn.disabled  = busy;
  const bar = gid('upload-progress-bar');
  if (bar) bar.hidden = !busy;
  if (busy) {
    const fill = gid('upload-progress-fill');
    if (fill) fill.style.width = '0%';
    const label = gid('upload-progress-label');
    if (label) label.textContent = 'Starting…';
  }
}

function _setUploadProgress(loaded, total) {
  const pct = total > 0 ? Math.min(100, Math.round((loaded / total) * 100)) : 0;
  const fill = gid('upload-progress-fill');
  if (fill) fill.style.width = pct + '%';
  const label = gid('upload-progress-label');
  if (label) label.textContent = fmtSize(loaded) + ' / ' + fmtSize(total) + ' (' + pct + '%)';
}

function handleUploadInputChange(inputEl, isDir) {
  if (_uploadBusy) return; // buttons are disabled while busy, but guard directly too
  const files = Array.from(inputEl.files || []);
  inputEl.value = ''; // allow re-selecting the same file(s) later
  if (!files.length) return;
  startUpload(files.map(f => ({ file: f, relPath: f.webkitRelativePath || f.name })));
}

// startUpload is the single entry point into the upload flow (conflict
// check → resolution → XHR), shared by the file/folder <input> pickers and
// the drag-and-drop zone below. items is [{file, relPath}].
function startUpload(items) {
  if (_uploadBusy || !items.length) return;
  _pendingUpload = items;
  _setUploadBusy(true);
  checkUploadConflicts();
}

/* ── Drag-and-drop upload ─────────────────── */
// Reads all entries out of a DirectoryReader — Chrome caps each readEntries()
// call at ~100 results, so it must be called repeatedly until it returns an
// empty batch, not just once.
function _readAllDirEntries(reader) {
  return new Promise((resolve, reject) => {
    let all = [];
    (function readBatch() {
      reader.readEntries(entries => {
        if (!entries.length) { resolve(all); return; }
        all = all.concat(entries);
        readBatch();
      }, reject);
    })();
  });
}

function _fileFromEntry(entry) {
  return new Promise((resolve, reject) => entry.file(resolve, reject));
}

// Recursively walks a FileSystemEntry (from DataTransferItem.webkitGetAsEntry),
// appending {file, relPath} pairs to out — same shape handleUploadInputChange
// builds from an <input> selection, so both feed into startUpload identically.
async function _walkDroppedEntry(entry, basePath, out) {
  if (entry.isFile) {
    const file = await _fileFromEntry(entry);
    out.push({ file, relPath: basePath + entry.name });
  } else if (entry.isDirectory) {
    const entries = await _readAllDirEntries(entry.createReader());
    for (const child of entries) {
      await _walkDroppedEntry(child, basePath + entry.name + '/', out);
    }
  }
}

const _dropZone = gid('file-list-wrap');
let _dragDepth = 0; // counts nested dragenter/dragleave so the highlight doesn't flicker over child rows

_dropZone.addEventListener('dragover', e => {
  if (_uploadBusy || !e.dataTransfer.types.includes('Files')) return;
  e.preventDefault(); // required to allow a drop
});

_dropZone.addEventListener('dragenter', e => {
  if (_uploadBusy || !e.dataTransfer.types.includes('Files')) return;
  e.preventDefault();
  _dragDepth++;
  _dropZone.classList.add('drag-over');
});

_dropZone.addEventListener('dragleave', () => {
  _dragDepth = Math.max(0, _dragDepth - 1);
  if (_dragDepth === 0) _dropZone.classList.remove('drag-over');
});

_dropZone.addEventListener('drop', async e => {
  e.preventDefault();
  _dragDepth = 0;
  _dropZone.classList.remove('drag-over');
  if (_uploadBusy) return;

  const dtItems = e.dataTransfer.items;
  if (!dtItems || !dtItems.length) return;
  const entries = Array.from(dtItems)
    .map(it => (it.webkitGetAsEntry ? it.webkitGetAsEntry() : null))
    .filter(Boolean);
  if (!entries.length) return;

  const items = [];
  for (const entry of entries) {
    await _walkDroppedEntry(entry, '', items);
  }
  startUpload(items);
});

async function checkUploadConflicts() {
  try {
    const r = await fetch('/api/fs/upload/check?dest=' + encodeURIComponent(currentPath), {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ paths: _pendingUpload.map(p => p.relPath) })
    });
    if (!r.ok) { alert('Could not check upload targets: ' + (await r.text())); _pendingUpload = null; _setUploadBusy(false); return; }
    const { conflicts } = await r.json();
    const skip = await resolveConflictsSequentially(conflicts || []);
    if (skip === null) { _pendingUpload = null; _setUploadBusy(false); return; } // cancelled
    doUpload(skip);
  } catch (e) {
    alert('Could not check upload targets: ' + e);
    _pendingUpload = null;
    _setUploadBusy(false);
  }
}

// Walks conflicting relative paths one at a time via #upload-conflict-modal,
// honoring Yes / Yes to All / No / No to All / Cancel (classic Explorer/Finder
// style). Returns a Set of relPaths to skip, or null if the user cancelled.
function resolveConflictsSequentially(conflicts) {
  if (!conflicts.length) return Promise.resolve(new Set());
  return new Promise(resolve => {
    const skip = new Set();
    let forceMode = null; // 'yes' | 'no' once a "to All" choice is made
    let i = 0;
    function next() {
      if (i >= conflicts.length) { resolve(skip); return; }
      const relPath = conflicts[i++];
      if (forceMode === 'yes') { next(); return; }
      if (forceMode === 'no')  { skip.add(relPath); next(); return; }
      const remaining = conflicts.length - i + 1;
      gid('upload-conflict-path').textContent = relPath;
      gid('upload-conflict-yes-all-btn').style.display = remaining > 1 ? '' : 'none';
      gid('upload-conflict-no-all-btn').style.display  = remaining > 1 ? '' : 'none';
      gid('upload-conflict-modal').classList.add('open');
      _uploadConflictResolve = mode => {
        gid('upload-conflict-modal').classList.remove('open');
        if (mode === 'cancel')  { resolve(null); return; }
        if (mode === 'yes-all') { forceMode = 'yes'; next(); return; }
        if (mode === 'no-all')  { forceMode = 'no'; skip.add(relPath); next(); return; }
        if (mode === 'no') skip.add(relPath);
        next();
      };
    }
    next();
  });
}

function resolveUploadConflict(mode) {
  if (_uploadConflictResolve) {
    const resolve = _uploadConflictResolve;
    _uploadConflictResolve = null;
    resolve(mode);
  }
}
// Clicking the backdrop cancels the conflict flow (the generic .modal-backdrop
// click handler in app.js only hides the modal — it doesn't resolve our promise).
gid('upload-conflict-modal').addEventListener('click', e => {
  if (e.target.id === 'upload-conflict-modal') resolveUploadConflict('cancel');
});

function doUpload(skipSet) {
  const items = _pendingUpload.filter(p => !skipSet.has(p.relPath));
  _pendingUpload = null;
  if (!items.length) { _setUploadBusy(false); return; }
  // relPath goes in the field *name*, not the filename: the server relies on
  // this because Go's mime/multipart strips any "/" from Part.FileName() but
  // leaves Part.FormName() untouched (see api_files_upload.go).
  const fd = new FormData();
  items.forEach(p => fd.append(p.relPath, p.file, p.file.name));
  const xhr = new XMLHttpRequest();
  xhr.open('POST', '/api/fs/upload?dest=' + encodeURIComponent(currentPath) + '&overwrite=true');
  xhr.upload.onprogress = e => {
    if (e.lengthComputable) _setUploadProgress(e.loaded, e.total);
  };
  xhr.onload = () => {
    _setUploadBusy(false);
    if (xhr.status >= 200 && xhr.status < 300) {
      navigateFiles(currentPath);
    } else {
      alert('Upload failed: ' + xhr.responseText);
    }
  };
  xhr.onerror = () => { _setUploadBusy(false); alert('Upload failed.'); };
  xhr.send(fd);
}

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
          '<span class="modal-row-icon">' + iconSvg('folder', 'icon-folder', true) + '</span>' +
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
