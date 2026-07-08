/* ── Files section ───────────────────────── */
let _navController  = null; // AbortController for the in-flight navigateFiles request
let showHiddenFiles = localStorage.getItem('showHiddenFiles') === 'true';
if (showHiddenFiles) { const cb = gid('show-hidden-files'); if (cb) cb.checked = true; }

async function navigateFiles(path) {
  if (_navController) _navController.abort();
  _navController = new AbortController();
  const url = path ? '/api/fs?path=' + encodeURIComponent(path) : '/api/fs';
  // listing a big directory on a network filesystem can take a while —
  // show feedback instead of a frozen listing (only when it's actually slow)
  const slowTimer = setTimeout(() => {
    gid('file-tbody').innerHTML = '<tr><td colspan="4" class="td-empty">Loading…</td></tr>';
  }, 300);
  try {
    const r = await fetch(url, { signal: _navController.signal });
    if (!r.ok) {
      clearTimeout(slowTimer);
      const msg = (await r.text()).trim();
      gid('file-tbody').innerHTML =
        '<tr><td colspan="4" class="td-error">' + escHtml(msg || 'Error loading path.') + '</td></tr>';
      renderFileTree();
      updateBookmarkStarBtn();
      return;
    }
    const entries = (await r.json()) || [];
    clearTimeout(slowTimer);
    const truncated = r.headers.get('X-FS-Truncated') === 'true';
    // Resolve actual path from first entry's parent (API returns full paths)
    const resolved = entries.length > 0
      ? entries[0].path.replace(/\/[^/]+$/, '') || '/'
      : (path || '/');
    const pathChanged = resolved !== currentPath;
    if (pathChanged) {
      _selected.clear();
      _selAnchor = -1;
      _fileFilter = ''; // the filter applies to one directory; entering another starts fresh
      _setVal('file-search', '');
    }
    currentPath = resolved;
    renderFileBreadcrumb(resolved);
    _setVal('path-input', resolved);
    // real directory changes push a history entry so Back/Forward walk the
    // directory trail; refreshes just canonicalize the current one
    setHash('#files' + encodeURI(resolved), pathChanged);
    renderFileListing(entries, truncated);
  } catch(err) {
    clearTimeout(slowTimer);
    if (err.name === 'AbortError') return;
    gid('file-tbody').innerHTML =
      '<tr><td colspan="4" class="td-error">Error loading path.</td></tr>';
  }
  renderFileTree();
  updateBookmarkStarBtn();
}

// Segment clicks navigate and stop propagation — a click anywhere else in
// the crumbs swaps in the path input (see _setupPathBox).
function _renderBreadcrumb(containerId, path, navFn) {
  const parts = (path || '/').split('/').filter(Boolean);
  const el = gid(containerId);
  el.innerHTML =
    '<span class="bc-seg" onclick="event.stopPropagation();' + navFn + '(\'/\')">/</span>' +
    parts.map((p, i) => {
      const full = '/' + parts.slice(0, i + 1).join('/');
      const isLast = i === parts.length - 1;
      if (isLast) return '<span class="bc-cur">' + escHtml(p) + '</span>';
      return '<span class="bc-seg" onclick="event.stopPropagation();' + navFn + '(\'' + escHtml(full) + '\')">' + escHtml(p) + '</span>' +
        '<span class="bc-sep">/</span>';
    }).join('');
  el.scrollLeft = el.scrollWidth; // keep the tail of a long path visible
}
function renderFileBreadcrumb(path) { _renderBreadcrumb('bc',    path, 'navigateFiles'); }

/* ── Path box (breadcrumb ⇄ input) ───────── */
// Windows-Explorer-style path control: a clickable breadcrumb; clicking
// anywhere that isn't a segment link swaps in a text input prefilled with
// the current path. Enter navigates, Escape or leaving the field switches
// back to the breadcrumb.
function _setupPathBox(crumbsId, inputId, navFn) {
  const crumbs = gid(crumbsId), input = gid(inputId);
  crumbs.addEventListener('click', () => {
    crumbs.hidden = true;
    input.hidden  = false;
    input.focus();
    input.setSelectionRange(input.value.length, input.value.length);
  });
  const exit = () => { input.hidden = true; crumbs.hidden = false; };
  input.addEventListener('keydown', e => {
    if (e.key === 'Enter') { exit(); navFn(input.value.trim() || '/'); }
    else if (e.key === 'Escape') { e.stopPropagation(); exit(); }
  });
  input.addEventListener('blur', () => { if (!input.hidden) exit(); });
}
_setupPathBox('bc',    'path-input',    p => navigateFiles(p));
_setupPathBox('fp-bc', 'fp-path-input', p => fpNavigate(p));

let _lastEntries   = [];
let _lastTruncated = false;
let _visEntries    = []; // entries as currently rendered (filtered + sorted); row handlers index into this

/* ── Name filter (search bar) ────────────── */
let _fileFilter = ''; // current search-bar query; cleared when the directory changes

// filterFiles re-renders the listing showing only entries whose name
// contains every whitespace-separated term of q.
function filterFiles(q) {
  _fileFilter = q;
  renderFileListing(_lastEntries, _lastTruncated);
}

// Escape clears the filter without bubbling up to the listing-wide Escape
// handler (which clears the selection).
gid('file-search').addEventListener('keydown', e => {
  if (e.key === 'Escape' && e.target.value) {
    e.stopPropagation();
    e.target.value = '';
    filterFiles('');
  }
});

/* ── Sorting ─────────────────────────────── */
let fileSortKey = localStorage.getItem('fileSortKey') || 'name';
let fileSortAsc = localStorage.getItem('fileSortAsc') !== 'false';

// setFileSort sorts by the given column, or flips the direction when it's
// already the active column.
function setFileSort(key) {
  if (fileSortKey === key) fileSortAsc = !fileSortAsc;
  else { fileSortKey = key; fileSortAsc = true; }
  localStorage.setItem('fileSortKey', fileSortKey);
  localStorage.setItem('fileSortAsc', fileSortAsc);
  renderFileListing(_lastEntries, _lastTruncated);
}

function _updateSortHeaders() {
  const labels = { name: 'Name', size: 'Size', modified: 'Modified' };
  for (const k in labels) {
    const th = gid('th-' + k);
    if (th) th.textContent = labels[k] + (fileSortKey === k ? (fileSortAsc ? ' ↑' : ' ↓') : '');
  }
}
_updateSortHeaders(); // reflect the persisted sort in the headers before the first render

function _sortEntries(entries) {
  const dir = fileSortAsc ? 1 : -1;
  const byName = (a, b) => a.name.localeCompare(b.name, undefined, { numeric: true, sensitivity: 'base' });
  return entries.slice().sort((a, b) => {
    if (a.is_dir !== b.is_dir) return a.is_dir ? -1 : 1; // folders always first
    if (fileSortKey === 'size' && a.is_dir) return byName(a, b); // dirs have no size; keep them name-ordered
    let cmp = 0;
    if (fileSortKey === 'size')          cmp = (a.size || 0) - (b.size || 0);
    else if (fileSortKey === 'modified') cmp = new Date(a.modified_at || 0) - new Date(b.modified_at || 0);
    else                                 cmp = byName(a, b);
    return dir * cmp || byName(a, b); // ties break by name A→Z
  });
}

/* ── Selection ───────────────────────────── */
let _selected  = new Set(); // paths of selected rows
let _selAnchor = -1;        // _visEntries index of the last plainly-clicked row (shift-range anchor)

function _selectedEntries() { return _visEntries.filter(e => _selected.has(e.path)); }

// _applySelection re-applies .selected classes in place (no re-render, so
// double-click still works).
function _applySelection() {
  document.querySelectorAll('#file-tbody tr[data-idx]').forEach(tr => {
    const e = _visEntries[+tr.dataset.idx];
    tr.classList.toggle('selected', !!e && _selected.has(e.path));
  });
}

// onRowClick selects row i: plain click selects just that row, Ctrl/Cmd
// toggles it, Shift extends a range from the anchor.
function onRowClick(ev, i) {
  closeRowMenu();
  const e = _visEntries[i];
  if (!e) return;
  if (ev.shiftKey && _selAnchor >= 0) {
    if (!ev.ctrlKey && !ev.metaKey) _selected.clear();
    const [lo, hi] = _selAnchor < i ? [_selAnchor, i] : [i, _selAnchor];
    for (let j = lo; j <= hi; j++) _selected.add(_visEntries[j].path);
  } else if (ev.ctrlKey || ev.metaKey) {
    if (_selected.has(e.path)) _selected.delete(e.path); else _selected.add(e.path);
    _selAnchor = i;
  } else {
    _selected.clear();
    _selected.add(e.path);
    _selAnchor = i;
  }
  _applySelection();
}

function clearFileSelection() {
  if (!_selected.size) return;
  _selected.clear();
  _selAnchor = -1;
  _applySelection();
}

// openEntry navigates into a directory or views a file (double-click on a
// row, or single click on its name). Debounced so a double-click on the
// name can't open twice.
let _lastOpen = { i: -1, t: 0 };
function openEntry(i) {
  const now = Date.now();
  if (_lastOpen.i === i && now - _lastOpen.t < 500) return;
  _lastOpen = { i, t: now };
  const e = _visEntries[i];
  if (!e) return;
  if (e.is_dir) navigateFiles(e.path);
  else viewPath(e.path, e.size || 0);
}

/* ── Row "More" / context menu ───────────── */
// A position:fixed dropdown appended to <body> so it isn't clipped by the
// scrollable table; closed on outside click, Escape, or scroll. Shared by
// the ⋮ button and right-click — both act on the current selection.
let _rowMenuEl = null;

function closeRowMenu() {
  if (_rowMenuEl) { _rowMenuEl.remove(); _rowMenuEl = null; }
}

// _menuHtml builds the dropdown items: the full per-item menu for a single
// selection, bulk actions only for a multi-selection.
function _menuHtml(sel) {
  const mi = (icon, label, fn, cls) =>
    '<div class="row-menu-item' + (cls || '') + '" onclick="closeRowMenu();' + fn + '">' +
      iconSvg(icon) + ' ' + escHtml(label) + '</div>';
  const sep = '<div class="row-menu-sep"></div>';
  const clip =
    mi('content_cut', 'Cut', 'ctxCut()') +
    mi('content_copy', 'Copy', 'ctxCopy()') +
    (_clipboard ? mi('content_paste', 'Paste (' + _clipboard.entries.length + ')', 'ctxPaste()') : '');
  const moveCopy =
    mi('folder', 'Move to…', 'ctxMoveTo()') +
    mi('content_copy', 'Copy to…', 'ctxCopyTo()');
  const del = sep + mi('delete', 'Delete', 'ctxDelete()', ' row-menu-danger');
  if (sel.length === 1) {
    const e = sel[0];
    return mi('download', e.is_dir ? 'Download as zip' : 'Download', 'ctxDownload()') + sep + clip + sep +
      mi('edit', 'Rename', 'ctxRename()') + moveCopy + del;
  }
  return '<div class="row-menu-hdr">' + sel.length + ' selected</div>' +
    mi('download', 'Download', 'ctxDownload()') + sep + clip + sep + moveCopy + del;
}

// _openMenuAt shows a dropdown with the given items at pos ({top,left} or
// {top,right}), clamped to the viewport.
function _openMenuAt(pos, html, anchorPath) {
  const menu = document.createElement('div');
  menu.className = 'row-menu';
  menu.dataset.path = anchorPath || '';
  menu.innerHTML = html;
  document.body.appendChild(menu); // append first so its size is measurable
  const r = menu.getBoundingClientRect();
  menu.style.top = Math.max(4, Math.min(pos.top, window.innerHeight - r.height - 8)) + 'px';
  if (pos.left != null) menu.style.left = Math.max(4, Math.min(pos.left, window.innerWidth - r.width - 8)) + 'px';
  else menu.style.right = pos.right + 'px';
  _rowMenuEl = menu;
}

// _openRowMenu shows the menu for the current selection.
function _openRowMenu(pos, anchorPath) {
  const sel = _selectedEntries();
  if (!sel.length) return;
  _openMenuAt(pos, _menuHtml(sel), anchorPath);
}

// _menuTargetRow collapses the selection to row i unless the row is
// already part of it.
function _menuTargetRow(i) {
  const e = _visEntries[i];
  if (!e) return null;
  if (!_selected.has(e.path)) {
    _selected  = new Set([e.path]);
    _selAnchor = i;
    _applySelection();
  }
  return e;
}

async function toggleRowMenu(ev, i) {
  ev.stopPropagation();
  const e = _visEntries[i];
  const reopening = _rowMenuEl && e && _rowMenuEl.dataset.path === e.path;
  closeRowMenu();
  if (reopening || !_menuTargetRow(i)) return; // clicking the same row's button again just closes it
  const rect = ev.currentTarget.getBoundingClientRect(); // before await: currentTarget nulls after the handler returns
  await _validateClipboard();
  _openRowMenu({ top: rect.bottom + 4, right: window.innerWidth - rect.right }, e.path);
}

async function onRowContext(ev, i) {
  ev.preventDefault();
  ev.stopPropagation();
  closeRowMenu();
  const e = _menuTargetRow(i);
  if (!e) return;
  const pos = { top: ev.clientY, left: ev.clientX };
  await _validateClipboard();
  _openRowMenu(pos, e.path);
}

/* ── Context-menu actions (operate on the selection) ── */
// downloadEntry triggers a download for an entry. Directories are checked
// against the server's zip size limits first (check=1), so a too-large
// folder surfaces as an error pill — the <a download> click of the real
// download gives no error feedback.
async function downloadEntry(e) {
  if (e.is_dir) {
    try {
      const r = await fetch('/api/fs/download?path=' + encodeURIComponent(e.path) + '&check=1');
      if (!r.ok) { showProgressError('Download ' + e.name, await r.text()); return; }
    } catch (err) {
      showProgressError('Download ' + e.name, String(err));
      return;
    }
  }
  downloadPath(e.path);
}

function ctxDownload() {
  // Staggered so the browser registers each as its own download.
  _selectedEntries().forEach((e, i) => setTimeout(() => downloadEntry(e), i * 300));
}

function ctxRename() {
  const sel = _selectedEntries();
  if (sel.length === 1) renamePath(sel[0].path, sel[0].name);
}

function _selLabel(sel) { return sel.length === 1 ? '"' + sel[0].name + '"' : sel.length + ' items'; }

function ctxMoveTo() {
  const sel = _selectedEntries();
  if (!sel.length) return;
  openFilePickerForCallback(currentPath, 'Move ' + _selLabel(sel) + ' to…',
    (destDir, newName) => sel.forEach(e =>
      doMove(e.path, e.name, destDir, newName !== e.name ? newName : undefined)),
    sel.length === 1 ? { name: sel[0].name } : null);
}

function ctxCopyTo() {
  const sel = _selectedEntries();
  if (!sel.length) return;
  openFilePickerForCallback(currentPath, 'Copy ' + _selLabel(sel) + ' to…',
    (destDir, newName) => sel.forEach(e =>
      doCopy(e.path, e.name, destDir, newName !== e.name ? newName : undefined)),
    sel.length === 1 ? { name: sel[0].name } : null);
}

/* ── Clipboard (cut/copy/paste) ──────────── */
// An app-internal clipboard, not the OS one. Cut marks the selection to be
// moved on paste (rows dim until pasted or Escape); copy pastes as copies
// and survives multiple pastes.
let _clipboard = null; // { mode: 'cut'|'copy', entries: [{path, name}] }

function _renderCutMarks() {
  const cutSet = _clipboard && _clipboard.mode === 'cut'
    ? new Set(_clipboard.entries.map(e => e.path)) : null;
  document.querySelectorAll('#file-tbody tr[data-idx]').forEach(tr => {
    const e = _visEntries[+tr.dataset.idx];
    tr.classList.toggle('cut', !!(cutSet && e && cutSet.has(e.path)));
  });
}

function ctxCut()  { _setClipboard('cut'); }
function ctxCopy() { _setClipboard('copy'); }
// _setClipboard also mirrors the selected paths onto the OS clipboard as
// text (best-effort). The paste handler compares against it: if the user
// copies something else afterwards, the mismatch marks the file clipboard
// stale — matching how Explorer/Finder supersede a file copy with a text
// copy. Bonus: the paths can be pasted into a terminal.
function _setClipboard(mode) {
  const sel = _selectedEntries();
  if (!sel.length) return;
  _clipboard = { mode, entries: sel.map(e => ({ path: e.path, name: e.name })), osText: null };
  if (navigator.clipboard && navigator.clipboard.writeText) {
    const text = sel.map(e => e.path).join('\n');
    navigator.clipboard.writeText(text)
      .then(() => { if (_clipboard) _clipboard.osText = text; })
      .catch(() => { /* keep internal-only behavior */ });
  }
  _renderCutMarks();
}

// _validateClipboard drops a stale internal clipboard before a menu shows
// Paste: stale when the OS clipboard no longer holds what cut/copy
// mirrored there (the user copied something else since). Reading the OS
// clipboard needs a permission the browser may deny — then the internal
// clipboard is trusted as-is.
async function _validateClipboard() {
  if (!_clipboard || _clipboard.osText == null) return;
  if (!(navigator.clipboard && navigator.clipboard.readText)) return;
  try {
    const txt = await navigator.clipboard.readText();
    if (txt !== _clipboard.osText) { _clipboard = null; _renderCutMarks(); }
  } catch { /* permission denied — keep internal-only behavior */ }
}

// ctxPaste moves (cut) or copies (copy) the clipboard entries into the
// current directory. Cut entries already in this directory are skipped;
// a cut clipboard is one-shot, a copy clipboard can be pasted again.
function ctxPaste() {
  if (!_clipboard) return;
  const { mode, entries } = _clipboard;
  for (const e of entries) {
    if (mode === 'cut') {
      if (_parentDir(e.path) === currentPath) continue; // already here
      doMove(e.path, e.name, currentPath);
    } else {
      doCopy(e.path, e.name, currentPath);
    }
  }
  if (mode === 'cut') { _clipboard = null; _renderCutMarks(); }
}

function ctxDelete() {
  const sel = _selectedEntries();
  if (!sel.length) return;
  deletePaths(sel.map(e => e.path), sel.length === 1 ? sel[0].name : sel.length + ' items');
}

document.addEventListener('click', e => {
  if (_rowMenuEl && !_rowMenuEl.contains(e.target)) closeRowMenu();
});
document.addEventListener('keydown', e => {
  if (e.key === 'Escape') {
    closeRowMenu();
    clearFileSelection();
    if (_clipboard) { _clipboard = null; _renderCutMarks(); }
  }
});
gid('file-list-wrap').addEventListener('scroll', closeRowMenu);
// Clicking the blank area below the rows clears the selection.
gid('file-list-wrap').addEventListener('click', e => {
  if (!e.target.closest('tr')) clearFileSelection();
});

// Right-clicking the blank area shows a menu for the directory itself:
// Paste (when the clipboard is non-empty), Copy path, Select all, Refresh.
gid('file-list-wrap').addEventListener('contextmenu', async e => {
  if (e.target.closest('tr[data-idx]')) return; // rows have their own menu
  e.preventDefault();
  closeRowMenu();
  const pos = { top: e.clientY, left: e.clientX };
  await _validateClipboard();
  const mi = (icon, label, fn) =>
    '<div class="row-menu-item" onclick="closeRowMenu();' + fn + '">' + iconSvg(icon) + ' ' + label + '</div>';
  _openMenuAt(pos,
    (_clipboard ? mi('content_paste', 'Paste (' + _clipboard.entries.length + ')', 'ctxPaste()') : '') +
    mi('content_copy', 'Copy path', 'ctxCopyPath()') +
    mi('check', 'Select all', 'ctxSelectAll()') +
    mi('refresh', 'Refresh', 'navigateFiles(currentPath)'));
});

// ctxCopyPath copies the current directory's path to the OS clipboard.
// This also supersedes any internal file clipboard, like any text copy.
function ctxCopyPath() {
  if (navigator.clipboard && navigator.clipboard.writeText) {
    navigator.clipboard.writeText(currentPath).catch(() => {});
  }
}

function ctxSelectAll() {
  _selected = new Set(_visEntries.map(en => en.path));
  _applySelection();
}

// Listing shortcuts (Ctrl on Linux/Windows, Cmd on Mac): A select all,
// X cut, C copy, V paste, Delete delete selection. Ignored while a modal
// is open or focus is in a form control.
document.addEventListener('keydown', e => {
  if (!gid('sec-files').classList.contains('active')) return;
  const t = e.target;
  if (t && (t.tagName === 'INPUT' || t.tagName === 'TEXTAREA' || t.tagName === 'SELECT' || t.isContentEditable)) return;
  if (document.querySelector('.modal-backdrop.open')) return;
  const mod = (e.ctrlKey || e.metaKey) && !e.shiftKey && !e.altKey;
  const key = e.key.toLowerCase();
  if (mod && key === 'a') {
    e.preventDefault();
    ctxSelectAll();
  } else if (mod && (key === 'x' || key === 'c') && _selected.size) {
    // a live text selection (e.g. in a log pill) wins over the file rows
    if (String(window.getSelection())) return;
    e.preventDefault();
    if (key === 'x') ctxCut(); else ctxCopy();
  } else if (e.key === 'Delete' && _selected.size) {
    e.preventDefault();
    ctxDelete();
  }
});

// Ctrl/Cmd+V arrives as a native paste event, whose clipboardData is
// readable without any permission. If the OS clipboard no longer holds
// what cut/copy mirrored there (see _setClipboard), the user copied
// something else in between — the file clipboard is stale, so drop it
// instead of pasting.
document.addEventListener('paste', e => {
  if (!gid('sec-files').classList.contains('active')) return;
  const t = e.target;
  if (t && (t.tagName === 'INPUT' || t.tagName === 'TEXTAREA' || t.tagName === 'SELECT' || t.isContentEditable)) return;
  if (document.querySelector('.modal-backdrop.open')) return;
  if (!_clipboard) return;
  if (_clipboard.osText != null) {
    const txt = e.clipboardData ? e.clipboardData.getData('text/plain') : '';
    if (txt !== _clipboard.osText) {
      _clipboard = null;
      _renderCutMarks();
      return;
    }
  }
  e.preventDefault();
  ctxPaste();
});

/* ── Drag (rubber-band) selection ────────── */
// Press-and-drag over the listing draws a marquee; rows it crosses become
// the selection (Ctrl/Cmd-drag adds to the existing selection). Becomes a
// drag after a 4px threshold — below that it stays a plain click.
const _mq = { didDrag: false, startX: 0, startY: 0, base: null, rows: [], el: null, lastX: 0, lastY: 0, timer: null };

gid('file-list-wrap').addEventListener('mousedown', e => {
  if (e.button !== 0) return;
  if (e.target.closest('button') || e.target.closest('.file-name') || e.target.closest('thead')) return;
  const wrap = gid('file-list-wrap');
  const wr = wrap.getBoundingClientRect();
  // A press on the wrap's own scrollbars must stay a scrollbar drag.
  if (e.clientX - wr.left > wrap.clientWidth || e.clientY - wr.top > wrap.clientHeight) return;
  _mq.didDrag = false;
  _mq.startX  = e.clientX - wr.left + wrap.scrollLeft;
  _mq.startY  = e.clientY - wr.top  + wrap.scrollTop;
  _mq.lastX   = e.clientX;
  _mq.lastY   = e.clientY;
  _mq.base    = (e.ctrlKey || e.metaKey) ? new Set(_selected) : new Set();
  _mq.rows    = Array.from(document.querySelectorAll('#file-tbody tr[data-idx]')).map(tr => {
    const r = tr.getBoundingClientRect();
    const top = r.top - wr.top + wrap.scrollTop;
    const entry = _visEntries[+tr.dataset.idx];
    return { path: entry && entry.path, top, bottom: top + r.height };
  });
  document.addEventListener('mousemove', _mqMove);
  document.addEventListener('mouseup', _mqUp);
  e.preventDefault(); // no native text/element drag while marqueeing (the click still fires)
});

function _mqMove(e) {
  _mq.lastX = e.clientX;
  _mq.lastY = e.clientY;
  if (!_mq.didDrag) {
    const wrap = gid('file-list-wrap');
    const wr = wrap.getBoundingClientRect();
    const x = e.clientX - wr.left + wrap.scrollLeft;
    const y = e.clientY - wr.top  + wrap.scrollTop;
    if (Math.abs(x - _mq.startX) < 4 && Math.abs(y - _mq.startY) < 4) return;
    _mq.didDrag = true;
    closeRowMenu();
    _mq.el = document.createElement('div');
    _mq.el.className = 'marquee';
    wrap.appendChild(_mq.el);
    // keeps edge auto-scroll going while the mouse sits still past the edge
    _mq.timer = setInterval(_mqUpdate, 50);
  }
  _mqUpdate();
}

function _mqUpdate() {
  const wrap = gid('file-list-wrap');
  const wr = wrap.getBoundingClientRect();
  if (_mq.lastY < wr.top)              wrap.scrollTop -= Math.min(30, (wr.top - _mq.lastY) / 2);
  else if (_mq.lastY > wr.bottom)      wrap.scrollTop += Math.min(30, (_mq.lastY - wr.bottom) / 2);
  // clamp the moving corner to the visible box, in content coordinates
  const x = Math.min(Math.max(_mq.lastX - wr.left, 0), wrap.clientWidth)  + wrap.scrollLeft;
  const y = Math.min(Math.max(_mq.lastY - wr.top,  0), wrap.clientHeight) + wrap.scrollTop;
  const left = Math.min(x, _mq.startX), top = Math.min(y, _mq.startY);
  const w = Math.abs(x - _mq.startX),   h = Math.abs(y - _mq.startY);
  Object.assign(_mq.el.style, { left: left + 'px', top: top + 'px', width: w + 'px', height: h + 'px' });
  // rows span the full width, so vertical overlap is the whole test
  _selected = new Set(_mq.base);
  for (const r of _mq.rows) {
    if (r.path && r.bottom > top && r.top < top + h) _selected.add(r.path);
  }
  _applySelection();
}

function _mqUp() {
  document.removeEventListener('mousemove', _mqMove);
  document.removeEventListener('mouseup', _mqUp);
  if (_mq.timer) { clearInterval(_mq.timer); _mq.timer = null; }
  if (_mq.el)    { _mq.el.remove(); _mq.el = null; }
  // _mq.didDrag stays set until the trailing click is swallowed below
}

// Swallow the click the browser fires after a real drag, so it can't
// rewrite the selection the marquee just made (capture phase, so it runs
// before the row/blank-area handlers).
document.addEventListener('click', e => {
  if (_mq.didDrag) { _mq.didDrag = false; e.stopPropagation(); }
}, true);

function rowDownload(ev, i) { ev.stopPropagation(); const e = _visEntries[i]; if (e) downloadEntry(e); }
function rowBookmark(ev, i) { ev.stopPropagation(); const e = _visEntries[i]; if (e) toggleBookmark(e.path); }

function renderFileListing(entries, truncated) {
  _lastEntries   = entries;
  _lastTruncated = truncated;
  _updateSortHeaders();
  closeRowMenu();
  const tbody = gid('file-tbody');
  let shown = showHiddenFiles ? entries : entries.filter(e => !e.name.startsWith('.'));
  const terms = searchTerms(_fileFilter);
  if (terms.length) shown = shown.filter(e => matchesAllTerms(e.name, terms));
  const vis = _sortEntries(shown);
  _visEntries = vis;
  // drop selected paths no longer listed (deleted, renamed, or hidden)
  _selected = new Set(vis.filter(e => _selected.has(e.path)).map(e => e.path));
  if (!vis.length) {
    tbody.innerHTML = '<tr><td colspan="4" class="td-empty">' +
      (terms.length ? 'No matches.' : 'Empty directory') + '</td></tr>';
    return;
  }
  const truncRow = truncated
    ? '<tr><td colspan="4" class="td-warn">' +
      iconSvg('warning') + ' Showing first 2000 entries — use the path bar to navigate deeper.</td></tr>'
    : '';
  tbody.innerHTML = truncRow + vis.map((e, i) => {
    const rowAttrs = ' data-idx="' + i + '"' + (_selected.has(e.path) ? ' class="selected"' : '') +
      ' onclick="onRowClick(event,' + i + ')" ondblclick="openEntry(' + i + ')"' +
      ' oncontextmenu="onRowContext(event,' + i + ')"';
    const nameCell = '<td><span class="file-name" onclick="event.stopPropagation();openEntry(' + i + ')"' +
      ' title="' + (e.is_dir ? 'Open folder' : 'Open file') + '">' +
      (e.is_dir ? iconSvg('folder', 'icon-folder', true) : iconSvg('draft')) + ' ' + escHtml(e.name) +
      '</span></td>';
    const dlTitle = e.is_dir ? 'Download as zip (symlinks inside the folder are skipped, not followed)' : 'Download';
    const dlBtn = '<button class="btn btn-sm btn-ghost btn-icon" onclick="rowDownload(event,' + i + ')" title="' + escHtml(dlTitle) + '">' + iconSvg('download') + '</button>';
    const moreBtn = '<button class="btn btn-sm btn-ghost btn-icon" onclick="toggleRowMenu(event,' + i + ')" title="More">' + iconSvg('more_vert') + '</button>';
    if (e.is_dir) {
      const isStarred = fileBookmarks.some(b => b.path === e.path);
      const bmBtn = '<button class="btn btn-sm btn-ghost btn-icon' + (isStarred ? ' starred' : '') + '" onclick="rowBookmark(event,' + i + ')" title="Bookmark">' +
        iconSvg('star', null, isStarred) + '</button>';
      return '<tr' + rowAttrs + '>' + nameCell +
        '<td class="td-muted">—</td>' +
        '<td class="td-muted">' + (e.modified_at ? fmtDate(e.modified_at) : '—') + '</td>' +
        '<td class="file-actions-cell">' + bmBtn + dlBtn + moreBtn + '</td>' +
        '</tr>';
    }
    return '<tr' + rowAttrs + '>' + nameCell +
      '<td class="td-muted">' + fmtSize(e.size) + '</td>' +
      '<td class="td-muted">' + (e.modified_at ? fmtDate(e.modified_at) : '—') + '</td>' +
      '<td class="file-actions-cell">' + dlBtn + moreBtn + '</td>' +
      '</tr>';
  }).join('');
  _renderCutMarks();
}

// viewLargeFileWarnBytes is the size above which viewPath asks for
// confirmation before opening a file.
const viewLargeFileWarnBytes = 1 << 30; // 1 GiB

// viewPath opens a file in a new tab with Content-Disposition: inline, so
// the browser renders it directly if it can (downloadPath, by contrast,
// always forces a save-as download).
async function viewPath(path, size) {
  if (size > viewLargeFileWarnBytes) {
    const ok = await askConfirm('This file is ' + fmtSize(size) + '. Open it anyway?',
      { title: 'Open large file?', okLabel: 'Open' });
    if (!ok) return;
  }
  window.open('/api/fs/download?path=' + encodeURIComponent(path) + '&inline=1', '_blank');
}

// downloadPath forces a save-as download via a temporary <a download>
// click, so several downloads kicked off together don't cancel each other.
function downloadPath(path) {
  const a = document.createElement('a');
  a.href = '/api/fs/download?path=' + encodeURIComponent(path);
  a.download = '';
  document.body.appendChild(a);
  a.click();
  a.remove();
}

/* ── Delete ──────────────────────────────── */
let _deleteConfirmPaths = null;

// deletePaths asks for confirmation once (label is the single entry's name
// or "N items"), then fires doDelete per path.
function deletePaths(paths, label) {
  _deleteConfirmPaths = paths;
  gid('delete-confirm-name').textContent = label;
  gid('delete-confirm-modal').classList.add('open');
}

function confirmDelete(mode) {
  gid('delete-confirm-modal').classList.remove('open');
  const paths = _deleteConfirmPaths;
  _deleteConfirmPaths = null;
  if (mode !== 'yes' || !paths) return;
  paths.forEach(p => doDelete(p));
}
// Clicking the backdrop cancels, matching the "cancel" button.
onBackdropClick(gid('delete-confirm-modal'), () => confirmDelete('cancel'));

// doDelete kicks off DELETE /api/fs as a background task and watches it
// via a bar progress pill.
async function doDelete(path) {
  const name = path.replace(/\/+$/, '').split('/').pop() || path;
  try {
    const r = await fetch('/api/fs?path=' + encodeURIComponent(path), { method: 'DELETE' });
    if (!r.ok) { showProgressError('Delete ' + name, await r.text()); return; }
    const { id } = await r.json();
    openBarProgress('Deleting ' + name, id, async msg => {
      if (!msg.ok) {
        // the pill itself shows the error / cancelled state
        navigateFiles(currentPath); // some entries may have been removed before it stopped
        return;
      }
      await _dropDanglingBookmark(path);
      navigateFiles(currentPath);
    }, { formatLabel: (cur, tot) => cur + ' / ' + tot + ' entries' });
  } catch (e) {
    showProgressError('Delete ' + name, String(e));
  }
}

// _dropDanglingBookmark removes bookmarks at path or inside its subtree —
// deleting/moving a folder strands those too. Best-effort: failures are
// swallowed.
async function _dropDanglingBookmark(path) {
  const prefix = path.replace(/\/+$/, '') + '/';
  const stale = fileBookmarks.filter(b => b.path === path || b.path.startsWith(prefix));
  for (const b of stale) {
    try {
      fileBookmarks = (await (await fetch('/api/fs/bookmarks?path=' + encodeURIComponent(b.path), { method: 'DELETE' })).json()) || [];
    } catch { /* best-effort */ }
  }
}

/* ── Rename ──────────────────────────────── */
let _renamePath = null;
let _renameOrig = null; // original name, so an unchanged submit is a no-op

function renamePath(path, name) {
  _renamePath = path;
  _renameOrig = name;
  const input = gid('rename-input');
  input.value = name;
  gid('rename-modal').classList.add('open');
  input.focus();
  // pre-select just the base name, not the extension
  const dot = name.lastIndexOf('.');
  if (dot > 0) input.setSelectionRange(0, dot);
  else input.select();
}

function confirmRename(mode) {
  gid('rename-modal').classList.remove('open');
  const path = _renamePath;
  _renamePath = null;
  if (mode !== 'yes' || !path) return;
  doRename(path);
}
onBackdropClick(gid('rename-modal'), () => confirmRename('cancel'));

gid('rename-input').addEventListener('keydown', e => {
  if (e.key === 'Enter') confirmRename('yes');
});

async function doRename(path) {
  const newName = gid('rename-input').value.trim();
  if (!newName || newName === _renameOrig) return;
  try {
    const r = await fetch('/api/fs?path=' + encodeURIComponent(path), {
      method: 'PATCH',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ new_name: newName })
    });
    if (!r.ok) { showProgressError('Rename ' + newName, await r.text()); return; }
  } catch (e) {
    showProgressError('Rename ' + newName, String(e));
    return;
  }
  await _dropDanglingBookmark(path);
  navigateFiles(currentPath);
}

/* ── Move ────────────────────────────────── */
// doMove PATCHes /api/fs with dest_dir (and new_name, when given, to move
// under a different name). A same-filesystem move responds synchronously
// with {"path": newPath}; a cross-filesystem move runs as a background
// task, responds with {"id": taskID}, and is watched via a bar progress
// pill.
async function doMove(path, name, destDir, newName) {
  let data;
  try {
    const body = { dest_dir: destDir };
    if (newName) body.new_name = newName;
    const r = await fetch('/api/fs?path=' + encodeURIComponent(path), {
      method: 'PATCH',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(body)
    });
    if (!r.ok) { showProgressError('Move ' + name, await r.text()); return; }
    data = await r.json();
  } catch (e) {
    showProgressError('Move ' + name, String(e));
    return;
  }

  if (data.id) {
    openBarProgress('Moving ' + name, data.id, async msg => {
      if (!msg.ok) return; // the pill itself shows the error / cancelled state
      await _dropDanglingBookmark(path);
      navigateFiles(currentPath);
    }, { formatLabel: (cur, tot) => cur + ' / ' + tot + ' entries' });
    return;
  }

  await _dropDanglingBookmark(path);
  navigateFiles(currentPath);
}

/* ── Copy ────────────────────────────────── */
// doCopy starts a background copy task (POST /api/fs/copy, with new_name
// when given to copy under a different name) and watches it via a bar
// progress pill.
async function doCopy(path, name, destDir, newName) {
  try {
    const body = { dest_dir: destDir };
    if (newName) body.new_name = newName;
    const r = await fetch('/api/fs/copy?path=' + encodeURIComponent(path), {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(body)
    });
    if (!r.ok) { showProgressError('Copy ' + name, await r.text()); return; }
    const { id } = await r.json();
    openBarProgress('Copying ' + name, id, msg => {
      if (!msg.ok) return; // the pill itself shows the error / cancelled state
      navigateFiles(currentPath);
    }, { formatLabel: (cur, tot) => cur + ' / ' + tot + ' entries' });
  } catch (e) {
    showProgressError('Copy ' + name, String(e));
  }
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
    if (!r.ok) { showProgressError('Bookmark', await r.text()); return; }
    fileBookmarks = (await r.json()) || [];
  } catch (e) {
    showProgressError('Bookmark', String(e));
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

/* ── Upload ──────────────────────────────── */
let _pendingUpload      = null; // [{file, relPath}] awaiting conflict resolution / upload
let _uploadConflictResolve = null; // resolves the promise inside resolveConflictsSequentially
let _uploadBusy  = false; // true from file selection through the end of the XHR (covers the conflict-resolution modal too)
let _uploadQueue = [];    // [{items, pill}] uploads picked/dropped while one was already in flight

// Only one upload runs at a time; a pick/drop while busy is queued with
// its own "Queued…" pill and starts automatically when the current upload
// finishes.
function _setUploadBusy(busy) {
  _uploadBusy = busy;
  if (!busy && _uploadQueue.length) {
    const next = _uploadQueue.shift();
    next.pill.el.remove();
    startUpload(next.items);
  }
}

function handleUploadInputChange(inputEl, isDir) {
  const files = Array.from(inputEl.files || []);
  inputEl.value = ''; // allow re-selecting the same file(s) later
  if (!files.length) return;
  startUpload(files.map(f => ({ file: f, relPath: f.webkitRelativePath || f.name })));
}

// startUpload runs the upload flow (conflict check → resolution → XHR) for
// items = [{file, relPath}]; queues if an upload is already in flight.
function startUpload(items) {
  if (!items.length) return;
  if (_uploadBusy) { _enqueueUpload(items); return; }
  _pendingUpload = items;
  _setUploadBusy(true);
  checkUploadConflicts();
}

function _enqueueUpload(items) {
  const label = items.length === 1 ? items[0].relPath : items.length + ' items';
  let entry;
  const pill = openBarProgress('Queued: ' + label, null, null, {
    onCancel: () => {
      _uploadQueue = _uploadQueue.filter(q => q !== entry);
      pill.el.remove();
    },
  });
  pill.el.querySelector('.prog-pill-label').textContent = 'Queued…';
  entry = { items, pill };
  _uploadQueue.push(entry);
}

/* ── Drag-and-drop upload ─────────────────── */
// _readAllDirEntries reads all entries out of a DirectoryReader, calling
// readEntries() until it returns an empty batch (Chrome caps each call at
// ~100 results).
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

// _walkDroppedEntry recursively walks a dropped FileSystemEntry, appending
// {file, relPath} pairs to out.
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

// dragover/dragenter must call preventDefault() to mark the element as a
// valid drop target — otherwise the browser handles the drop itself and
// navigates to the dropped file.
_dropZone.addEventListener('dragover', e => {
  if (!e.dataTransfer.types.includes('Files')) return;
  e.preventDefault();
});

_dropZone.addEventListener('dragenter', e => {
  if (!e.dataTransfer.types.includes('Files')) return;
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
    if (!r.ok) { showProgressError('Upload', await r.text()); _pendingUpload = null; _setUploadBusy(false); return; }
    const { conflicts } = await r.json();
    const skip = await resolveConflictsSequentially(conflicts || []);
    if (skip === null) { _pendingUpload = null; _setUploadBusy(false); return; } // cancelled
    doUpload(skip);
  } catch (e) {
    showProgressError('Upload', String(e));
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
onBackdropClick(gid('upload-conflict-modal'), () => resolveUploadConflict('cancel'));

function doUpload(skipSet) {
  const items = _pendingUpload.filter(p => !skipSet.has(p.relPath));
  _pendingUpload = null;
  if (!items.length) { _setUploadBusy(false); return; }
  // relPath goes in the field *name*, not the filename: the server relies on
  // this because Go's mime/multipart strips any "/" from Part.FileName() but
  // leaves Part.FormName() untouched (see api_files_upload.go).
  const fd = new FormData();
  items.forEach(p => fd.append(p.relPath, p.file, p.file.name));

  const title = items.length === 1 ? 'Uploading ' + items[0].relPath : 'Uploading ' + items.length + ' items';
  let xhr; // declared before openBarProgress so the onCancel closure can reach it
  // No taskId: progress comes from the XHR's upload events, so the pill is
  // driven manually via bar.update()/bar.done(); onCancel aborts the XHR.
  const bar = openBarProgress(title, null, null, {
    formatLabel: (cur, tot) => fmtSize(cur) + ' / ' + fmtSize(tot),
    onCancel: () => xhr && xhr.abort(),
  });

  xhr = new XMLHttpRequest();
  xhr.open('POST', '/api/fs/upload?dest=' + encodeURIComponent(currentPath) + '&overwrite=true');
  xhr.upload.onprogress = e => {
    if (e.lengthComputable) bar.update(e.loaded, e.total);
  };
  xhr.onload = () => {
    _setUploadBusy(false);
    if (xhr.status >= 200 && xhr.status < 300) {
      bar.done(true);
      navigateFiles(currentPath);
    } else {
      bar.done(false, xhr.responseText || 'Upload failed');
    }
  };
  xhr.onerror = () => { _setUploadBusy(false); bar.done(false, 'Network error'); };
  xhr.onabort = () => { _setUploadBusy(false); bar.done(false, 'Cancelled'); navigateFiles(currentPath); };
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
  _fpOnSelect = null; // this is form-fill mode, not callback mode
  _fpSetNameRow(null);
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
  _renderFpBookmarkSelect();
  loadFileBookmarks().then(_renderFpBookmarkSelect); // picks up bookmarks even if the Files tab was never visited this session
  fpNavigate(fpPath);
  gid('fp-modal').classList.add('open');
}

// openFilePickerForCallback opens the directory-picker modal and invokes
// onSelect(path, newName) with the chosen directory, instead of filling a
// form field like openFilePicker. Always dir mode. opts.name prefills an
// editable Name field (rename on move/copy); newName is undefined when the
// field is hidden or left empty.
function openFilePickerForCallback(initialPath, title, onSelect, opts) {
  fpTargetId  = '';
  fpMode      = 'dir';
  fpSuffix    = '';
  _fpOnSelect = onSelect;
  _fpSetNameRow(opts && opts.name);
  fpPath = initialPath || srvScratch || srvHome || '/';
  gid('fp-title').textContent = title || 'Select Directory';
  const selDirBtn = gid('fp-select-dir-btn');
  if (selDirBtn) selDirBtn.style.display = '';
  _renderFpBookmarkSelect();
  loadFileBookmarks().then(_renderFpBookmarkSelect);
  fpNavigate(fpPath);
  gid('fp-modal').classList.add('open');
}

// _renderFpBookmarkSelect populates the file picker's "Go to…" dropdown
// with the permanent shortcuts (Root/$HOME/$SCRATCH) and the bookmarks.
function _renderFpBookmarkSelect() {
  const sel = gid('fp-bookmark-select');
  if (!sel) return;
  const permanent = [['/', '/ Root']];
  if (srvHome)    permanent.push([srvHome,    '$HOME']);
  if (srvScratch) permanent.push([srvScratch, '$SCRATCH']);
  let html = '<option value="">Go to…</option>' +
    '<optgroup label="Shortcuts">' +
      permanent.map(([p, label]) => '<option value="' + escHtml(p) + '">' + escHtml(label) + '</option>').join('') +
    '</optgroup>';
  if (fileBookmarks.length) {
    html += '<optgroup label="Bookmarks">' +
      fileBookmarks.map(b => '<option value="' + escHtml(b.path) + '">' + escHtml(b.label || b.path) + '</option>').join('') +
      '</optgroup>';
  }
  sel.innerHTML = html;
}

function fpJumpToBookmark(selectEl) {
  const path = selectEl.value;
  selectEl.value = '';
  if (path) fpNavigate(path);
}

function _renderFpBreadcrumb(path) { _renderBreadcrumb('fp-bc', path, 'fpNavigate'); }

let _fpController = null; // AbortController for the in-flight fpNavigate request

async function fpNavigate(path) {
  if (_fpController) _fpController.abort();
  _fpController = new AbortController();
  fpPath = path;
  const url = path ? '/api/fs?path=' + encodeURIComponent(path) : '/api/fs';
  try {
    const r = await fetch(url, { signal: _fpController.signal });
    if (!r.ok) {
      const msg = (await r.text()).trim();
      gid('fp-body').innerHTML =
        '<div class="modal-error">' + escHtml(msg || 'Error loading directory') + '</div>';
      return;
    }
    const entries = (await r.json()) || [];
    // Resolve actual path from API response
    if (entries.length > 0) {
      fpPath = entries[0].path.replace(/\/[^/]+$/, '') || '/';
    }
    _renderFpBreadcrumb(fpPath);
    _setVal('fp-path-input', fpPath);
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

// _fpSetNameRow shows the picker's Name field prefilled with name, or
// hides it when name is falsy.
function _fpSetNameRow(name) {
  gid('fp-name-row').hidden = !name;
  gid('fp-name-input').value = name || '';
}

gid('fp-name-input').addEventListener('keydown', e => {
  if (e.key === 'Enter') fpSelect();
});

function fpSelectDir(path) {
  if (_fpOnSelect) {
    const cb = _fpOnSelect;
    _fpOnSelect = null;
    const newName = gid('fp-name-row').hidden ? '' : gid('fp-name-input').value.trim();
    closeModal('fp-modal');
    cb(path, newName || undefined);
    return;
  }
  _setVal(fpTargetId, path);
  const e = gid(fpTargetId);
  if (e) e.dispatchEvent(new Event('input'));
  closeModal('fp-modal');
}
function fpSelect() { fpSelectDir(fpPath); }
