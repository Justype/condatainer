/* ── Shared utilities ─────────────────────── */
function escHtml(s) {
  return String(s == null ? '' : s)
    .replace(/&/g, '&amp;').replace(/</g, '&lt;')
    .replace(/>/g, '&gt;').replace(/"/g, '&quot;');
}
function fmtSize(bytes) {
  if (!bytes) return '—';
  const u = ['B', 'KB', 'MB', 'GB', 'TB'];
  let i = 0, b = bytes;
  while (b >= 1024 && i < u.length - 1) { b /= 1024; i++; }
  return b.toFixed(1) + ' ' + u[i];
}
function fmtDate(iso) {
  if (!iso) return '—';
  try { return new Date(iso).toLocaleString(); } catch { return iso; }
}
function fmtRelTime(iso) {
  if (!iso) return '—';
  const diff = Date.now() - new Date(iso).getTime();
  const m = Math.floor(diff / 60000);
  if (m < 1) return 'just now';
  if (m < 60) return m + 'm ago';
  const h = Math.floor(m / 60);
  if (h < 24) return h + 'h ago';
  return Math.floor(h / 24) + 'd ago';
}
function pathTailHtml(path, maxWidth, cls) {
  const s = String(path || '');
  const style = maxWidth ? ' style="max-width:' + escHtml(maxWidth) + '"' : '';
  const classes = 'path-tail' + (cls ? ' ' + cls : '');
  return '<span class="' + classes + '"' + style + ' title="' + escHtml(s) + '"><bdi>' + escHtml(s || '—') + '</bdi></span>';
}
function calcDuration(startIso, endIso, status, walltimeStr) {
  if (!startIso) return '—';
  if (!endIso && status === 'done') return walltimeStr || '—';
  if (!endIso && !['pending', 'starting', 'running'].includes(status)) return '—';
  const s = new Date(startIso).getTime();
  const e = endIso ? new Date(endIso).getTime() : Date.now();
  if (!Number.isFinite(s) || !Number.isFinite(e) || e < s) return '—';
  const diff = e - s;
  const h = Math.floor(diff / 3600000);
  const m = Math.floor((diff % 3600000) / 60000);
  return h ? h + 'h ' + m + 'm' : m + 'm';
}
function badgeClass(status) {
  if (status === 'running') return 'badge-green';
  if (status === 'pending' || status === 'starting') return 'badge-amber';
  if (status === 'failed')  return 'badge-red';
  return 'badge-grey';
}
function gid(id) { return document.getElementById(id); }
function _setText(id, val) { const e = gid(id); if (e) e.textContent = val || '—'; }
function _setVal(id, val)  { const e = gid(id); if (e) e.value = val; }
// Sets cfg-cwd value and placeholder. Placeholder is always $SCRATCH (or $HOME).
// Value is set to val if non-empty, otherwise defaults to $SCRATCH/$HOME.
function _setCwd(val) {
  const e = gid('cfg-cwd');
  if (!e) return;
  const def = srvScratch || srvHome || '';
  e.placeholder = def;
  e.value = val || def;
}
// fill: true → i-fill-{name} (filled variant); default → i-{name} (outline)
function iconSvg(name, cls, fill) {
  const id = fill ? 'fill-' + name : name;
  return '<svg class="icon' + (cls ? ' ' + cls : '') + '" aria-hidden="true"><use href="#i-' + id + '"></use></svg>';
}

/* ── Global state ────────────────────────── */
let allJobs     = [];
let allHelpers  = [];
let allOverlays = [];
let allHistory  = [];

// start section
let selectedHelper   = null;
let selectedModules          = [];
let selectedExternalOverlays = [];
let _helperParamKeys = [];

// detail panel
let detailJobId     = null;
let detailTabActive = 'messages';
let detailSSE       = null;

// file browser
let currentPath = '';
let srvHome    = '';
let srvScratch = '';
let srvNotification = '';

// file picker modal
let fpTargetId = '', fpMode = 'dir', fpPath = '', fpSuffix = '';

// overlay picker modal
let opTargetType  = 'module';

/* ── Theme ───────────────────────────────── */
const _sysMq      = matchMedia('(prefers-color-scheme: dark)');
const _savedTheme = localStorage.getItem('conda-theme');
let theme = _savedTheme || 'system';

function _resolvedTheme() {
  return theme === 'system' ? (_sysMq.matches ? 'dark' : 'light') : theme;
}
function _applyTheme() {
  document.documentElement.setAttribute('data-theme', _resolvedTheme());
}
_applyTheme();
_updateThemeBtn();

_sysMq.addEventListener('change', () => { if (theme === 'system') { _applyTheme(); _updateThemeBtn(); } });

function toggleTheme() {
  theme = theme === 'light' ? 'dark' : theme === 'dark' ? 'system' : 'light';
  localStorage.setItem('conda-theme', theme);
  _applyTheme();
  _updateThemeBtn();
}
function _updateThemeBtn() {
  const icon  = gid('theme-icon');
  const label = document.querySelector('#theme-toggle .sb-btn-label');
  if (icon)  icon.innerHTML  = iconSvg(theme === 'dark' ? 'dark_mode' : theme === 'light' ? 'light_mode' : 'computer');
  if (label) label.textContent = theme === 'dark' ? 'Dark Mode' : theme === 'light' ? 'Light Mode' : ('System (' + (_sysMq.matches ? 'Dark' : 'Light') + ')');
}
gid('theme-toggle').addEventListener('click', toggleTheme);

/* ── Sidebar collapse ───────────────────── */
let sidebarCollapsed = localStorage.getItem('sidebarCollapsed') === 'true';
const collapseTab = gid('collapse-tab');
function _updateCollapseTab() {
  collapseTab.innerHTML = iconSvg(sidebarCollapsed ? 'chevron_right' : 'chevron_left');
  collapseTab.title = sidebarCollapsed ? 'Expand sidebar' : 'Collapse sidebar';
}
collapseTab.addEventListener('click', () => {
  sidebarCollapsed = !sidebarCollapsed;
  localStorage.setItem('sidebarCollapsed', sidebarCollapsed);
  gid('sidebar').classList.toggle('collapsed', sidebarCollapsed);
  _updateCollapseTab();
});
const sidebar = gid('sidebar');
sidebar.classList.add('no-transition');
sidebar.classList.toggle('collapsed', sidebarCollapsed);
requestAnimationFrame(() => requestAnimationFrame(() => sidebar.classList.remove('no-transition')));
_updateCollapseTab();

/* ── Navigation ──────────────────────────── */
document.querySelectorAll('.nav-item').forEach(el => {
  el.addEventListener('click', () => navigate(el.dataset.section));
});
function navigate(id) {
  document.querySelectorAll('.nav-item').forEach(e =>
    e.classList.toggle('active', e.dataset.section === id));
  document.querySelectorAll('.section').forEach(s =>
    s.classList.toggle('active', s.id === 'sec-' + id));
  document.querySelectorAll('.sb-btn[data-section]').forEach(b =>
    b.classList.toggle('active-section', b.dataset.section === id));
  closeDetail();
  if      (id === 'jobs')     loadJobs();
  else if (id === 'start')    loadHelpers();
  else if (id === 'overlays') loadOverlays();
  else if (id === 'history')  loadHistory();
  else if (id === 'files')    { renderFileTree(); navigateFiles(currentPath); }
}

/* ── Status polling ──────────────────────── */
async function pollStatus() {
  try {
    const r = await fetch('/api/status', { signal: AbortSignal.timeout(2000) });
    const d = await r.json();
    _setStatus(true, d);
  } catch {
    _setStatus(false, null);
  }
}
function _setStatus(alive, d) {
  const dot  = gid('status-dot');
  const text = gid('status-text');
  if (dot)  dot.classList.toggle('dead', !alive);
  if (text) text.textContent = alive ? (d.running || 0) + ' running' : 'unreachable';
  if (d) {
    _setText('srv-port',    d.port);
    _setText('srv-pid',     d.pid);
    _setText('srv-uptime',  d.uptime);
    _setText('srv-version', d.version);
    _setText('ssh-hint', 'LocalForward ' + d.port + ' localhost:' + d.port);
    if (d.home)    srvHome    = d.home;
    if (d.scratch) srvScratch = d.scratch;
    const cwdEl = gid('cfg-cwd');
    if (cwdEl) {
      const def = srvScratch || srvHome || '';
      cwdEl.placeholder = def;
      if (!cwdEl.value) cwdEl.value = def;
    }
    const notif = d.notification || '';
    srvNotification = notif;
    _updateNotifCard();
    renderFileTree();
  }
}

/* ── Browser notifications ───────────────── */
function _updateNotifCard() {
  const card = gid('notif-card');
  const statusEl = gid('notif-status');
  const btn = gid('notif-btn');
  if (!card) return;
  const webEnabled = srvNotification === 'web' || srvNotification === 'both';
  card.style.display = webEnabled ? '' : 'none';
  if (!webEnabled) return;
  if (typeof Notification === 'undefined') {
    statusEl.textContent = 'Browser notifications are not supported.';
    btn.style.display = 'none';
    return;
  }
  if (Notification.permission === 'granted') {
    statusEl.textContent = 'Browser notifications are enabled.';
    btn.style.display = 'none';
  } else if (Notification.permission === 'denied') {
    statusEl.textContent = 'Notifications blocked. Allow them in your browser site settings.';
    btn.style.display = 'none';
  } else {
    statusEl.textContent = 'Click below to receive a notification when a helper job starts.';
    btn.style.display = '';
  }
}

function requestNotifPermission() {
  if (typeof Notification === 'undefined') return;
  Notification.requestPermission().then(() => _updateNotifCard());
}

/* ── Modal helpers ───────────────────────── */
function closeModal(id) { gid(id).classList.remove('open'); }
document.querySelectorAll('.modal-backdrop').forEach(el => {
  el.addEventListener('click', e => { if (e.target === el) closeModal(el.id); });
});

/* ── Keyboard shortcuts ──────────────────── */
document.addEventListener('keydown', e => {
  if (e.key === 'Escape') {
    closeDetail();
    closeModal('fp-modal');
    closeModal('op-modal');
  }
});

/* ── Outside-click closes detail ────────── */
gid('app').addEventListener('click', e => {
  const detail = gid('detail');
  if (detail.classList.contains('open')
      && !detail.contains(e.target)
      && !e.target.closest('.job-card')
      && !e.target.closest('tr[data-id]')) {
    closeDetail();
  }
});

/* ── Copy helpers ────────────────────────── */
function initCodeBlocks() {
  document.querySelectorAll('.code-block').forEach(el => {
    if (el.querySelector('button')) return;
    const btn = document.createElement('button');
    btn.className = 'btn btn-sm';
    btn.innerHTML = iconSvg('content_copy');
    btn.onclick = function() { copyStr(el.querySelector(':not(button)')?.textContent || '', this); };
    el.appendChild(btn);
  });
}
function copyStr(str, btn) {
  navigator.clipboard.writeText(str).then(() => {
    const o = btn.innerHTML;
    btn.innerHTML = iconSvg('check');
    setTimeout(() => { btn.innerHTML = o; }, 1500);
  });
}
initCodeBlocks();
