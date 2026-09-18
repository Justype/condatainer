/* ── Project section ─────────────────────── */
// No cwd field of its own: this panel follows the same shared cfg-cwd input
// every other cwd-dependent panel already defers to (helpers.js's launch
// payload, overlays.js's file-picker default) — one "which directory"
// concept for the user to manage, not a second one just for this panel.
let _projectStatus = null;

async function loadProjectStatus() {
  const panel = gid('project-status-panel');
  if (panel && !_projectStatus) panel.innerHTML = '<div class="modal-empty">Loading…</div>';
  const cwd = gid('cfg-cwd')?.value || '';
  try {
    const r = await fetch('/api/project/status?cwd=' + encodeURIComponent(cwd));
    _projectStatus = (await r.json()) || null;
  } catch {
    _projectStatus = null;
  }
  renderProjectStatus();
}

function renderProjectStatus() {
  const panel = gid('project-status-panel');
  const lockBtn = gid('project-lock-btn');
  if (!panel) return;
  if (!_projectStatus) {
    panel.innerHTML = '<div class="modal-error">Could not load project status.</div>';
    return;
  }
  if (!_projectStatus.has_project) {
    panel.innerHTML = '<div class="modal-empty">No project here — no cnt-lock/ found for the working directory set on Start.</div>';
    if (lockBtn) lockBtn.innerHTML = iconSvg('folder') + ' Create project (lock)';
    return;
  }
  if (lockBtn) lockBtn.innerHTML = iconSvg('folder') + ' Lock project';

  let html = '<div class="param-grid">' +
    '<span class="param-name">Root</span><div class="param-cell mono">' + escHtml(_projectStatus.root) + '</div>' +
    '<span class="param-name">Pins</span><div class="param-cell">' + (_projectStatus.pin_count || 0) + '</div>' +
    '</div>';

  const unpinned = _projectStatus.unpinned || [];
  if (unpinned.length) {
    html += '<div class="f-divider">Declared but not pinned (' + unpinned.length + ')</div>' +
      '<div class="oc-note">' + iconSvg('warning') +
      '<span>Run <code>condatainer project lock</code> or <code>project pin</code> to pin these.</span></div>' +
      '<div class="mono" style="font-size:12px;">' + unpinned.map(escHtml).join('<br>') + '</div>';
  }

  const usedNotPinned = _projectStatus.unpinned_helper_overlays || [];
  if (usedNotPinned.length) {
    html += '<div class="f-divider">Overlays used by helpers but not pinned</div>' +
      '<div class="mono" style="font-size:12px;">' + usedNotPinned.map(escHtml).join('<br>') + '</div>';
  }

  const usage = _projectStatus.manual_pin_usage || {};
  const usageKeys = Object.keys(usage).sort();
  if (usageKeys.length) {
    html += '<div class="f-divider">Manual pin usage</div>' +
      usageKeys.map(k => {
        const names = usage[k] || [];
        const label = names.length
          ? 'used by: ' + names.map(escHtml).join(', ')
          : 'no recorded helper usage';
        return '<div class="mono" style="font-size:12px;">' + escHtml(k) +
          ' <span class="td-muted">— ' + label + '</span></div>';
      }).join('');
  }

  panel.innerHTML = html;
}

async function lockProject() {
  const btn = gid('project-lock-btn');
  if (btn) btn.disabled = true;
  const cwd = gid('cfg-cwd')?.value || '';
  try {
    const r = await fetch('/api/project/lock', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ cwd }),
    });
    if (!r.ok) {
      showProgressError('Locking project', await r.text());
      if (btn) btn.disabled = false;
      return;
    }
    const { id } = await r.json();
    openLogProgress('Locking project', id, async () => {
      if (btn) btn.disabled = false;
      await loadProjectStatus();
    });
  } catch (e) {
    if (btn) btn.disabled = false;
    showProgressError('Locking project', String(e));
  }
}

// Re-reads on the same input event cfg-cwd's other consumers already listen
// for, so switching directories on Start updates this panel too, without a
// manual refresh — but only while the panel is visible, to avoid a fetch on
// every keystroke elsewhere.
gid('cfg-cwd')?.addEventListener('input', () => {
  if (gid('sec-project')?.classList.contains('active')) loadProjectStatus();
});
